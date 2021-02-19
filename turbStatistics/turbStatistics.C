#include "turbStatistics.H"

#include "fvCFD.H"
#include "fvm.H"
#include "fvcGrad.H"
#include "tensorField.H"
#include "vectorField.H"
#include "dimensionedTypes.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace functionObjects {
    // * * * * * * * * * * * * * * * * Runtime selection   * * * * * * * * * //
    defineTypeNameAndDebug(turbStatistics, 0);

    addToRunTimeSelectionTable
    (
       functionObject,
       turbStatistics,
       dictionary
    );

    // * * * * * * * * * * * * * * * * Global properites and utilities   * * * * * * * * * //

    const dimensionSet dimStat( dimLength*dimLength/pow3(dimTime) );
    const dimensionSet dimIncPressure( dimVelocity*dimVelocity );

    const word DIV_SCHEME_NAME = "statistics";


    volScalarField &turbStatistics::newScalarField(const word &name, const fvMesh &mesh, dimensionSet dim, IOobject::readOption read, IOobject::writeOption write) {
        dimensioned<scalar> dimValue("value", dim, 0.);
        return mesh.objectRegistry::store(new volScalarField
                                                (
                                                    IOobject
                                                    (
                                                        name,
                                                        mesh.time().timeName(),
                                                        mesh,
                                                        read,
                                                        write
                                                    ),
                                                    mesh,
                                                    dimValue
                                                )
                                          );
    }

    volVectorField &turbStatistics::newVectorField(const word &name, const fvMesh &mesh, dimensionSet dim, IOobject::readOption read, IOobject::writeOption write) {
        dimensioned<vector> dimValue("value", dim, vector::zero);
        return mesh.objectRegistry::store(new volVectorField
                                                (
                                                    IOobject
                                                    (
                                                        name,
                                                        mesh.time().timeName(),
                                                        mesh,
                                                        read,
                                                        write
                                                    ),
                                                    mesh,
                                                    dimValue
                                                )
                                          );
    }

    volTensorField &turbStatistics::newTensorField(const word &name, const fvMesh &mesh, dimensionSet dim, IOobject::readOption read, IOobject::writeOption write) {
        dimensioned<tensor> dimValue("value", dim, tensor::zero);
        return mesh.objectRegistry::store(new volTensorField
                                            (
                                                IOobject
                                                (
                                                    name,
                                                    mesh.time().timeName(),
                                                    mesh,
                                                    read,
                                                    write
                                                ),
                                                mesh,
                                                dimValue
                                            )
                                          );
    }

    const scalarField &turbStatistics::density() const {
        return mesh_.lookupObject<scalarField>(rhoFieldName);
    }


    const dimensionSet& findPressureDim(const fvMesh& mesh) {
        const volScalarField& p = mesh.lookupObjectRef<volScalarField>("p");
        return p.dimensions();
    }

    dimensionSet findRhoDim(const fvMesh& mesh) {
        const volScalarField& p = mesh.lookupObjectRef<volScalarField>("p");
        return p.dimensions() == dimIncPressure ? dimless : dimMass / pow3(dimLength);
    }

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
    turbStatistics::turbStatistics
    (
        const word& name,
        const Time& time,
        const dictionary& dict
    )
    :
        fvMeshFunctionObject(name, time, dict),
        name_(name),
        obr_(mesh_),
        timeStep_(0),
        nu_("viscosity",dimViscosity,1),
        rhoRef_(1.0),
        meanU_( newVectorField("meanU", mesh_, dimVelocity, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE) ),
        up_( newVectorField("Up", mesh_, dimVelocity)  ),
        meanP_( newScalarField("meanP", mesh_, findPressureDim(mesh_), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE) ),
        pp_( newScalarField("pp", mesh_, findPressureDim(mesh_)) ),
        meanUU_( newTensorField("meanUU", mesh_, dimVelocity*dimVelocity, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE) ),
        tau_( newTensorField("tau", mesh_, dimIncPressure) ),
        k_( newScalarField("turbKinEnergy", mesh_, dimVelocity*dimVelocity, IOobject::NO_READ, IOobject::AUTO_WRITE) ),
        convK_( newScalarField("convK", mesh_, dimStat) ),
        prodK_( newScalarField("prodK", mesh_, dimStat) ),
        meanGradU_ddot_GradU_(newScalarField("meanGradUGradU", mesh_, dimless/(dimTime*dimTime), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE)  ),
        eps_( newScalarField("eps", mesh_, dimStat) ),
        diffK_( newScalarField("diffK", mesh_, dimStat) ),
        meanUUU_( newVectorField("menUUU", mesh_, pow3(dimVelocity), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE) ),
        triPointCorrel_( newScalarField("triPointCorellation", mesh_, dimStat) ),
        meanPU_( newVectorField("meanPU", mesh_, findPressureDim(mesh_)*dimVelocity, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE) ),
        pUFlucts_( newScalarField("pUFlucts", mesh_, dimStat) ),
        tauW_( newTensorField("tauW", mesh_, dimVelocity*dimVelocity) ),
        statisticTime_
        (
            IOobject
            (
                "statisticTime",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            )
        )
    {
        read(dict);
        timeStep_ = statisticTime_.lookupOrDefault<label>("statisticTimeStep", 0);
    }

    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    turbStatistics::~turbStatistics()
    {}


    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    bool turbStatistics::read(const dictionary& dict) {

//        const IOdictionary &trProp = obr_.lookupObject<IOdictionary>("transportProperties");
//        trProp.lookup("nu") >> nu_;
        dict.lookup("nu") >> nu_;

        if(dict.found("rhoRef")) {
            dict.lookup("rhoRef") >> rhoRef_;
        }

        rhoFieldName = dict.lookupOrDefault<word>("rhoName", "rho");

        return true;
    }

    bool turbStatistics::execute() {
        ++timeStep_;

        statisticTime_.add("statisticTimeStep",timeStep_,true);

        updateMeanAndFluctFields();
        updateStressTensor();
        updateKineticEnergy();
        updateConvectionStatistics();
        updateProductionStatistics();
        updateDyssypationStatisticsAndLaminarStress();
        updateDiffusiveStatistics();
        updateTriPointsCorellationStatistics();
        updatePressureVelocityFluctStatistics();

        return true;
    }

    bool turbStatistics::write() {
        return true;
    }


    // * * * * * * * * * * * * * * * Internal Fields Handling  * * * * * * * * * //

    void turbStatistics::updateMeanAndFluctFields() {
        const volVectorField & U = obr_.lookupObject<volVectorField>("U");
        const volScalarField & p = obr_.lookupObject<volScalarField>("p");

        meanU_ = meanU_*oldFrac() + U*newFrac();
        meanP_ = meanP_*oldFrac() + p*newFrac();

        up_ = U - meanU_;
        pp_ = p - meanP_;

        meanUU_= meanUU_*oldFrac() + ( U*U )*newFrac() ;
    }

    void turbStatistics::updateStressTensor() {
        tau_ = -(meanUU_ - meanU_*meanU_);
//        if(tau_.dimensions() == dimIncPressure)
//            tau_ = -(meanUU_ - meanU_*meanU_);
//        else {
//            tau_ = - density() * (meanUU_ - meanU_*meanU_);
//        }
    }

    void turbStatistics::updateKineticEnergy() {
         k_ =  - 0.5 * tr(tau_);
    }

    void turbStatistics::updateConvectionStatistics() {
        //convK_ = fvc::div(linearInterpolate(meanU_) & mesh.Sf(), k_,"div(phi,k)");
        convK_ = fvc::div(meanU_ * k_, DIV_SCHEME_NAME);//, "div(phi,k)"
    }

    void turbStatistics::updateProductionStatistics() {
        prodK_ = tau_ && fvc::grad(meanU_);
    }

    // Warning: here is calculated grad of average -
    // but it was stated that we should calculate average in
    // time from each time gradient, which in my opinion is
    // the same, and this approach do not require additionl
    // temporal field of velocity gradient
    void turbStatistics::updateDyssypationStatisticsAndLaminarStress() {
        volTensorField gradMeanU(fvc::grad(meanU_));
        const volVectorField & U = obr_.lookupObject<volVectorField>("U");
        volTensorField gradU(fvc::grad(U));

        meanGradU_ddot_GradU_ = meanGradU_ddot_GradU_ *oldFrac() + ( gradU && gradU )*newFrac();

        eps_ =  - nu_ * ( meanGradU_ddot_GradU_ - ( gradMeanU && gradMeanU ) );

        tauW_ = nu_ * gradMeanU;

        //eps_ = (eps_*(timeStep_ - 1) - nu_* ( (gradMeanU - dimensionedTensor("one",dimless/dimTime, tensor::one) ) && gradMeanU ) )/ timeStep_;
    }

    void turbStatistics::updateDiffusiveStatistics() {
        diffK_ = fvc::laplacian(nu_, k_);
    }

//    old approach
//    void turbStatistics::updateTriPointsCorellationStatistics()
//    {
//        const volVectorField & U = obr_.lookupObject<volVectorField>("U");

//        meanUUU_ = ( meanUUU_*(timeStep_-1) + (U & U)*U )/timeStep_;

//       //tmp<volVectorField> uuu = meanUUU_ - (meanU_ & meanU_)*meanU_ + 2.*(meanU_*tau_); //....-uj(uiui)
//        volVectorField uuu(meanUUU_);

//        forAll(uuu, cellI)
//        {
//            for(label j=0; j<3; ++j)
//            {
//                for(label i=0; i<3; ++i)
//                {
//                    scalar mUi = meanU_.internalField()[cellI][i];
//                    scalar mUj = meanU_.internalField()[cellI][j];
//                    scalar mUiUj = - tau_.internalField()[cellI][3*i+j];
//                    scalar mUiUi = - tau_.internalField()[cellI][3*i+i];
//                    uuu.internalField()[cellI][j]+= -mUi*mUi*mUj - 2*mUi*mUiUj - mUj*mUiUi;
//                }
//            }
//        }

//        triPointCorrel_ = -0.5*fvc::div(uuu,"div(phi,k)");
//    }

    void turbStatistics::updateTriPointsCorellationStatistics() {
        const volVectorField & U = obr_.lookupObject<volVectorField>("U");

        meanUUU_ = meanUUU_*oldFrac() + ( (U & U)*U )*newFrac();

        triPointCorrel_ = -0.5*fvc::div( meanUUU_ - (meanU_ & meanU_)*meanU_  + 2*(meanU_ & tau_) -2*meanU_*k_, DIV_SCHEME_NAME);//, "div(phi,k)");
    }

    void turbStatistics::updatePressureVelocityFluctStatistics() {
        const volVectorField & U = obr_.lookupObject<volVectorField>("U");
        const volScalarField & p = obr_.lookupObject<volScalarField>("p");

        meanPU_ = meanPU_*oldFrac() + (p*U)*newFrac();

        if(p.dimensions() == dimIncPressure) {
            pUFlucts_ = -fvc::div( meanPU_ - meanP_*meanU_, DIV_SCHEME_NAME);//, "div(phi,k)"); //1./rhoRef_*
        }
        else {
            const volScalarField rho = lookupObjectRef<volScalarField>(rhoFieldName);
            pUFlucts_ = -fvc::div( (meanPU_ - meanP_*meanU_) / rho, DIV_SCHEME_NAME);//, "div(phi,k)");
        }
    }


}
// ************************************************************************* //
}
// ************************************************************************* //
