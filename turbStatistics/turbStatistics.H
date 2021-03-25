#ifndef TURB_STATISTICS_H
#define TURB_STATISTICS_H

#include "word.H"
#include "volMesh.H"
#include "volFieldsFwd.H"
#include "autoPtr.H"
//#include "OutputFilterFunctionObject.H"
#include "dimensionedTypes.H"
#include "Time.H"
#include "fvMeshFunctionObject.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes
class objectRegistry;
class dictionary;
class word;
class mapPolyMesh;
class polyMesh;
class fvMesh;

namespace functionObjects {
/*---------------------------------------------------------------------------*\
                         Class velGrad Declaration
\*---------------------------------------------------------------------------*/

extern const dimensionSet dimStat;
extern const dimensionSet dimIncPressure;

class turbStatistics : public fvMeshFunctionObject
{
protected:

    // Private data

        //- Name of this calculator
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        label timeStep_;

        dimensionedScalar nu_;

        scalar rhoRef_;
        word rhoFieldName;
        //---------------------------------------------------
        // u - velocity field
        // p - pressure field
        // up - velocity fluctuation field
        // pp - pressure fluctuation field
        // mean(..) - time averaging

        // - mean u
        volVectorField& meanU_;
        // - mean u fluctuation
        volVectorField& meanUp_;
        // - mean voritcity
        volVectorField& meanVorticity_;

        // - mean p
        volScalarField& meanP_;

        // - p fluctuation
//        volScalarField& pp_;

        // - meanUU - real velocity multiplication
        volTensorField& meanUU_;
        // - reynolds stress tensor: mean(up_i*up_j)
        volTensorField& tau_;
        // - turbulent kinetic energy: 1/2*(sum(-tau_ii))
        volScalarField& k_;
        // - convection statistic : div(mean(u_j)*k)
        volScalarField& convK_;
        // - production statisctics: tau:grad(mean(u))
        volScalarField& prodK_;
        // temporal mean( sum_ij( (grad(U))^2 ) )
        volScalarField& meanGradU_ddot_GradU_;
        // - dysypation statistics: -nu*sum_ij(mean( grad(up) outer grad(up) ))
        volScalarField& eps_;
        // - diffusion statistics: -nu * laplacian(k) = -nu * div(grad(k))
        volScalarField& diffK_;
        // - mean 3 u multiplication vector:meanUUU[j]=mean(sum_i(u_i*u_i)*u_j)
        volVectorField& meanUUU_;
        // - 3 points correlation statistics: -1/2 * div(sum_i( mean(up_i*up_i*up_j) ) )
        volScalarField& triPointCorrel_;
        //mean (pressure mul vel)
        volVectorField& meanPU_;
        // - pressure-velocity fluctuations: -1/rho*div(mean(pp*up))
        volScalarField& pUFlucts_;
        // laminar tress nu*grad(mean(U))
        volTensorField& tauW_;
        //------------------------------------------------------

        IOdictionary statisticTime_;



    // Private Member Functions

        //- Disallow default bitwise copy construct
        turbStatistics(const turbStatistics&);

        //- Disallow default bitwise assignment
        void operator=(const turbStatistics&);

        inline scalar oldFrac() const;
        inline scalar newFrac() const;


        void updateMeanFields();
        void updateStressTensor();
        void updateKineticEnergy();
        void updateConvectionStatistics();
        void updateProductionStatistics();
        void updateDyssypationStatisticsAndLaminarStress();
        void updateDiffusiveStatistics();
        void updateTriPointsCorellationStatistics();
        void updatePressureVelocityFluctStatistics();


        static volScalarField& newScalarField(const word &name, const fvMesh& mesh, dimensionSet dim, IOobject::readOption read=IOobject::NO_READ, IOobject::writeOption write=IOobject::NO_WRITE);
        static volVectorField& newVectorField(const word &name, const fvMesh& mesh, dimensionSet dim, IOobject::readOption read=IOobject::NO_READ, IOobject::writeOption write=IOobject::NO_WRITE);
        static volTensorField& newTensorField(const word &name, const fvMesh& mesh, dimensionSet dim, IOobject::readOption read=IOobject::NO_READ, IOobject::writeOption write=IOobject::NO_WRITE);


        const scalarField& density() const;

public:

    //- Runtime type information
    TypeName("turbStatistics");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
    turbStatistics
    (
        const word& name,
        const Time& time,
        const dictionary& dict
    );


    //- Destructor
    virtual ~turbStatistics();

    // Member Functions

        //- Return name of the system call set
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the system calls
        virtual bool read(const dictionary&);

        //- Execute the "executeCalls" at each time-step
        virtual bool execute();

        //- Execute the "endCalls" at the final time-loop
//        virtual void end() {}

        //- Called when time was set at the end of the Time::operator++
//        virtual void timeSet() {}

        //- Write, execute the "writeCalls"
        virtual bool write();

//        //- Update for changes of mesh
//        virtual void updateMesh(const mapPolyMesh&)
//        {}

//        //- Update for changes of mesh
//        virtual void movePoints(const polyMesh&)
//        {}

};

scalar turbStatistics::oldFrac() const
{
    return (1./timeStep_)*(timeStep_-1);
}

scalar turbStatistics::newFrac() const
{
    return 1./timeStep_;
}

//typedef OutputFilterFunctionObject<turbStatistics> turbStatisticsFunctionObject;







} // End namespace functionObjects

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif //TURB_STATISTICS_H
