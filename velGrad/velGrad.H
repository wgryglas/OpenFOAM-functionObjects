#ifndef velGrad_H
#define velGrad_H

#include "fvCFD.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes


//class objectRegistry;
//class dictionary;
//class word;
//class mapPolyMesh;
//class polyMesh;
//class volTensorField;

/*---------------------------------------------------------------------------*\
                         Class velGrad Declaration
\*---------------------------------------------------------------------------*/



class velGrad
{
protected:

    // Private data

        //- Name of this calculator
        word name_;

        tmp<volTensorField> gradU;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        velGrad(const velGrad&);

        //- Disallow default bitwise assignment
        void operator=(const velGrad&);


public:

    //- Runtime type information
    TypeName("velGrad");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        velGrad
        (
            const word& name,
            const objectRegistry& unused,
            const dictionary&,
            const bool loadFromFilesUnused = false
        );


    //- Destructor
    virtual ~velGrad();

    // Member Functions

        //- Return name of the system call set
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the system calls
        virtual void read(const dictionary&);

        //- Execute the "executeCalls" at each time-step
        virtual void execute();

        //- Execute the "endCalls" at the final time-loop
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write, execute the "writeCalls"
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
