#include "velGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(velGrad, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::velGrad::velGrad
(
    const word& name,
    const objectRegistry&,
    const dictionary& dict,
    const bool
)
:
    name_(name)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velGrad::~velGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velGrad::read(const dictionary& dict)
{
}


void Foam::velGrad::execute()
{

}


void Foam::velGrad::end()
{
   
}


void Foam::velGrad::timeSet()
{
    // Do nothing
}


void Foam::velGrad::write()
{
   
}


// ************************************************************************* //
