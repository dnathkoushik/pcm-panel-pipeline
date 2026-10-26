/*---------------------------------------------------------------------------*\
    boussinesqBuoyancySource -- see header for the model description.
\*---------------------------------------------------------------------------*/

#include "boussinesqBuoyancySource.H"
#include "basicThermo.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(boussinesqBuoyancySource, 0);
    addToRunTimeSelectionTable(option, boussinesqBuoyancySource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::boussinesqBuoyancySource::boussinesqBuoyancySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    beta_(0),
    TRef_(0),
    g_(Zero),
    UName_("U")
{
    read(dict);

    fieldNames_.resize(1);
    fieldNames_[0] = UName_;
    fv::option::resetApplied();

    Info<< "    boussinesqBuoyancySource: beta " << beta_ << " 1/K, TRef "
        << TRef_ << " K, g " << g_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::boussinesqBuoyancySource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);
    const scalarField& T = thermo.T().primitiveField();
    const scalarField& rhoc = rho.primitiveField();

    volVectorField::Internal Su
    (
        IOobject
        (
            name_ + ":Su",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector(dimForce/dimVolume, Zero)
    );

    for (const label celli : cells_)
    {
        Su[celli] = rhoc[celli]*beta_*(TRef_ - T[celli])*g_;
    }

    eqn += Su;
}


bool Foam::fv::boussinesqBuoyancySource::read(const dictionary& dict)
{
    if (!fv::cellSetOption::read(dict))
    {
        return false;
    }

    coeffs_.readEntry("beta", beta_);
    coeffs_.readEntry("TRef", TRef_);
    coeffs_.readEntry("g", g_);
    UName_ = coeffs_.getOrDefault<word>("U", "U");

    return true;
}


// ************************************************************************* //
