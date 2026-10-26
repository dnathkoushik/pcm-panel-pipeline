/*---------------------------------------------------------------------------*\
    pcmEnthalpyPorositySource -- see header for the model description.
\*---------------------------------------------------------------------------*/

#include "pcmEnthalpyPorositySource.H"
#include "basicThermo.H"
#include "fvmSup.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(pcmEnthalpyPorositySource, 0);
    addToRunTimeSelectionTable(option, pcmEnthalpyPorositySource, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

inline Foam::scalar Foam::fv::pcmEnthalpyPorositySource::gamma
(
    const scalar T
) const
{
    return min(max((T - Tsol_)/(Tliq_ - Tsol_), scalar(0)), scalar(1));
}


inline Foam::scalar Foam::fv::pcmEnthalpyPorositySource::cpTrue
(
    const scalar T
) const
{
    return CpS_ + (CpL_ - CpS_)*gamma(T);
}


inline Foam::scalar Foam::fv::pcmEnthalpyPorositySource::kappaTrue
(
    const scalar T
) const
{
    return kappaS_ + (kappaL_ - kappaS_)*gamma(T);
}


inline Foam::scalar Foam::fv::pcmEnthalpyPorositySource::Ecell
(
    const scalar T,
    const scalar CpRef
) const
{
    // E(T) = [h_sen(T) - h_sen(Tliq)] - CpRef (T - Tliq) + gamma L
    const scalar dT = Tliq_ - Tsol_;

    if (T >= Tliq_)
    {
        return (CpL_ - CpRef)*(T - Tliq_) + L_;
    }
    else if (T > Tsol_)
    {
        const scalar g = (T - Tsol_)/dT;
        return
            (CpRef - CpS_)*(Tliq_ - T)
          - 0.5*(CpL_ - CpS_)*dT*(1 - sqr(g))
          + g*L_;
    }

    return
        CpRef*(Tliq_ - T)
      - CpS_*(Tsol_ - T)
      - 0.5*(CpS_ + CpL_)*dT;
}


void Foam::fv::pcmEnthalpyPorositySource::update()
{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    const scalarField& T = thermo.T().primitiveField();
    const tmp<volScalarField> tCp(thermo.Cp());
    const tmp<volScalarField> tkappa(thermo.kappa());
    const scalarField& Cp = tCp().primitiveField();
    const scalarField& kappa = tkappa().primitiveField();

    scalarField& g = liquidFraction_.primitiveFieldRef();
    scalarField& E = E_.primitiveFieldRef();
    scalarField& dEdT = dEdT_.primitiveFieldRef();
    scalarField& kEx = kappaExcess_.primitiveFieldRef();

    const scalar latentSlope = L_/(Tliq_ - Tsol_);

    for (const label celli : cells_)
    {
        const scalar Tc = T[celli];
        g[celli] = gamma(Tc);
        E[celli] = Ecell(Tc, Cp[celli]);
        dEdT[celli] =
            cpTrue(Tc) - Cp[celli]
          + ((Tc > Tsol_ && Tc < Tliq_) ? latentSlope : scalar(0));
        kEx[celli] = kappaTrue(Tc) - kappa[celli];
    }

    liquidFraction_.correctBoundaryConditions();
    E_.correctBoundaryConditions();
    // kappaExcess_ keeps zero boundary values: no extra interface flux

    // First evaluation in a new time step: T is still the converged
    // solution of the previous step, so this is E at the old time level
    const label timeIndex = mesh_.time().timeIndex();
    if (E0TimeIndex_ != timeIndex)
    {
        E0_ = E;
        E0TimeIndex_ = timeIndex;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::pcmEnthalpyPorositySource::pcmEnthalpyPorositySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    Tsol_(0),
    Tliq_(0),
    L_(0),
    Amush_(0),
    eps_(0),
    CpS_(0),
    CpL_(0),
    kappaS_(0),
    kappaL_(0),
    latentConvection_(true),
    energyName_("h"),
    UName_("U"),
    phiName_("phi"),
    liquidFraction_
    (
        IOobject
        (
            "liquidFraction",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 1.0),
        "zeroGradient"
    ),
    E_
    (
        IOobject
        (
            "pcmExtraEnthalpy",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass, Zero),
        "zeroGradient"
    ),
    dEdT_
    (
        IOobject
        (
            "pcmExtraEnthalpyDT",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero),
        "zeroGradient"
    ),
    kappaExcess_
    (
        IOobject
        (
            "pcmKappaExcess",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPower/dimLength/dimTemperature, Zero),
        "calculated"
    ),
    E0_(mesh.nCells(), Zero),
    E0TimeIndex_(-1),
    kinkLimiter_(true),
    TLin_(mesh.nCells(), Zero),
    heLin_(mesh.nCells(), Zero)
{
    read(dict);

    fieldNames_.resize(2);
    fieldNames_[0] = energyName_;
    fieldNames_[1] = UName_;
    fv::option::resetApplied();

    update();

    Info<< "    pcmEnthalpyPorositySource: Tsolidus " << Tsol_
        << " K, Tliquidus " << Tliq_ << " K, L " << L_
        << " J/kg, Amush " << Amush_ << ", epsilon " << eps_ << nl
        << "    initial liquid fraction "
        << gSum(liquidFraction_.primitiveField()*mesh.V().field())
          /gSum(mesh.V().field())
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::pcmEnthalpyPorositySource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    update();

    const volScalarField& he = eqn.psi();
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);
    const tmp<volScalarField> tCp(thermo.Cp());
    const scalarField& Cp = tCp().primitiveField();

    const scalar rDeltaT = 1.0/mesh_.time().deltaTValue();
    const scalarField& rhoNew = rho.primitiveField();
    const scalarField& rhoOld = rho.oldTime().primitiveField();

    volScalarField::Internal Su
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
        dimensionedScalar(dimPower/dimVolume, Zero)
    );
    volScalarField::Internal Sp
    (
        IOobject
        (
            name_ + ":Sp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    );

    const scalarField& T = thermo.T().primitiveField();

    // -d(rho E)/dt with E^{k+1} = E^k + dE/dT (h - h^k)/CpRef  (Euler)
    for (const label celli : cells_)
    {
        const scalar spc = rhoNew[celli]*dEdT_[celli]*rDeltaT/Cp[celli];
        Su[celli] =
          - (rhoNew[celli]*E_[celli] - rhoOld[celli]*E0_[celli])*rDeltaT
          + spc*he[celli];
        Sp[celli] = spc;

        TLin_[celli] = T[celli];
        heLin_[celli] = he[celli];
    }

    eqn += Su;
    eqn -= fvm::Sp(Sp, he);

    if (latentConvection_)
    {
        const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>(phiName_);
        eqn -= fvc::div(phi, E_);
    }

    eqn += fvc::laplacian(kappaExcess_, thermo.T());
}


void Foam::fv::pcmEnthalpyPorositySource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    update();

    volScalarField::Internal D
    (
        IOobject
        (
            name_ + ":D",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    );

    for (const label celli : cells_)
    {
        const scalar g = liquidFraction_[celli];
        D[celli] = Amush_*sqr(1 - g)/(pow3(g) + eps_);
    }

    // Eq. (21) as a sink: S = -D u
    eqn -= fvm::Sp(D, eqn.psi());
}


void Foam::fv::pcmEnthalpyPorositySource::correct(volScalarField& he)
{
    if (!kinkLimiter_ || he.name() != energyName_)
    {
        return;
    }

    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);
    const tmp<volScalarField> tCp(thermo.Cp());
    const scalarField& Cp = tCp().primitiveField();

    scalarField& hef = he.primitiveFieldRef();

    // Within tol of a kink a cell moves freely; a limited cell is placed
    // tol past the kink so the next linearisation sees the new phase
    const scalar tol = 1e-3*(Tliq_ - Tsol_);
    const scalar kinks[2] = {Tsol_, Tliq_};

    label nLimited = 0;

    for (const label celli : cells_)
    {
        const scalar T0 = TLin_[celli];
        const scalar Tn = T0 + (hef[celli] - heLin_[celli])/Cp[celli];
        scalar Tlim = Tn;

        for (const scalar kink : kinks)
        {
            if (T0 < kink - tol && Tlim > kink)
            {
                Tlim = kink + tol;
            }
            else if (T0 > kink + tol && Tlim < kink)
            {
                Tlim = kink - tol;
            }
        }

        if (Tlim != Tn)
        {
            hef[celli] = heLin_[celli] + Cp[celli]*(Tlim - T0);
            ++nLimited;
        }
    }

    if (returnReduce(nLimited, sumOp<label>()))
    {
        he.correctBoundaryConditions();

        if (debug)
        {
            Info<< "    pcmEnthalpyPorositySource: limited "
                << returnReduce(nLimited, sumOp<label>())
                << " cells at the solidus/liquidus" << endl;
        }
    }
}


bool Foam::fv::pcmEnthalpyPorositySource::read(const dictionary& dict)
{
    if (!fv::cellSetOption::read(dict))
    {
        return false;
    }

    coeffs_.readEntry("Tsolidus", Tsol_);
    coeffs_.readEntry("Tliquidus", Tliq_);
    coeffs_.readEntry("L", L_);
    coeffs_.readEntry("Amush", Amush_);
    coeffs_.readEntry("epsilon", eps_);
    coeffs_.readEntry("CpSolid", CpS_);
    coeffs_.readEntry("CpLiquid", CpL_);
    coeffs_.readEntry("kappaSolid", kappaS_);
    coeffs_.readEntry("kappaLiquid", kappaL_);
    latentConvection_ = coeffs_.getOrDefault<bool>("latentConvection", true);
    kinkLimiter_ = coeffs_.getOrDefault<bool>("kinkLimiter", true);
    energyName_ = coeffs_.getOrDefault<word>("energyField", "h");
    UName_ = coeffs_.getOrDefault<word>("U", "U");
    phiName_ = coeffs_.getOrDefault<word>("phi", "phi");

    if (Tliq_ <= Tsol_)
    {
        FatalIOErrorInFunction(coeffs_)
            << "Tliquidus (" << Tliq_ << ") must exceed Tsolidus ("
            << Tsol_ << ")" << exit(FatalIOError);
    }

    return true;
}


// ************************************************************************* //
