/*---------------------------------------------------------------------------*\
    naturalConvectionTemperature -- see header for the model description.
\*---------------------------------------------------------------------------*/

#include "naturalConvectionTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar
Foam::naturalConvectionTemperatureFvPatchScalarField::layerResistance() const
{
    scalar R = 0;
    forAll(thicknessLayers_, i)
    {
        R += thicknessLayers_[i]/kappaLayers_[i];
    }
    return R;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::naturalConvectionTemperatureFvPatchScalarField::
naturalConvectionTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    Ta_(nullptr),
    C_(0),
    Cstable_(0),
    n_(0.25),
    thicknessLayers_(),
    kappaLayers_(),
    roomHeatCapacity_(0),
    energy_(0),
    Q_(0),
    timeIndex_(-1)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1.0;
}


Foam::naturalConvectionTemperatureFvPatchScalarField::
naturalConvectionTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    Ta_(Function1<scalar>::New("Ta", dict, &db())),
    C_(dict.get<scalar>("C")),
    Cstable_(dict.getOrDefault<scalar>("Cstable", 0.5*C_)),
    n_(dict.getOrDefault<scalar>("n", 0.25)),
    thicknessLayers_(dict.getOrDefault<scalarList>("thicknessLayers", scalarList())),
    kappaLayers_(dict.getOrDefault<scalarList>("kappaLayers", scalarList())),
    roomHeatCapacity_(dict.getOrDefault<scalar>("roomHeatCapacity", 0)),
    energy_(dict.getOrDefault<scalar>("energy", 0)),
    Q_(0),
    timeIndex_(-1)
{
    if (thicknessLayers_.size() != kappaLayers_.size())
    {
        FatalIOErrorInFunction(dict)
            << "thicknessLayers and kappaLayers differ in length"
            << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchScalarField::operator=(patchInternalField());
    }

    if (dict.found("refValue"))
    {
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        refValue() = *this;
        refGrad() = Zero;
        valueFraction() = 1.0;
    }
}


Foam::naturalConvectionTemperatureFvPatchScalarField::
naturalConvectionTemperatureFvPatchScalarField
(
    const naturalConvectionTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    Ta_(ptf.Ta_.clone()),
    C_(ptf.C_),
    Cstable_(ptf.Cstable_),
    n_(ptf.n_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    roomHeatCapacity_(ptf.roomHeatCapacity_),
    energy_(ptf.energy_),
    Q_(ptf.Q_),
    timeIndex_(ptf.timeIndex_)
{}


Foam::naturalConvectionTemperatureFvPatchScalarField::
naturalConvectionTemperatureFvPatchScalarField
(
    const naturalConvectionTemperatureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    Ta_(ptf.Ta_.clone()),
    C_(ptf.C_),
    Cstable_(ptf.Cstable_),
    n_(ptf.n_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    roomHeatCapacity_(ptf.roomHeatCapacity_),
    energy_(ptf.energy_),
    Q_(ptf.Q_),
    timeIndex_(ptf.timeIndex_)
{}


Foam::naturalConvectionTemperatureFvPatchScalarField::
naturalConvectionTemperatureFvPatchScalarField
(
    const naturalConvectionTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    Ta_(ptf.Ta_.clone()),
    C_(ptf.C_),
    Cstable_(ptf.Cstable_),
    n_(ptf.n_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    roomHeatCapacity_(ptf.roomHeatCapacity_),
    energy_(ptf.energy_),
    Q_(ptf.Q_),
    timeIndex_(ptf.timeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::naturalConvectionTemperatureFvPatchScalarField::Ta() const
{
    const scalar Ta0 = Ta_->value(db().time().timeOutputValue());
    return roomHeatCapacity_ > 0 ? Ta0 - energy_/roomHeatCapacity_ : Ta0;
}


void Foam::naturalConvectionTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // First evaluation of a new time step: book the last step's heat
    const label ti = db().time().timeIndex();
    if (ti != timeIndex_)
    {
        if (timeIndex_ >= 0)
        {
            energy_ += Q_*db().time().deltaT0Value();
        }
        timeIndex_ = ti;
    }

    const basicThermo& thermo =
        db().lookupObject<basicThermo>(basicThermo::dictName);
    const scalarField kappaw(thermo.kappa(patch().index()));

    const scalar Ta = this->Ta();
    const scalar R = layerResistance();
    const scalarField& Tw = *this;

    scalarField U(size());
    forAll(U, facei)
    {
        const scalar dTw = Ta - Tw[facei];
        const scalar Cdir = dTw >= 0 ? C_ : Cstable_;

        // Air-side difference: the layers take part of Ta - Tw
        scalar hc = Cdir*pow(mag(dTw), n_);
        for (label iter = 0; iter < 3; ++iter)
        {
            hc = Cdir*pow(mag(dTw)/(1 + hc*R), n_);
        }
        U[facei] = 1/(1/max(hc, SMALL) + R);
    }

    refValue() = Ta;
    refGrad() = Zero;
    valueFraction() = U/(U + kappaw*patch().deltaCoeffs());

    Q_ = gSum(U*(Ta - Tw)*patch().magSf());

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::naturalConvectionTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    Ta_->writeData(os);
    os.writeEntry("C", C_);
    os.writeEntry("Cstable", Cstable_);
    os.writeEntry("n", n_);
    if (thicknessLayers_.size())
    {
        // List::writeEntry tags binary lists (List<scalar>) so that a dictionary
        // can read them back; os.writeEntry does not (decomposePar writes binary)
        thicknessLayers_.writeEntry("thicknessLayers", os);
        kappaLayers_.writeEntry("kappaLayers", os);
    }
    os.writeEntry("roomHeatCapacity", roomHeatCapacity_);
    os.writeEntry("energy", energy_);
}


// * * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        naturalConvectionTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
