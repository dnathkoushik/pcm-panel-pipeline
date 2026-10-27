/*---------------------------------------------------------------------------*\
Class
    Foam::naturalConvectionTemperatureFvPatchScalarField

Description
    Room air on the panel surface as a natural-convection coefficient, for the
    PCM region when the room air is not solved (screening mode).  The heat
    flux into the region is
        q = U (Ta - Tw),   1/U = 1/hc + sum(thicknessLayers/kappaLayers)
        hc = C |Ta - Ts|^n            (air colder than the room air: unstable)
        hc = Cstable |Ta - Ts|^n      (surface warmer: stable, warm face down)
    with Ts the air-side surface temperature (behind the layers).  It is
    applied implicitly as a mixed condition, as externalWallHeatFluxTemperature
    does (valueFraction = U/(U + kappa deltaCoeffs)); hc lags one evaluation.

    Ta is a Function1 of time.  With roomHeatCapacity > 0 the room is closed
    and adiabatic (the paper's test room) and gives up the heat the panel
    takes: Ta(t) = Ta0(t) - int Q dt/roomHeatCapacity, Q = sum(q |Sf|).  The
    running integral is written as `energy`, so restarts continue it.

    kappa is the region's thermo kappa (laminar region).

Usage
    \verbatim
    pcm_to_air
    {
        type            naturalConvectionTemperature;
        Ta              constant 298.15;  // Function1 [K]
        C               4.4;              // [W/m2/K^(1+n)]
        Cstable         2.2;              // optional, default C/2
        n               0.25;             // optional, default 0.25 (laminar)
        thicknessLayers (0.001);          // optional
        kappaLayers     (0.163);
        roomHeatCapacity 30770;           // optional [J/K], default 0
        value           uniform 278.15;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef naturalConvectionTemperatureFvPatchScalarField_H
#define naturalConvectionTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "Function1.H"

namespace Foam
{

class naturalConvectionTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        autoPtr<Function1<scalar>> Ta_;
        scalar C_;
        scalar Cstable_;
        scalar n_;
        scalarList thicknessLayers_;
        scalarList kappaLayers_;
        scalar roomHeatCapacity_;

        //- Heat taken from the room so far [J]
        scalar energy_;

        //- Heat flow into the region at the last evaluation [W]
        scalar Q_;

        label timeIndex_;


    // Private Member Functions

        scalar layerResistance() const;


public:

    TypeName("naturalConvectionTemperature");


    // Constructors

        naturalConvectionTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        naturalConvectionTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        naturalConvectionTemperatureFvPatchScalarField
        (
            const naturalConvectionTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        naturalConvectionTemperatureFvPatchScalarField
        (
            const naturalConvectionTemperatureFvPatchScalarField&
        );

        naturalConvectionTemperatureFvPatchScalarField
        (
            const naturalConvectionTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new naturalConvectionTemperatureFvPatchScalarField(*this)
            );
        }

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new naturalConvectionTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Room air temperature now [K]
        scalar Ta() const;

        virtual void updateCoeffs();

        virtual void write(Ostream&) const;
};


} // End namespace Foam

#endif
