/*---------------------------------------------------------------------------*\
Class
    Foam::fv::pcmEnthalpyPorositySource

Description
    Enthalpy-porosity phase change for a Boussinesq fluid region of
    chtMultiRegionFoam, following Yang et al., Energy & Buildings 285 (2023)
    112898, Sec. 3.2-3.3.  Written for OpenFOAM (openfoam.com) v2306+.

    Liquid fraction, Eq. (1):
        gamma = clamp((T - Tsolidus)/(Tliquidus - Tsolidus), 0, 1)

    Energy, Eqs. (15), (18)-(20):
        d(rho H)/dt + div(rho u H) = div(k grad T),
        H = h_sen(T) + gamma L,  h_sen = int cp(T) dT
    cp and k are constant in each phase and vary linearly in the mushy zone
    (assumption 1 of Sec. 3.2).  The thermo package solves for
    h = CpRef (T - Tstd) with constant CpRef and kappaRef, so this source adds
        S_h = - d(rho E)/dt - div(phi E) + div((k(T) - kappaRef) grad T)
    with E(T) = h_sen(T) - CpRef T + gamma L (+ const). d(rho E)/dt is
    Newton-linearised in T (Voller & Swaminathan) so the PIMPLE outer
    correctors converge the latent heat; time discretisation is Euler.
    The excess conductivity is zero on boundary faces so the coupled
    interface flux uses kappaRef (set it to the liquid value).

    Momentum sink, Eq. (21) (applied as a sink, i.e. with a minus sign):
        S_U = - Amush (1 - gamma)^2/(gamma^3 + eps) U

Usage
    \verbatim
    pcmPhaseChange
    {
        type            pcmEnthalpyPorositySource;
        selectionMode   all;
        energyField     h;
        Tsolidus        286;
        Tliquidus       293;
        L               115000;
        Amush           1e6;
        epsilon         0.001;
        CpSolid         1830;
        CpLiquid        2210;
        kappaSolid      0.183;
        kappaLiquid     0.169;
        latentConvection true;
    }
    \endverbatim
    Written fields: liquidFraction.

\*---------------------------------------------------------------------------*/

#ifndef pcmEnthalpyPorositySource_H
#define pcmEnthalpyPorositySource_H

#include "cellSetOption.H"
#include "volFields.H"

namespace Foam
{
namespace fv
{

class pcmEnthalpyPorositySource
:
    public fv::cellSetOption
{
    // Private Data

        scalar Tsol_;
        scalar Tliq_;
        scalar L_;
        scalar Amush_;
        scalar eps_;
        scalar CpS_;
        scalar CpL_;
        scalar kappaS_;
        scalar kappaL_;
        bool latentConvection_;
        word energyName_;
        word UName_;
        word phiName_;

        //- Liquid fraction gamma, Eq. (1)
        volScalarField liquidFraction_;

        //- Extra enthalpy E(T) [J/kg]
        volScalarField E_;

        //- dE/dT [J/kg/K]
        volScalarField dEdT_;

        //- k(T) - kappaRef [W/m/K], zero on boundaries
        volScalarField kappaExcess_;

        //- E at the end of the previous time step
        scalarField E0_;

        //- Time index of E0_
        label E0TimeIndex_;


    // Private Member Functions

        //- Evaluate gamma, E, dE/dT and the excess conductivity from T
        void update();

        //- Cell-wise model functions
        inline scalar gamma(const scalar T) const;
        inline scalar cpTrue(const scalar T) const;
        inline scalar kappaTrue(const scalar T) const;
        inline scalar Ecell(const scalar T, const scalar CpRef) const;


public:

    TypeName("pcmEnthalpyPorositySource");

    pcmEnthalpyPorositySource
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~pcmEnthalpyPorositySource() = default;


    // Member Functions

        //- Energy source (compressible form)
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Momentum sink (compressible form)
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        virtual bool read(const dictionary& dict);
};


} // End namespace fv
} // End namespace Foam

#endif
