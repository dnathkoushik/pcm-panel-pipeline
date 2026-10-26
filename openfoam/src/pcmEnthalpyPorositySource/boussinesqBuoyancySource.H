/*---------------------------------------------------------------------------*\
Class
    Foam::fv::boussinesqBuoyancySource

Description
    Boussinesq buoyancy as a momentum source, for a fluid region of
    chtMultiRegionFoam whose thermo uses a constant density (rhoConst):
        S_U = rho0 beta (TRef - T) g
    Density is then constant in every equation except this force, which is
    the Boussinesq approximation of Yang et al. (2023), Sec. 3 (as in
    Fluent).  The Boussinesq equation of state instead puts rho(T) into the
    continuity and energy equations too; in a closed region (sealed PCM
    capsule, closed room) that leaves the total mass unconserved, and the
    resulting continuity errors appear as spurious heat sources (CI showed
    one PCM cell cooling ~40 K/s with nothing colder around it).

Usage
    \verbatim
    boussinesq
    {
        type            boussinesqBuoyancySource;
        selectionMode   all;
        beta            8.5e-4;         // [1/K]
        TRef            289.5;          // [K]
        g               (0 0 -9.81);    // [m/s2]
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef boussinesqBuoyancySource_H
#define boussinesqBuoyancySource_H

#include "cellSetOption.H"

namespace Foam
{
namespace fv
{

class boussinesqBuoyancySource
:
    public fv::cellSetOption
{
    // Private Data

        scalar beta_;
        scalar TRef_;
        vector g_;
        word UName_;


public:

    TypeName("boussinesqBuoyancySource");

    boussinesqBuoyancySource
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~boussinesqBuoyancySource() = default;


    // Member Functions

        //- Momentum source (compressible form)
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
