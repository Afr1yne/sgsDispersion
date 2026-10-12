/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Afr1yne
-------------------------------------------------------------------------------
License
    This file is part of sgsDispersion, an OpenFOAM user library.

    sgsDispersion is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    sgsDispersion is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with sgsDispersion.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::SGSDispersion

Description
    SGS velocity fluctuations for Lagrangian particles in incompressible LES.
    Comola et al. (2019), DOI:10.1029/2018JD029614, Eqs. (4)-(8).

SourceFiles
    SGSDispersion.C

\*---------------------------------------------------------------------------*/

#ifndef SGSDispersion_H
#define SGSDispersion_H

#include "DispersionModel.H"
#include "volFields.H"
#include "SGSTimeStatistics.H"
#include "SGSStochasticStep.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SGSDispersion Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SGSDispersion
:
    public DispersionModel<CloudType>
{
    // Private Data

        //- Lagrangian timescale coefficient
        scalar C0_;

        //- Particle crossing-trajectory coefficient
        scalar beta_;

        //- SGS dissipation coefficient
        scalar cEpsilon_;

        //- Unit vector in the vertical direction
        vector verticalDirection_;

        //- Initialise SGS velocity from the stationary distribution
        bool stationaryInitial_;

        //- Internal time-weighted statistics
        SGSTimeStatistics statistics_;

        //- Time of the last completed flow sample
        scalar lastStatisticsTime_;

        //- Last time index processed by cacheFields
        label lastTimeIndex_;

        //- A previous SGS variance is available
        bool havePreviousVariance_;

        //- The full averaging duration has been reported
        bool fullWindowReported_;

        //- SGS variance from the previous sample
        scalarField previousVariance_;

        //- SGS velocity variance
        autoPtr<volScalarField> sigma2_;

        //- SGS fraction of the turbulent kinetic energy
        autoPtr<volScalarField> alpha_;

        //- Mean energy dissipation
        autoPtr<volScalarField> epsilon_;

        //- Resolved turbulent kinetic energy
        autoPtr<volScalarField> resolvedK_;

        //- Averaging duration represented by the current statistics
        autoPtr<volScalarField> windowDuration_;

        //- Partial time derivative of the SGS variance
        autoPtr<volScalarField> dSigma2_;

        //- Spatial gradient of the SGS variance
        autoPtr<volVectorField> gradSigma2_;


    // Private Member Functions

        //- Construct a diagnostic scalar field
        autoPtr<volScalarField> makeScalarField
        (
            const word& name,
            const dimensionSet& dims,
            const IOobject::writeOption write
        ) const;


public:

    //- Runtime type information
    TypeName("sgsDispersion");


    // Constructors

        //- Construct from components
        SGSDispersion(const dictionary& dict, CloudType& owner);

        //- Construct copy
        SGSDispersion(const SGSDispersion<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DispersionModel<CloudType>> clone() const
        {
            return autoPtr<DispersionModel<CloudType>>
            (
                new SGSDispersion<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SGSDispersion();


    // Member Functions

        //- Cache flow statistics for particle tracking
        virtual void cacheFields(const bool store);

        //- Update the SGS velocity and return the particle-seen flow velocity
        virtual vector update
        (
            const scalar dt,
            const label celli,
            const vector& U,
            const vector& Uc,
            vector& UTurb,
            scalar& tTurb
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SGSDispersion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
