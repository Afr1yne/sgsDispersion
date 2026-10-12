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
    Foam::SGSTimeStatistics

Description
    Time-weighted velocity and production statistics over a sliding window.

SourceFiles
    SGSTimeStatistics.C

\*---------------------------------------------------------------------------*/

#ifndef SGSTimeStatistics_H
#define SGSTimeStatistics_H

#include "vectorField.H"
#include "scalarField.H"
#include "FIFOStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class SGSTimeStatistics Declaration
\*---------------------------------------------------------------------------*/

class SGSTimeStatistics
{
    // Private Data Types

        //- A velocity and production sample with its time weight
        struct Sample
        {
            scalar weight;
            vectorField U;
            scalarField P;

            Sample
            (
                const scalar dt,
                const vectorField& u,
                const scalarField& p
            );
        };


    // Private Data

        //- Requested averaging duration
        scalar window_;

        //- Duration currently represented by the samples
        scalar duration_;

        //- Time-weighted mean velocity
        vectorField meanU_;

        //- Time-weighted central second moment
        scalarField central2_;

        //- Time-weighted mean production
        scalarField meanP_;

        //- Samples ordered from oldest to newest
        FIFOStack<Sample> samples_;


    // Private Member Functions

        //- Remove a complete or partial sample weight
        void removeWeight(const Sample& sample, const scalar weight);


public:

    // Constructors

        //- Construct for a fixed number of cells and averaging duration
        SGSTimeStatistics(const label size, const scalar window);


    // Member Functions

        //- Add a sample representing dt seconds
        void add
        (
            const vectorField& U,
            const scalarField& production,
            const scalar dt
        );

        //- Duration represented by the current samples
        scalar duration() const
        {
            return duration_;
        }

        //- Number of stored samples
        label sampleCount() const
        {
            return samples_.size();
        }

        //- Mean velocity
        const vectorField& meanU() const
        {
            return meanU_;
        }

        //- Mean production
        const scalarField& meanProduction() const
        {
            return meanP_;
        }

        //- Resolved turbulent kinetic energy in a cell
        scalar resolvedK(const label celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
