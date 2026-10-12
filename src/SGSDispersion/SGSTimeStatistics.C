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

\*---------------------------------------------------------------------------*/

#include "SGSTimeStatistics.H"
#include "error.H"
#include <cmath>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::SGSTimeStatistics::removeWeight
(
    const Sample& sample,
    const scalar weight
)
{
    const scalar remaining = duration_ - weight;

    forAll(meanU_, celli)
    {
        const vector delta = sample.U[celli] - meanU_[celli];
        const vector newMean = meanU_[celli] - weight/remaining*delta;

        central2_[celli] = max
        (
            scalar(0),
            central2_[celli]
          - weight*(delta & (sample.U[celli] - newMean))
        );

        meanU_[celli] = newMean;
        meanP_[celli] -= weight/remaining*(sample.P[celli] - meanP_[celli]);
    }

    duration_ = remaining;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SGSTimeStatistics::Sample::Sample
(
    const scalar dt,
    const vectorField& u,
    const scalarField& p
)
:
    weight(dt),
    U(u),
    P(p)
{}


Foam::SGSTimeStatistics::SGSTimeStatistics
(
    const label size,
    const scalar window
)
:
    window_(window),
    duration_(0),
    meanU_(size, Zero),
    central2_(size, 0),
    meanP_(size, 0)
{
    if (!(window_ > 0) || !std::isfinite(window_))
    {
        FatalErrorInFunction
            << "averagingTime must be finite and positive"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SGSTimeStatistics::add
(
    const vectorField& U,
    const scalarField& production,
    const scalar dt
)
{
    if
    (
        U.size() != meanU_.size()
     || production.size() != meanP_.size()
     || !(dt > 0)
     || !std::isfinite(dt)
    )
    {
        FatalErrorInFunction
            << "Invalid statistics sample; fixed cell topology is required"
            << exit(FatalError);
    }

    const scalar total = duration_ + dt;

    forAll(U, celli)
    {
        const vector delta = U[celli] - meanU_[celli];
        const vector newMean = meanU_[celli] + dt/total*delta;

        central2_[celli] += dt*(delta & (U[celli] - newMean));
        meanU_[celli] = newMean;
        meanP_[celli] += dt/total*(production[celli] - meanP_[celli]);
    }

    duration_ = total;
    samples_.push(Sample(dt, U, production));

    // Retain exactly window_ seconds, including a partial oldest interval
    while (duration_ > window_)
    {
        Sample& oldest = samples_.front();
        const scalar removed = min(oldest.weight, duration_ - window_);

        removeWeight(oldest, removed);
        oldest.weight -= removed;

        if (oldest.weight <= SMALL)
        {
            samples_.pop_front();
        }
    }
}


Foam::scalar Foam::SGSTimeStatistics::resolvedK(const label celli) const
{
    return duration_ > VSMALL
      ? max(scalar(0), 0.5*central2_[celli]/duration_)
      : scalar(0);
}

// ************************************************************************* //
