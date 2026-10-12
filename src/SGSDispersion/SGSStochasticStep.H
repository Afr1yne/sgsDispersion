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

Namespace
    Foam::sgsDispersionDetail

Description
    Local coefficients and exact linear transition for SGS particle dispersion.

SourceFiles
    SGSStochasticStep.C

\*---------------------------------------------------------------------------*/

#ifndef SGSStochasticStep_H
#define SGSStochasticStep_H

#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sgsDispersionDetail
{

//- Exact transition for du = (-a*u + b)*dt + B*dW with frozen coefficients
struct Transition
{
    scalar decay;
    scalar forcing;
    scalar noiseStd;
};


//- Local coefficients of the SGS velocity equation
struct Coefficients
{
    scalar Tp;
    scalar a;
    scalar B2;
    vector forcing;
};


//- Evaluate (exp(x) - 1)/x, including the limit at zero
scalar exprel(const scalar x);

//- Calculate the exact frozen-coefficient transition
Transition transition
(
    const scalar a,
    const scalar b,
    const scalar B2,
    const scalar dt
);

//- Calculate the correlation time along a particle trajectory
scalar particleTime
(
    const scalar variance,
    const scalar epsilon,
    const scalar verticalSpeed,
    const scalar C0 = 4,
    const scalar beta = 2
);

//- Spherical-particle drag correction, Comola et al. (2019), Eq. (11)
scalar dragCorrection(const scalar Re);

//- Assemble Comola et al. (2019), Eqs. (4) and (8)
Coefficients coefficients
(
    const scalar variance,
    const scalar epsilon,
    const scalar fraction,
    const vector& velocity,
    const scalar verticalSpeed,
    const vector& gradient,
    const scalar partialVarianceDt,
    const scalar C0 = 4,
    const scalar beta = 2
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sgsDispersionDetail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
