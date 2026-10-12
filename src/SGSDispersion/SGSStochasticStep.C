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

#include "SGSStochasticStep.H"
#include "error.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::sgsDispersionDetail::exprel(const scalar x)
{
    if (mag(x) < 1e-5)
    {
        return 1 + x/2 + x*x/6 + x*x*x/24 + x*x*x*x/120;
    }

    return std::expm1(x)/x;
}


Foam::sgsDispersionDetail::Transition Foam::sgsDispersionDetail::transition
(
    const scalar a,
    const scalar b,
    const scalar B2,
    const scalar dt
)
{
    if (!(dt >= 0 && B2 >= 0) || !std::isfinite(a + b + B2 + dt))
    {
        FatalErrorInFunction
            << "Invalid SGS transition coefficients" << exit(FatalError);
    }

    const scalar x = -a*dt;

    if (mag(x) > 300)
    {
        FatalErrorInFunction
            << "Variance drift unresolved: reduce time step"
            << exit(FatalError);
    }

    return
    {
        std::exp(x),
        b*dt*exprel(x),
        std::sqrt(B2*dt*exprel(2*x))
    };
}


Foam::scalar Foam::sgsDispersionDetail::particleTime
(
    const scalar variance,
    const scalar epsilon,
    const scalar verticalSpeed,
    const scalar C0,
    const scalar beta
)
{
    return
        (2*variance/(C0*epsilon))
       /std::sqrt(1 + std::pow(beta*verticalSpeed, 2)/variance);
}


Foam::scalar Foam::sgsDispersionDetail::dragCorrection(const scalar Re)
{
    if (!(Re >= 0) || !std::isfinite(Re))
    {
        FatalErrorInFunction
            << "Invalid particle Reynolds number" << exit(FatalError);
    }

    return 1 + 0.15*std::pow(Re, 0.687);
}


Foam::sgsDispersionDetail::Coefficients
Foam::sgsDispersionDetail::coefficients
(
    const scalar variance,
    const scalar epsilon,
    const scalar fraction,
    const vector& velocity,
    const scalar verticalSpeed,
    const vector& gradient,
    const scalar partialVarianceDt,
    const scalar C0,
    const scalar beta
)
{
    if
    (
        !(variance > 0 && epsilon > 0 && fraction > 0 && fraction <= 1)
     || !(C0 > 0 && beta >= 0)
    )
    {
        FatalErrorInFunction
            << "Invalid SGS model coefficients" << exit(FatalError);
    }

    scalar material = partialVarianceDt;
    vector forcing(Zero);

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        if
        (
            !std::isfinite(velocity[cmpt])
         || !std::isfinite(gradient[cmpt])
        )
        {
            FatalErrorInFunction
                << "Nonfinite SGS derivative or particle velocity"
                << exit(FatalError);
        }

        material += velocity[cmpt]*gradient[cmpt];
        forcing[cmpt] = 0.5*gradient[cmpt];
    }

    const scalar Tp = particleTime(variance, epsilon, verticalSpeed, C0, beta);
    const scalar a = fraction/Tp - 0.5*material/variance;
    const scalar B2 = 2*fraction*variance/Tp;

    if (!(Tp > 0) || !std::isfinite(Tp + a + B2))
    {
        FatalErrorInFunction
            << "Nonfinite SGS model coefficients" << exit(FatalError);
    }

    return {Tp, a, B2, forcing};
}

// ************************************************************************* //
