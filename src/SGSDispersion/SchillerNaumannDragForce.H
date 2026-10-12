/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::SchillerNaumannDragForce

Description
    Spherical-particle drag with f(Re) = 1 + 0.15*Re^0.687.
    Comola et al. (2019), DOI:10.1029/2018JD029614, Eq. (11).

SourceFiles
    SchillerNaumannDragForce.C

\*---------------------------------------------------------------------------*/

#ifndef SchillerNaumannDragForce_H
#define SchillerNaumannDragForce_H

#include "ParticleForce.H"
#include "SGSStochasticStep.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class SchillerNaumannDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SchillerNaumannDragForce
:
    public ParticleForce<CloudType>
{
public:

    //- Runtime type information
    TypeName("schillerNaumannDrag");


    // Constructors

        //- Construct from mesh
        SchillerNaumannDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        SchillerNaumannDragForce
        (
            const SchillerNaumannDragForce<CloudType>& df
        );

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new SchillerNaumannDragForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SchillerNaumannDragForce();


    // Member Functions

        //- Calculate the coupled force
        virtual forceSuSp calcCoupled
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td,
            const scalar dt,
            const scalar mass,
            const scalar Re,
            const scalar muc
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SchillerNaumannDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
