/*---------------------------------------------------------------------------*\
    ParticleForce interface follows OpenFOAM's SphereDragForce:
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
    SPDX-License-Identifier: GPL-3.0-or-later

    Optional companion force: Comola et al. (2019), equation (11), spherical case.
\*---------------------------------------------------------------------------*/
#ifndef SchillerNaumannDragForce_H
#define SchillerNaumannDragForce_H

#include "ParticleForce.H"
#include "SGSStochasticStep.H"

namespace Foam
{
template<class CloudType>
class SchillerNaumannDragForce : public ParticleForce<CloudType>
{
public:
    TypeName("schillerNaumannDrag");

    SchillerNaumannDragForce
    (
        CloudType& owner, const fvMesh& mesh, const dictionary& dict
    )
    : ParticleForce<CloudType>(owner, mesh, dict, typeName, false)
    {}

    SchillerNaumannDragForce(const SchillerNaumannDragForce& rhs)
    : ParticleForce<CloudType>(rhs)
    {}

    virtual autoPtr<ParticleForce<CloudType>> clone() const
    {
        return autoPtr<ParticleForce<CloudType>>
        (
            new SchillerNaumannDragForce(*this)
        );
    }

    virtual forceSuSp calcCoupled
    (
        const typename CloudType::parcelType& p,
        const typename CloudType::parcelType::trackingData&,
        const scalar, const scalar mass, const scalar Re, const scalar muc
    ) const
    {
        // m/tp, tp = rho_p*d_p^2/(18*mu*f(Re)); Eq. (11).
        // Keep the paper expression; do not silently substitute high-Re Cd=0.44.
        return forceSuSp
        (
            Zero,
            mass*18*muc*sgsDispersionDetail::dragCorrection(Re)/(p.rho()*sqr(p.d()))
        );
    }
};
}
#endif
