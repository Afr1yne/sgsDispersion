/*---------------------------------------------------------------------------*\
    Runtime selection registration for the classic basicKinematicCloud.
    SPDX-License-Identifier: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "basicKinematicCloud.H"
#include "SGSDispersion.H"
#include "SchillerNaumannDragForce.H"

namespace Foam
{
    makeDispersionModelType(SGSDispersion, basicKinematicCloud);
    makeParticleForceModelType(SchillerNaumannDragForce, basicKinematicCloud);
}
