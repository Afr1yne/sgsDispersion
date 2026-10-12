/*---------------------------------------------------------------------------*\
    SGS dispersion model for Lagrangian particles in incompressible LES.
    Model equations: Comola et al. (2019), doi:10.1029/2018JD029614, (4)-(8).
    SPDX-License-Identifier: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#ifndef SGSDispersion_H
#define SGSDispersion_H

#include "DispersionModel.H"
#include "volFields.H"
#include "SGSTimeStatistics.H"
#include "SGSStochasticStep.H"

namespace Foam
{

template<class CloudType>
class SGSDispersion : public DispersionModel<CloudType>
{
    scalar C0_;
    scalar beta_;
    scalar cEpsilon_;
    vector verticalDirection_;
    bool stationaryInitial_;

    SGSTimeStatistics statistics_;
    scalar lastStatisticsTime_;
    label lastTimeIndex_;
    bool havePreviousVariance_;
    bool fullWindowReported_;
    scalarField previousVariance_;

    autoPtr<volScalarField> sigma2_;
    autoPtr<volScalarField> alpha_;
    autoPtr<volScalarField> epsilon_;
    autoPtr<volScalarField> resolvedK_;
    autoPtr<volScalarField> windowDuration_;
    autoPtr<volScalarField> dSigma2_;
    autoPtr<volVectorField> gradSigma2_;

public:

    TypeName("sgsDispersion");

    SGSDispersion(const dictionary&, CloudType&);
    SGSDispersion(const SGSDispersion&);

    virtual autoPtr<DispersionModel<CloudType>> clone() const
    {
        return autoPtr<DispersionModel<CloudType>>
        (
            new SGSDispersion(*this)
        );
    }

    virtual void cacheFields(const bool store);

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

} // End namespace Foam

#ifdef NoRepository
    #include "SGSDispersion.C"
#endif

#endif
