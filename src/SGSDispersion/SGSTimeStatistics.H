/*---------------------------------------------------------------------------*\
    Internal time-window statistics for the SGS dispersion model.
    SPDX-License-Identifier: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#ifndef SGSTimeStatistics_H
#define SGSTimeStatistics_H

#include "vectorField.H"
#include "scalarField.H"
#include "error.H"
#include <deque>
#include <cmath>

namespace Foam
{

class SGSTimeStatistics
{
    struct Sample
    {
        scalar weight;
        vectorField U;
        scalarField P;

        Sample(const scalar dt, const vectorField& u, const scalarField& p)
        : weight(dt), U(u), P(p)
        {}
    };

    scalar window_;
    scalar duration_;
    vectorField meanU_;
    scalarField central2_;
    scalarField meanP_;
    std::deque<Sample> samples_;

    void removeWeight(const Sample& sample, const scalar weight)
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
            meanP_[celli] -=
                weight/remaining*(sample.P[celli] - meanP_[celli]);
        }

        duration_ = remaining;
    }

public:

    SGSTimeStatistics(const label size, const scalar window)
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

    void add
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
        samples_.emplace_back(dt, U, production);

        // Retain exactly window_ seconds, including a partial oldest interval.
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

    scalar duration() const { return duration_; }
    label sampleCount() const { return label(samples_.size()); }
    const vectorField& meanU() const { return meanU_; }
    const scalarField& meanProduction() const { return meanP_; }

    scalar resolvedK(const label celli) const
    {
        return duration_ > VSMALL
            ? max(scalar(0), 0.5*central2_[celli]/duration_)
            : scalar(0);
    }
};

} // End namespace Foam

#endif
