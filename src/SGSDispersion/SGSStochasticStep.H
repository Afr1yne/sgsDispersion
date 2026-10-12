// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef SGS_STOCHASTIC_STEP_H
#define SGS_STOCHASTIC_STEP_H
#include <cmath>
#include <stdexcept>
#include <array>

namespace sgsDispersionDetail
{
// Exact transition for frozen local coefficients: du=(-a*u+b)dt+B*dW.
inline double exprel(const double x)
{
    if (std::abs(x)<1e-5)
        return 1+x/2+x*x/6+x*x*x/24+x*x*x*x/120;
    return std::expm1(x)/x;
}
struct Transition { double decay, forcing, noiseStd; };
inline Transition transition(double a, double b, double B2, double dt)
{
    if (!(dt>=0 && B2>=0) || !std::isfinite(a+b+B2+dt))
        throw std::domain_error("Invalid SGS transition coefficients");
    const double x=-a*dt;
    if (std::abs(x)>300)
        throw std::domain_error("Variance drift unresolved: reduce time step");
    return {std::exp(x), b*dt*exprel(x),
            std::sqrt(B2*dt*exprel(2*x))};
}
inline double particleTime(double variance, double epsilon, double verticalSpeed,
                           double C0=4, double beta=2)
{
    return (2*variance/(C0*epsilon)) /
           std::sqrt(1+std::pow(beta*verticalSpeed,2)/variance);
}

struct Coefficients
{
    double Tp, a, B2;
    std::array<double, 3> forcing;
};

inline double dragCorrection(double Re)
{
    if (!(Re>=0) || !std::isfinite(Re))
        throw std::domain_error("Invalid particle Reynolds number");
    return 1+0.15*std::pow(Re,0.687);
}

// Comola et al. (2019), equations (4) and (8), evaluated on a heavy-particle path.
// The spatial/material derivatives are inputs, independently of their FV discretisation.
inline Coefficients coefficients
(
    double variance, double epsilon, double fraction,
    const std::array<double, 3>& velocity, double verticalSpeed,
    const std::array<double, 3>& gradient, double partialVarianceDt,
    double C0=4, double beta=2
)
{
    if (!(variance>0 && epsilon>0 && fraction>0 && fraction<=1 && C0>0 && beta>=0))
        throw std::domain_error("Invalid SGS model coefficients");
    double material=partialVarianceDt;
    std::array<double,3> forcing;
    for (int i=0; i<3; ++i)
    {
        if (!std::isfinite(velocity[i]) || !std::isfinite(gradient[i]))
            throw std::domain_error("Nonfinite SGS derivative or particle velocity");
        material += velocity[i]*gradient[i];
        forcing[i] = 0.5*gradient[i];
    }
    const double Tp=particleTime(variance,epsilon,verticalSpeed,C0,beta);
    const double a=fraction/Tp-0.5*material/variance;
    const double B2=2*fraction*variance/Tp;
    if (!(Tp>0) || !std::isfinite(Tp+a+B2))
        throw std::domain_error("Nonfinite SGS model coefficients");
    return {Tp,a,B2,forcing};
}
}
#endif
