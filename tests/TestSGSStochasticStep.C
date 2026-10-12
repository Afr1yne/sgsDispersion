// SPDX-License-Identifier: GPL-3.0-or-later
#include "../src/SGSDispersion/SGSStochasticStep.H"
#include <random>
#include <iostream>
#include <cstdlib>
void require(bool condition,const char* label)
{ if(!condition) { std::cerr<<"FAIL "<<label<<'\n'; std::exit(1); } }
int main()
{
    const double var=0.4,eps=0.08,alpha=0.2;
    const double Tp=sgsDispersionDetail::particleTime(var,eps,0);
    require(std::abs(Tp-2*var/(4*eps))<1e-14,"fluid timescale");
    require(sgsDispersionDetail::particleTime(var,eps,0.5)<Tp,"settling reduces timescale");
    const double a=alpha/Tp,B2=2*alpha*var/Tp;
    for(double dt:{1e-7,0.01,1.,100.})
    {
        auto q=sgsDispersionDetail::transition(a,0,B2,dt);
        require(std::abs(q.decay*q.decay*var+q.noiseStd*q.noiseStd-var)<1e-12,
                "stationary variance independent of timestep");
    }
    auto brownian=sgsDispersionDetail::transition(0,0.3,0.8,0.25);
    require(std::abs(brownian.forcing-0.075)<1e-14,"zero-drift limit");
    require(std::abs(brownian.noiseStd*brownian.noiseStd-0.2)<1e-14,"Brownian variance");
    auto growing=sgsDispersionDetail::transition(-0.1,0.2,0.8,1);
    require(growing.decay>1 && std::isfinite(growing.noiseStd),"variance growth drift");
    std::mt19937_64 rng(7419); std::normal_distribution<double> normal;
    const int N=200000,steps=50; const double dt=0.2;
    auto q=sgsDispersionDetail::transition(a,0,B2,dt);
    double sum=0,sq=0,cross=0,startSq=0;
    for(int i=0;i<N;++i)
    {
        double u=std::sqrt(var)*normal(rng),u0=u;
        for(int k=0;k<steps;++k)u=q.decay*u+q.noiseStd*normal(rng);
        sum+=u; sq+=u*u; cross+=u*u0; startSq+=u0*u0;
    }
    const double mean=sum/N,variance=sq/N-mean*mean,corr=cross/startSq;
    const double expected=std::exp(-a*steps*dt);
    require(std::abs(mean)<0.006,"ensemble mean");
    require(std::abs(variance/var-1)<0.015,"ensemble variance");
    require(std::abs(corr-expected)<0.008,"velocity autocorrelation");
    std::cout<<"PASS exact transition, limits, settling, variance and correlation\n"
             <<"samples="<<N<<" mean="<<mean<<" variance="<<variance
             <<" correlation="<<corr<<" expected="<<expected<<'\n';
}
