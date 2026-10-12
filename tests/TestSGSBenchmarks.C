// SPDX-License-Identifier: GPL-3.0-or-later
// Independent analytic/distribution benchmarks of the production SGS coefficient kernel.
#include "../src/SGSDispersion/SGSStochasticStep.H"
#include <array>
#include <random>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>

using V = std::array<double,3>;
void require(bool ok,const char* name)
{ if (!ok) { std::cerr<<"FAIL "<<name<<'\n'; std::exit(1); } }

// Independent direct reading of paper equation (4), not a second call to coefficients().
void coefficientsOracle()
{
    std::mt19937_64 rng(982341);
    std::uniform_real_distribution<double> positive(.01,2), sign(-1,1);
    double worst=0;
    for (int n=0;n<10000;++n)
    {
        double s=positive(rng), eps=positive(rng), alpha=.01+.98*positive(rng)/2;
        V velocity={sign(rng),sign(rng),sign(rng)}, gradient={sign(rng),sign(rng),sign(rng)};
        double partial=sign(rng), beta=positive(rng), C0=positive(rng)+1;
        auto c=sgsDispersionDetail::coefficients(s,eps,alpha,velocity,velocity[2],gradient,partial,C0,beta);
        double Tf=2*s/(C0*eps), Tp=Tf/std::sqrt(1+beta*beta*velocity[2]*velocity[2]/s);
        double derivative=partial;
        for(int i=0;i<3;++i) derivative+=velocity[i]*gradient[i];
        for(int i=0;i<3;++i)
        {
            double fluctuation=sign(rng);
            double paper=-alpha*fluctuation/Tp+.5*(derivative*fluctuation/s+gradient[i]);
            double implementation=-c.a*fluctuation+c.forcing[i];
            worst=std::max(worst,std::abs(paper-implementation)/(1+std::abs(paper)));
        }
        require(std::abs(c.B2-2*alpha*s/Tp)<1e-12,"paper noise amplitude");
    }
    require(worst<1e-13,"paper drift oracle");
    std::cout<<"PAPER_COEFFICIENTS samples=10000 max_scaled_drift_error="<<worst<<'\n';
}

// Taylor displacement variance for a stationary OU velocity, integrated over [0,T].
double taylor(double variance,double tau,double T)
{ return 2*variance*tau*(T-tau*(1-std::exp(-T/tau))); }

void homogeneous()
{
    const int N=100000;
    const double s=.4,eps=.08,alpha=.2,T=10;
    const double tau=2*s/(4*eps*alpha);
    const double reference=taylor(s,tau,T);
    double preceding=1e100;
    for(double h:{.5,.25,.125})
    {
        int steps=int(T/h);
        auto c=sgsDispersionDetail::coefficients(s,eps,alpha,V{0,0,0},0,V{0,0,0},0);
        auto q=sgsDispersionDetail::transition(c.a,0,c.B2,h);
        // Exact covariance sum of the discrete, explicit position integral.
        double discrete=s*h*h*steps;
        for(int k=1;k<steps;++k) discrete+=2*s*h*h*(steps-k)*std::exp(-k*h/tau);
        double bias=std::abs(discrete/reference-1);
        require(bias<preceding*.26,"position integral converges on halving dt");
        preceding=bias;
        std::mt19937_64 rng(44197);
        std::normal_distribution<double> normal;
        double mx=0,x2=0,u2=0;
        for(int n=0;n<N;++n)
        {
            double u=std::sqrt(s)*normal(rng),x=0;
            for(int k=0;k<steps;++k)
            { x+=u*h; u=q.decay*u+q.noiseStd*normal(rng); }
            mx+=x; x2+=x*x; u2+=u*u;
        }
        double variance=x2/N-std::pow(mx/N,2);
        require(std::abs(variance/discrete-1)<.025,"homogeneous trajectory variance vs independent covariance");
        require(std::abs(u2/N/s-1)<.02,"homogeneous velocity variance");
        std::cout<<"HOMOGENEOUS dt="<<h<<" N="<<N<<" measured_x_variance="<<variance
                 <<" discrete_reference="<<discrete<<" continuous_reference="<<reference
                 <<" position_bias="<<bias<<'\n';
    }
}

struct Mixing { double chiSquare, normalisedVelocityVariance, firstHarmonic, secondHarmonic; };
Mixing mixing(double h,bool correction)
{
    const int N=60000,bins=16;
    const double L=10,pi=std::acos(-1.),T=20,eps=.08,alpha=.25;
    std::vector<double> z(N);
    std::vector<V> u(N);
    std::mt19937_64 rng(853129);
    std::uniform_real_distribution<double> uniform(0,L);
    std::normal_distribution<double> normal;
    auto variance=[&](double x){return .4*(1+.6*std::sin(2*pi*x/L));};
    for(int n=0;n<N;++n)
    { z[n]=uniform(rng);for(int i=0;i<3;++i)u[n][i]=std::sqrt(variance(z[n]))*normal(rng); }
    for(int k=0;k<int(T/h+.5);++k)
    {
        for(int n=0;n<N;++n)
        {
            // Match the native ordering: explicit transport using old velocity,
            // followed by the SGS velocity update at the reached position.
            z[n]=std::fmod(z[n]+u[n][2]*h,L);
            if(z[n]<0)z[n]+=L;
            const double s=variance(z[n]);
            V gradient={0,0,correction?.4*.6*(2*pi/L)*std::cos(2*pi*z[n]/L):0};
            // Fluid-tracer consistency test of Eq.(4): beta=0, hence Tp=Tf.
            // Heavy particles are not required to remain uniformly distributed.
            auto c=sgsDispersionDetail::coefficients(s,eps,alpha,u[n],u[n][2],gradient,0,4,0);
            for(int i=0;i<3;++i)
            {
                auto q=sgsDispersionDetail::transition(c.a,c.forcing[i],c.B2,h);
                u[n][i]=q.decay*u[n][i]+q.forcing+q.noiseStd*normal(rng);
            }
        }
    }
    std::array<int,bins> counts={};
    double v2=0,s1=0,s2=0;
    for(int n=0;n<N;++n)
    {
        ++counts[std::min(bins-1,int(z[n]/L*bins))];
        for(int i=0;i<3;++i)v2+=u[n][i]*u[n][i]/variance(z[n]);
        s1+=std::sin(2*pi*z[n]/L);s2+=std::cos(2*pi*z[n]/L);
    }
    double chi=0,expected=double(N)/bins;
    for(int count:counts)chi+=(count-expected)*(count-expected)/expected;
    return {chi,v2/(3*N),s1/N,s2/N};
}

void nonuniform()
{
    // For uniform spatial density and Gaussian conditional velocity of variance s(z),
    // the streaming and gradient drift parts of the stationary Fokker-Planck
    // equation cancel analytically. The OU part separately leaves that PDF invariant.
    const double s=.37,g=.14,w=-.8,r2=1.1;
    double streaming=-w*g/(2*s)*(r2/s-3);
    double velocityDivergence=-g*w/(2*s)*(3-r2/s);
    require(std::abs(streaming+velocityDivergence)<1e-14,"inhomogeneous stationary PDF residual");
    for(double h:{.04,.02,.01})
    {
        auto m=mixing(h,true);
        std::cout<<"WELL_MIXED dt="<<h<<" N=60000 bins=16 chi_square="<<m.chiSquare
                 <<" normalised_velocity_variance="<<m.normalisedVelocityVariance
                 <<" harmonic_sin="<<m.firstHarmonic<<" harmonic_cos="<<m.secondHarmonic<<'\n';
        // Chi-square(15) 99.9% quantile ~37.7; 40 is a conservative fixed threshold.
        require(m.chiSquare<40,"nonuniform fluid-tracer spatial well-mixed benchmark");
        require(std::abs(m.normalisedVelocityVariance-1)<.035,"nonuniform conditional variance");
    }
    auto control=mixing(.02,false);
    std::cout<<"NEGATIVE_CONTROL dt=.02 chi_square="<<control.chiSquare
             <<" normalised_velocity_variance="<<control.normalisedVelocityVariance<<'\n';
    require(control.chiSquare>100,"well-mixed test detects omitted gradient drift");
}

// Heavy-particle filter in homogeneous turbulence: independently known stationary
// covariance and continuous displacement variance, with constant relaxation time.
void heavyParticles()
{
    const double s=.4, eps=.08, alpha=.2, T=10, response=.5;
    const double lambda=4*eps*alpha/(2*s), q=lambda*response;
    const double referenceVelocity=s/(1+q);
    const double referencePosition=2*s/(1-q*q)*
      ((T/lambda-(1-std::exp(-lambda*T))/(lambda*lambda))
        -q*(T*response-response*response*(1-std::exp(-T/response))));
    for (double h:{.1,.05,.025})
    {
        const int N=100000, steps=int(T/h);
        auto c=sgsDispersionDetail::coefficients(s,eps,alpha,V{0,0,0},0,V{0,0,0},0,4,0);
        auto transition=sgsDispersionDetail::transition(c.a,0,c.B2,h);
        const double attenuation=std::exp(-h/response);
        std::mt19937_64 rng(518027);
        std::normal_distribution<double> normal;
        double sumX=0,sumX2=0,sumV=0,sumV2=0;
        for(int n=0;n<N;++n)
        {
            double w=std::sqrt(s)*normal(rng);
            double v=referenceVelocity/s*w
              +std::sqrt(referenceVelocity-referenceVelocity*referenceVelocity/s)*normal(rng);
            double x=0;
            for(int k=0;k<steps;++k)
            {
                x+=h*v;
                w=transition.decay*w+transition.noiseStd*normal(rng);
                v=attenuation*v+(1-attenuation)*w;
            }
            sumX+=x;sumX2+=x*x;sumV+=v;sumV2+=v*v;
        }
        double varX=sumX2/N-std::pow(sumX/N,2),varV=sumV2/N-std::pow(sumV/N,2);
        std::cout<<"HEAVY dt="<<h<<" N="<<N<<" velocity_variance="<<varV
          <<" velocity_reference="<<referenceVelocity<<" x_variance="<<varX
          <<" x_reference="<<referencePosition<<'\n';
        require(std::abs(varV/referenceVelocity-1)<.025,"heavy stationary velocity filtering");
        require(std::abs(varX/referencePosition-1)<.025,"heavy displacement vs analytic inertial filter");
    }
    for(double Re:{0.,.01,.1,1.,10.,100.,1000.})
    {
        const double implementation=sgsDispersionDetail::dragCorrection(Re);
        const double oracle=1+.15*std::exp(.687*std::log(Re==0?1:Re));
        require(Re==0?implementation==1:std::abs(implementation-oracle)<1e-12,"paper sphere drag factor");
    }
    std::cout<<"PAPER_DRAG seven Reynolds numbers checked; Stokes limit exactly one\n";
}

int main()
{
    std::cout<<std::setprecision(12);
    coefficientsOracle();homogeneous();heavyParticles();nonuniform();
    std::cout<<"PASS paper drift/drag oracles, analytic tracer/heavy dispersion, time-step bias, nonuniform well-mixed and negative control\n";
}
