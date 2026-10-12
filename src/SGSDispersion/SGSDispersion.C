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

#include "SGSDispersion.H"
#include "turbulenceModel.H"
#include "fvcGrad.H"
#include "zeroGradientFvPatchFields.H"
#include <cmath>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::volScalarField>
Foam::SGSDispersion<CloudType>::makeScalarField
(
    const word& name,
    const dimensionSet& dims,
    const IOobject::writeOption write
) const
{
    const fvMesh& mesh = this->owner().mesh();

    return autoPtr<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                this->owner().name() + name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                write
            ),
            mesh,
            dimensionedScalar("zero", dims, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SGSDispersion<CloudType>::SGSDispersion
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionModel<CloudType>(dict, owner, typeName),
    C0_(this->coeffDict().template getOrDefault<scalar>("C0", 4)),
    beta_(this->coeffDict().template getOrDefault<scalar>("beta", 2)),
    cEpsilon_
    (
        this->coeffDict().template getOrDefault<scalar>("cEpsilon", 0.93)
    ),
    verticalDirection_
    (
        this->coeffDict().template getOrDefault<vector>
        ("verticalDirection", vector(0, 0, 1))
    ),
    stationaryInitial_
    (
        this->coeffDict().template getOrDefault<bool>
        ("stationaryInitial", false)
    ),
    statistics_
    (
        owner.mesh().nCells(),
        this->coeffDict().template get<scalar>("averagingTime")
    ),
    lastStatisticsTime_(owner.mesh().time().value()),
    lastTimeIndex_(-1),
    havePreviousVariance_(false),
    fullWindowReported_(false)
{
    if
    (
        !(C0_ > 0) || !(beta_ >= 0) || !(cEpsilon_ > 0)
     || !std::isfinite(C0_ + beta_ + cEpsilon_)
     || !(mag(verticalDirection_) > SMALL)
     || !std::isfinite(mag(verticalDirection_))
    )
    {
        FatalIOErrorInFunction(this->coeffDict())
            << "Invalid SGS model coefficients" << exit(FatalIOError);
    }

    verticalDirection_ /= mag(verticalDirection_);

    const fvMesh& mesh = owner.mesh();
    const turbulenceModel& turbulence = mesh.lookupObject<turbulenceModel>
    (
        IOobject::groupName(turbulenceModel::propertiesName, owner.U().group())
    );

    if (turbulence.get<word>("simulationType") != "LES")
    {
        FatalErrorInFunction
            << "sgsDispersion requires simulationType LES"
            << exit(FatalError);
    }

    sigma2_ = makeScalarField
    (
        "SGSVariance", sqr(dimVelocity), IOobject::AUTO_WRITE
    );
    alpha_ = makeScalarField
    (
        "SGSFraction", dimless, IOobject::AUTO_WRITE
    );
    epsilon_ = makeScalarField
    (
        "SGSDissipation", sqr(dimVelocity)/dimTime, IOobject::AUTO_WRITE
    );
    resolvedK_ = makeScalarField
    (
        "SGSResolvedK", sqr(dimVelocity), IOobject::AUTO_WRITE
    );
    windowDuration_ = makeScalarField
    (
        "SGSWindowDuration", dimTime, IOobject::AUTO_WRITE
    );
    dSigma2_ = makeScalarField
    (
        "SGSVarianceDt", sqr(dimVelocity)/dimTime, IOobject::NO_WRITE
    );

    Info<< "sgsDispersion: averagingTime="
        << this->coeffDict().template get<scalar>("averagingTime")
        << " s." << nl;
}


template<class CloudType>
Foam::SGSDispersion<CloudType>::SGSDispersion
(
    const SGSDispersion<CloudType>& rhs
)
:
    DispersionModel<CloudType>(rhs),
    C0_(rhs.C0_),
    beta_(rhs.beta_),
    cEpsilon_(rhs.cEpsilon_),
    verticalDirection_(rhs.verticalDirection_),
    stationaryInitial_(rhs.stationaryInitial_),
    statistics_(rhs.statistics_),
    lastStatisticsTime_(rhs.lastStatisticsTime_),
    lastTimeIndex_(rhs.lastTimeIndex_),
    havePreviousVariance_(rhs.havePreviousVariance_),
    fullWindowReported_(rhs.fullWindowReported_),
    previousVariance_(rhs.previousVariance_)
{
    // Deep-copy statistics for cloud snapshots; copied fields are unregistered.
    if (rhs.sigma2_)
    {
        sigma2_.reset(new volScalarField(*rhs.sigma2_));
    }
    if (rhs.alpha_)
    {
        alpha_.reset(new volScalarField(*rhs.alpha_));
    }
    if (rhs.epsilon_)
    {
        epsilon_.reset(new volScalarField(*rhs.epsilon_));
    }
    if (rhs.resolvedK_)
    {
        resolvedK_.reset(new volScalarField(*rhs.resolvedK_));
    }
    if (rhs.windowDuration_)
    {
        windowDuration_.reset(new volScalarField(*rhs.windowDuration_));
    }
    if (rhs.dSigma2_)
    {
        dSigma2_.reset(new volScalarField(*rhs.dSigma2_));
    }
    if (rhs.gradSigma2_)
    {
        gradSigma2_.reset(new volVectorField(*rhs.gradSigma2_));
    }
}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SGSDispersion<CloudType>::~SGSDispersion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SGSDispersion<CloudType>::cacheFields(const bool store)
{
    if (!store)
    {
        return;
    }

    const fvMesh& mesh = this->owner().mesh();
    const Time& time = mesh.time();

    if (lastTimeIndex_ == time.timeIndex())
    {
        return;
    }
    lastTimeIndex_ = time.timeIndex();

    const turbulenceModel& turbulence = mesh.lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName, this->owner().U().group()
        )
    );

    const volScalarField& delta = mesh.lookupObject<volScalarField>
    (
        IOobject::groupName("delta", this->owner().U().group())
    );

    const tmp<volScalarField> production =
        2*turbulence.nut()*magSqr(symm(fvc::grad(this->owner().U())));

    // kinematicParcelFoam evolves the cloud before solving the carrier step.
    // The available U/nut therefore represent the previously completed time.
    const scalar fieldTime = time.value() - time.deltaTValue();
    const scalar sampleDt = fieldTime - lastStatisticsTime_;

    if (sampleDt > SMALL)
    {
        statistics_.add
        (
            this->owner().U().primitiveField(), production().primitiveField(),
            sampleDt
        );
        lastStatisticsTime_ = fieldTime;
    }

    const scalar window =
        this->coeffDict().template get<scalar>("averagingTime");

    if (!fullWindowReported_ && statistics_.duration() >= window - SMALL)
    {
        Info<< "sgsDispersion: internal window filled; duration="
            << statistics_.duration() << " s; samples="
            << statistics_.sampleCount() << nl;
        fullWindowReported_ = true;
    }

    forAll(mesh.V(), celli)
    {
        const scalar eps = max
        (
            statistics_.duration() > VSMALL
          ? statistics_.meanProduction()[celli]
          : production()[celli],
            scalar(0)
        );
        const scalar e = pow(eps*delta[celli]/cEpsilon_, 2.0/3.0);
        const scalar k = statistics_.resolvedK(celli);

        epsilon_()[celli] = eps;
        resolvedK_()[celli] = k;
        windowDuration_()[celli] = statistics_.duration();
        sigma2_()[celli] = 2*e/3;
        alpha_()[celli] = e + k > VSMALL ? e/(e + k) : 0;
        dSigma2_()[celli] = havePreviousVariance_
          ? (sampleDt > SMALL
            ? (sigma2_()[celli] - previousVariance_[celli])/sampleDt
            : scalar(0))
          : 0;
    }

    sigma2_().correctBoundaryConditions();
    alpha_().correctBoundaryConditions();
    epsilon_().correctBoundaryConditions();
    resolvedK_().correctBoundaryConditions();
    windowDuration_().correctBoundaryConditions();
    gradSigma2_.reset(new volVectorField(fvc::grad(sigma2_())));
    previousVariance_ = sigma2_().primitiveField();
    havePreviousVariance_ = true;
}


template<class CloudType>
Foam::vector Foam::SGSDispersion<CloudType>::update
(
    const scalar dt,
    const label celli,
    const vector& U,
    const vector& Uc,
    vector& UTurb,
    scalar& tTurb
)
{
    if (dt <= 0)
    {
        return Uc + UTurb;
    }

    const scalar variance = sigma2_()[celli];
    const scalar eps = epsilon_()[celli];
    const scalar alpha = alpha_()[celli];

    if (variance <= VSMALL || eps <= VSMALL || alpha <= VSMALL)
    {
        UTurb = Zero;
        tTurb += dt;
        return Uc;
    }

    Random& random = this->owner().rndGen();

    if (tTurb == 0 && stationaryInitial_)
    {
        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            UTurb[cmpt] = sqrt(variance)*random.template GaussNormal<scalar>();
        }
    }

    // Comola et al. (2019), Eq. (8): absolute vertical particle velocity.
    const vector gradient = gradSigma2_()[celli];
    const auto local = sgsDispersionDetail::coefficients
    (
        variance,
        eps,
        alpha,
        U,
        U & verticalDirection_,
        gradient,
        dSigma2_()[celli],
        C0_,
        beta_
    );
    const scalar count = std::ceil(max(mag(local.a), alpha/local.Tp)*dt/0.25);

    if (!std::isfinite(count) || count > 10000)
    {
        FatalErrorInFunction
            << "Invalid or excessive SGS drift; reduce time step"
            << " or inspect variance"
            << exit(FatalError);
    }

    const label nSubSteps = max(label(1), label(count));

    // Exact linear transition with coefficients frozen within this parcel call.
    for (label step = 0; step < nSubSteps; ++step)
    {
        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            const auto transition = sgsDispersionDetail::transition
            (
                local.a, local.forcing[cmpt], local.B2, dt/nSubSteps
            );

            UTurb[cmpt] = transition.decay*UTurb[cmpt] + transition.forcing
              + transition.noiseStd*random.template GaussNormal<scalar>();
        }
    }

    tTurb += dt;
    return Uc + UTurb;
}

// ************************************************************************* //
