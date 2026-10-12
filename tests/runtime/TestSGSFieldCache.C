/* Native model cache oracle with irregular carrier time steps.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "basicKinematicCloud.H"
#include <vector>
#include <cmath>

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    volVectorField U(IOobject("U", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    #include "createPhi.H"
    singlePhaseTransportModel transport(U, phi);
    dimensionedScalar density("rhoInf", dimDensity, transport);
    volScalarField rho(IOobject("rho", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE), mesh, density);
    volScalarField muc(IOobject("muc", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE), rho*transport.nu());
    autoPtr<incompressible::turbulenceModel> turbulence
        (incompressible::turbulenceModel::New(U, phi, transport));
    basicKinematicCloud cloud("kinematicCloud", rho, U, muc, g);
    turbulence->validate();

    const scalar nut = turbulence->nut()()[0];
    const scalar window = 0.5, delta = 10, Ce = 0.93;
    struct Snapshot { double duration, shear; };
    std::vector<Snapshot> history;
    scalar lastFieldTime = 0, previousVariance = 0, oldShear = 0.1;
    bool havePrevious = false;
    scalar maxError = 0;
    auto check = [&](scalar got, scalar expected)
    {
        const scalar error = std::abs(got-expected)/(1+std::abs(expected));
        maxError = std::max(maxError, error);
        if (!std::isfinite(error) || error > 2e-10)
        {
            FatalErrorInFunction << "Cache mismatch: " << got << " expected "
                << expected << " at " << runTime.value() << exit(FatalError);
        }
    };
    for (scalar dt : {0.2, 0.1, 0.35, 0.05, 0.5, 0.1, 0.025})
    {
        runTime.setDeltaT(dt);
        ++runTime;
        const scalar fieldTime = runTime.value()-dt;
        const scalar sampleDt = fieldTime-lastFieldTime;
        if (sampleDt > SMALL)
        {
            history.push_back({sampleDt, oldShear});
            lastFieldTime = fieldTime;
        }
        // Independent brute-force window: trim durations, then recompute sums.
        scalar duration = 0;
        for (const auto& item : history) duration += item.duration;
        while (duration > window+1e-14 && !history.empty())
        {
            const scalar remove = std::min(duration-window, history.front().duration);
            history.front().duration -= remove;
            duration -= remove;
            if (history.front().duration < 1e-14) history.erase(history.begin());
        }
        scalar meanG = oldShear, meanG2 = sqr(oldShear), varG = 0;
        if (duration > SMALL)
        {
            meanG = meanG2 = 0;
            for (const auto& item : history)
            {
                meanG += item.duration*item.shear/duration;
                meanG2 += item.duration*sqr(item.shear)/duration;
            }
            for (const auto& item : history)
                varG += item.duration*sqr(item.shear-meanG)/duration;
        }
        const scalar epsilon = nut*meanG2;
        const scalar energy = std::pow(epsilon*delta/Ce, 2.0/3.0);
        const scalar variance = 2*energy/3;
        const scalar varianceDt = havePrevious && sampleDt > SMALL
            ? (variance-previousVariance)/sampleDt : 0;
        cloud.dispersion().cacheFields(true);
        auto lookup = [&](const word& suffix) -> const volScalarField&
        { return mesh.lookupObject<volScalarField>("kinematicCloudSGS"+suffix); };
        forAll(mesh.C(), celli)
        {
            const scalar kResolved = 0.5*sqr(mesh.C()[celli].z())*varG;
            check(lookup("Dissipation")[celli], epsilon);
            check(lookup("Variance")[celli], variance);
            check(lookup("VarianceDt")[celli], varianceDt);
            check(lookup("ResolvedK")[celli], kResolved);
            check(lookup("Fraction")[celli], energy/(energy+kResolved));
            check(lookup("WindowDuration")[celli], duration);
        }
        const scalar saved = lookup("VarianceDt")[0];
        cloud.dispersion().cacheFields(true);
        check(lookup("VarianceDt")[0], saved);
        Info<< "Cache oracle: carrier time=" << fieldTime << " current dt="
            << dt << " field interval=" << sampleDt << " duration=" << duration
            << " dVariance/dt=" << saved << nl;
        previousVariance = variance;
        havePrevious = true;
        const scalar newShear = 0.1+0.02*runTime.value();
        forAll(mesh.C(), celli)
            U[celli] = vector(newShear*mesh.C()[celli].z(), 0, 0);
        forAll(U.boundaryField(), patchi)
        {
            vectorField values(mesh.Cf().boundaryField()[patchi].size(), Zero);
            forAll(values, facei)
                values[facei].x() = newShear*mesh.Cf().boundaryField()[patchi][facei].z();
            U.boundaryFieldRef()[patchi] == values;
        }
        U.correctBoundaryConditions();
        oldShear = newShear;
    }
    Info<< "PASS native irregular-time cache, production, weighted variance, "
        << "partial expiry, derivative and duplicate-call checks; max error="
        << maxError << nl;
    return 0;
}
