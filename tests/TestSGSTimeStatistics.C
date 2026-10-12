#include "SGSTimeStatistics.H"
#include <iostream>
#include <cstdlib>

void check(bool condition, const char* name)
{
    if (!condition)
    {
        std::cerr << "FAIL " << name << '\n';
        std::exit(1);
    }
}

int main()
{
    using namespace Foam;
    SGSTimeStatistics statistics(1, 4);
    vectorField velocity(1, vector(1, 0, 0));
    scalarField production(1, 2);
    statistics.add(velocity, production, 1);
    velocity[0] = vector(3, 0, 0);
    production[0] = 4;
    statistics.add(velocity, production, 2);

    check(mag(statistics.meanU()[0] - vector(7.0/3, 0, 0)) < 1e-12, "weighted mean");
    check(mag(statistics.resolvedK(0) - 4.0/9) < 1e-12, "weighted central variance");
    check(mag(statistics.meanProduction()[0] - 10.0/3) < 1e-12, "production mean");

    velocity[0] = vector(5, 0, 0);
    production[0] = 8;
    statistics.add(velocity, production, 2);
    check(mag(statistics.duration() - 4) < 1e-12, "full interval removal");
    check(mag(statistics.meanU()[0].x() - 4) < 1e-12, "window mean");
    check(mag(statistics.resolvedK(0) - 0.5) < 1e-12, "window variance");

    velocity[0] = vector(9, 0, 0);
    production[0] = 10;
    statistics.add(velocity, production, 1);
    check(mag(statistics.meanU()[0].x() - 5.5) < 1e-12, "partial interval removal");
    check(mag(statistics.resolvedK(0) - 2.375) < 1e-12, "partial window variance");
    check(mag(statistics.meanProduction()[0] - 7.5) < 1e-12, "partial production mean");

    SGSTimeStatistics copied(statistics);
    velocity[0] = vector(7, 0, 0);
    production[0] = 3;
    statistics.add(velocity, production, 8);
    check(statistics.sampleCount() == 1, "sample exceeding entire window");
    check(statistics.resolvedK(0) < 1e-12, "constant sample zero variance");
    check(mag(copied.meanU()[0].x() - 5.5) < 1e-12, "independent deep copy");

    SGSTimeStatistics stable(1, 1);
    for (int i=0; i<10000; ++i)
    {
        velocity[0] = vector(1000000 + (i%2 ? 0.01 : -0.01), 0, 0);
        production[0] = 2;
        stable.add(velocity, production, 0.1);
    }
    check(mag(stable.resolvedK(0) - 0.00005) < 1e-10, "small variance on large mean");
    check(stable.sampleCount() <= 11, "bounded window storage");
    std::cout << "PASS weighted mean/variance, exact and partial expiry, deep copy, stable variance and bounded storage\n";
}
