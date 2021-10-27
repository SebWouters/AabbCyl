// Sebastian Wouters
// Copyright (c) 2021
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Alternative AABB-CYL collision detection

#include "IntrAabb3Cyl3.h"

#include <Mathematics/IntrAlignedBox3Cylinder3.h>
#include <Mathematics/Line.h>

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <math.h>

int main()
{
    using namespace gte;
    using Real = double;

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<Real> cen(-1.0000, 1.0000);
    std::uniform_real_distribution<Real> ext( 0.0001, 0.2500);
    std::uniform_real_distribution<Real> uni( 0.0000, 1.0000);

    constexpr size_t numItems = 512;
    std::vector<AlignedBox3<Real>> boxes;
    boxes.reserve(numItems);
    std::vector<Cylinder3<Real>> cylinders;
    cylinders.reserve(numItems);

    for (size_t idx = 0; idx < numItems; ++idx)
    {
        const Vector<3, Real> center = { cen(gen), cen(gen), cen(gen) };
        const Vector<3, Real> extent = { ext(gen), ext(gen), ext(gen) };
        boxes.emplace_back(center - extent, center + extent);
    }

    for (size_t idx = 0; idx < numItems; ++idx)
    {
        const Vector<3, Real> center = { cen(gen), cen(gen), cen(gen) };
        const Real phi   = 2.0 * M_PI * uni(gen);
        const Real theta = M_PI * uni(gen);
        const Real sinTheta = sin(theta);
        const Vector<3, Real> direction = { cos(phi) * sinTheta, sin(phi) * sinTheta, cos(theta) };
        const Line<3, Real> axis(center, direction);
        const Real radius = ext(gen);
        const Real height = 2.0 * ext(gen);
        cylinders.emplace_back(axis, radius, height);
    }

    auto start = std::chrono::system_clock::now();
    size_t countLcpIntr = 0;
    for (const AlignedBox3<Real>& box : boxes)
        for (const Cylinder3<Real>& cyl : cylinders)
        {
            const auto result1 = TIQuery<Real, AlignedBox3<Real>, Cylinder3<Real>>()(box, cyl);
            countLcpIntr += result1.intersect;
        }
    auto end = std::chrono::system_clock::now();
    double timeS = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * 1e-9;
    std::cout << "LCP:  Number of intersects / Number of cases = " << countLcpIntr << " / " << boxes.size() * cylinders.size() << std::endl;
    std::cout << "LCP:  Time [seconds] = " << timeS << std::endl;

    start = std::chrono::system_clock::now();
    size_t countProjIntr = 0;
    for (const AlignedBox3<Real>& box : boxes)
        for (const Cylinder3<Real>& cyl : cylinders)
        {
            const auto result2 = Solver<Real>()(box, cyl);
            countProjIntr += result2.intersect;
        }
    end = std::chrono::system_clock::now();
    timeS = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * 1e-9;
    std::cout << "Proj: Number of intersects / Number of cases = " << countProjIntr << " / " << boxes.size() * cylinders.size() << std::endl;
    std::cout << "Proj: Time [seconds] = " << timeS << std::endl;

    for (const AlignedBox3<Real>& box : boxes)
        for (const Cylinder3<Real>& cyl : cylinders)
        {
            const auto result1 = TIQuery<Real, AlignedBox3<Real>, Cylinder3<Real>>()(box, cyl);
            const auto result2 = Solver<Real>()(box, cyl);
            if (result1.intersect != result2.intersect)
            {
                std::cout << "Difference:\n";
                std::cout << "    box.min            = { " << box.min[0] << ", " << box.min[1] << ", " << box.min[2] << " }\n";
                std::cout << "    box.max            = { " << box.max[0] << ", " << box.max[1] << ", " << box.max[2] << " }\n";
                std::cout << "    cyl.axis.origin    = { " << cyl.axis.origin[0] << ", " << cyl.axis.origin[1] << ", " << cyl.axis.origin[2] << " }\n";
                std::cout << "    cyl.axis.direction = { " << cyl.axis.direction[0] << ", " << cyl.axis.direction[1] << ", " << cyl.axis.direction[2] << " }\n";
                std::cout << "    cyl.radius         = " << cyl.radius << "\n";
                std::cout << "    cyl.height         = " << cyl.height << "\n";
                std::cout << "LCP:  intersect        = " << result1.intersect << "\n";
                std::cout << "LCP:  numIterations    = " << result1.numLCPIterations << "\n";
                std::cout << "Proj: intersect        = " << result2.intersect << std::endl;
                return 0;
            }
        }

    return 0;
}


