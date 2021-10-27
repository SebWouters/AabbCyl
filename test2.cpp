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

    TIQuery<Real, AlignedBox3<Real>, Cylinder3<Real>> solver1;
    Solver<Real> solver2;

    auto start = std::chrono::system_clock::now();
    size_t countLcpIntr = 0;
    for (const AlignedBox3<Real>& box : boxes)
        for (const Cylinder3<Real>& cyl : cylinders)
        {
            const auto result = solver1(box, cyl);
            countLcpIntr += result.intersect;
        }
    auto end = std::chrono::system_clock::now();
    double timeS = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * 1e-9;
    std::cout << "LCP: Number of intersects / Number of cases = " << countLcpIntr << " / " << boxes.size() * cylinders.size() << std::endl;
    std::cout << "LCP: Time [seconds] = " << timeS << std::endl;

    start = std::chrono::system_clock::now();
    size_t countProjIntr = 0;
    for (const AlignedBox3<Real>& box : boxes)
        for (const Cylinder3<Real>& cyl : cylinders)
        {
            const auto result = solver2(box, cyl);
            countProjIntr += result.intersect;
        }
    end = std::chrono::system_clock::now();
    timeS = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * 1e-9;
    std::cout << "Projection: Number of intersects / Number of cases = " << countProjIntr << " / " << boxes.size() * cylinders.size() << std::endl;
    std::cout << "Projection: Time [seconds] = " << timeS << std::endl;

    for (const AlignedBox3<Real>& box : boxes)
        for (const Cylinder3<Real>& cyl : cylinders)
        {
            const auto result1 = solver1(box, cyl);
            const auto result2 = solver2(box, cyl);
            if (result1.intersect != result2.intersect)
            {
                std::cout << "Difference\n";
                std::cout << "    box.min = { " << box.min[0] << ", " << box.min[1] << ", " << box.min[2] << " }\n";
                std::cout << "    box.max = { " << box.max[0] << ", " << box.max[1] << ", " << box.max[2] << " }\n";
                std::cout << "    cyl.axis.origin    = { " << cyl.axis.origin[0] << ", " << cyl.axis.origin[1] << ", " << cyl.axis.origin[2] << " }\n";
                std::cout << "    cyl.axis.direction = { " << cyl.axis.direction[0] << ", " << cyl.axis.direction[1] << ", " << cyl.axis.direction[2] << " }\n";
                std::cout << "    cyl.radius = " << cyl.radius << "\n";
                std::cout << "    cyl.height = " << cyl.height << "\n";
                std::cout << "LCP: intersect        = " << result1.intersect << "\n";
                std::cout << "LCP: numIterations    = " << result1.numLCPIterations << "\n";
                std::cout << "Projection: intersect = " << result2.intersect << std::endl;
                return 0;
            }
        }

    /*
    Example failure:
        Difference
            box.min = { -0.572262, -0.591058, 0.83576 }
            box.max = { -0.206206, -0.129254, 1.13998 }
            cyl.axis.origin    = { 0.448691, -0.588225, 0.176887 }
            cyl.axis.direction = { 0.224273, -0.115889, -0.967611 }
            cyl.radius = 0.155937
            cyl.height = 0.441871
        LCP: intersect        = 1
        LCP: numIterations    = 3
        Projection: intersect = 0
    However, cyl.axis.origin[2] + 0.5 * cyl.height + cyl.radius = 0.5537595 < box.min[2]
    */

    return 0;
}


