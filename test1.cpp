// Sebastian Wouters
// Copyright (c) 2021
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Alternative AABB-CYL collision detection

#include "IntrAabb3Cyl3.h"

#include <Mathematics/IntrAlignedBox3Cylinder3.h>
#include <Mathematics/Line.h>

#include <iostream>
#include <chrono>


int main()
{
    using namespace gte;
    using Real = double;

    // Test case box1 - cyl1 : disjoint
    // Test case box1 - cyl2 : intersect
    const Vector<3, Real> box1MinPoint = { 0.0, 0.0, 0.0 };
    const Vector<3, Real> box1MaxPoint = { 1000.0, 1000.0, 1000.0 };
    const AlignedBox3<Real> box1(box1MinPoint, box1MaxPoint);

    const Vector<3, Real> point1 = { 1136.0,  952.0, 600.0 };
    const Vector<3, Real> point2 = { 1560.0, 1252.0, 900.0 };
    constexpr Real cyl1Radius = 200.0;
    constexpr Real cyl2Radius = 210.0;
    constexpr Real cyl12Height = 600.0;
    Vector<3, Real> direction12 = (point2 - point1);
    Normalize(direction12);
    Vector<3, Real> center12 = 0.5 * (point1 + point2);
    const Line<3, Real> axis12(center12, direction12);
    const Cylinder3<Real> cyl1(axis12, cyl1Radius, cyl12Height);
    const Cylinder3<Real> cyl2(axis12, cyl2Radius, cyl12Height);

    // Test case box2 - cyl3 : disjoint
    // Test case box2 - cyl4 : intersect
    const Vector<3, Real> box2MinPoint = { 0.0, 0.0, 0.0 };
    const Vector<3, Real> box2MaxPoint = { 1.0, 1.0, 1.0 };
    const AlignedBox3<Real> box2(box2MinPoint, box2MaxPoint);

    constexpr Real y = 0.25;
    const Vector<3, Real> center34 = { 1.0 + y, 1.0 + y, 1.0 };
    const Real sqrt3 = std::sqrt(3.0);
    const Vector<3, Real> direction34 = { 1.0 / sqrt3, 1.0 / sqrt3, -1.0 / sqrt3 };
    const Line<3, Real> axis34(center34, direction34);
    const Real cyl34Height = 2.0 * sqrt3;
    constexpr Real cyl3Radius = 0.20;
    constexpr Real cyl4Radius = 0.21;
    const Cylinder3<Real> cyl3(axis34, cyl3Radius, cyl34Height);
    const Cylinder3<Real> cyl4(axis34, cyl4Radius, cyl34Height);

    // Test case box2 - cyl5 : disjoint
    // Test case box2 - cyl6 : intersect
    const Vector<3, Real> center56 = { 1.5, 1.5, 1.5 };
    const Vector<3, Real> direction56 = { 1.0, 0.0, 0.0 };
    const Line<3, Real> axis56(center56, direction56);
    const Real cyl56Height = 1.001;
    const Real cyl5Radius = 0.707;
    const Real cyl6Radius = 0.708;
    const Cylinder3<Real> cyl5(axis56, cyl5Radius, cyl56Height);
    const Cylinder3<Real> cyl6(axis56, cyl6Radius, cyl56Height);

    auto start = std::chrono::system_clock::now();
    TIQuery<Real, AlignedBox3<Real>, Cylinder3<Real>> solver1;
    const auto result1 = solver1(box1, cyl1);
    const auto result2 = solver1(box1, cyl2);
    const auto result3 = solver1(box2, cyl3);
    const auto result4 = solver1(box2, cyl4);
    const auto result5 = solver1(box2, cyl5);
    const auto result6 = solver1(box2, cyl6);
    auto end = std::chrono::system_clock::now();
    double timeS = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * 1e-9;
    const bool lcpOK =
        result1.intersect == 0 &&
        result2.intersect == 1 &&
        result3.intersect == 0 &&
        result4.intersect == 1 &&
        result5.intersect == 0 &&
        result6.intersect == 1;
    std::cout << "LCP: Result = " << (lcpOK ? "ok" : "not ok") << std::endl;
    std::cout << "LCP: Time [seconds] = " << timeS << std::endl;

    start = std::chrono::system_clock::now();
    Solver<Real> solver2;
    const auto c1 = solver2(box1, cyl1);
    const auto c2 = solver2(box1, cyl2);
    const auto c3 = solver2(box2, cyl3);
    const auto c4 = solver2(box2, cyl4);
    const auto c5 = solver2(box2, cyl5);
    const auto c6 = solver2(box2, cyl6);
    end = std::chrono::system_clock::now();
    timeS = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * 1e-9;
    const bool projOK =
        c1 == Solver<Real>::Result::disjoint   &&
        c2 == Solver<Real>::Result::intersects &&
        c3 == Solver<Real>::Result::disjoint   &&
        c4 == Solver<Real>::Result::intersects &&
        c5 == Solver<Real>::Result::disjoint   &&
        c6 == Solver<Real>::Result::intersects;
    std::cout << "Proj: Result = " << (projOK ? "ok" : "not ok") << std::endl;
    std::cout << "Proj: Time [seconds] = " << timeS << std::endl;
    if (lcpOK && projOK)
        return 0;
    else
        return 127;
}


