// Sebastian Wouters
// Copyright (c) 2021
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Alternative AABB-CYL collision detection

#include <Mathematics/IntrAlignedBox3Cylinder3.h>
#include <Mathematics/Line.h>

#include <iostream>
#include <string>
#include <chrono>
#include <math.h>

namespace gte
{
    template <typename Real>
    class Solver
    {
    public:

        enum Result : uint8_t
        {
            error = 0,
            intersects = 1,
            disjoint = 2
        };

        std::string str(const Result result) const
        {
            switch (result)
            {
            case Result::intersects:
                return "intersects";
                break;
            case Result::disjoint:
                return "disjoint";
                break;
            case Result::error:
            default:
                return "error";
            }
        }

    private:

        template <size_t size>
        uint8_t convexHull2dSimple(std::array<uint8_t, size>& hull, const uint8_t numPoints, const std::array<Vector<2, Real>, size>& points) const
        {
            if (numPoints <= 2)
            {
                for (uint8_t idx = 0; idx < numPoints; ++idx) { hull[idx] = idx; }
                return numPoints;
            }

            // Compute the convex hull counter clockwise
            // Use Andrew's algorithm (see section 3.9 of Christer, Real-time collision detection (2005), ISBN 1-55860-732-3)
            std::array<uint8_t, size> sorted = {};
            for (uint8_t idx = 0; idx < numPoints; ++idx) { sorted[idx] = idx; }
            std::sort(sorted.begin(), sorted.begin() + numPoints,
                [&points](const uint8_t left, const uint8_t right)
                {
                    return (points[left][0] <  points[right][0]) ||
                          ((points[left][0] == points[right][0]) && (points[left][1] < points[right][1]));
                });

            const uint8_t jMin = sorted[0];
            const uint8_t jMax = sorted[numPoints - 1];

            // Lower half: points below line(points[jMin], points[jMax])
            hull[0] = jMin;
            uint8_t last = 0;
            for (uint8_t next = 1; next < numPoints - 1; ++next)
            {
                const uint8_t jNext = sorted[next];
                const uint8_t jLast = hull[last];
                Vector<2, Real> v1 = points[jMax]  - points[jLast];
                Vector<2, Real> v2 = points[jNext] - points[jLast];
                // if points[jNext] below line(points[jLast], points[jMax]): insert jNext in hull
                if (v1[0] * v2[1] < v1[1] * v2[0])
                {
                    // points[jMin] must belong to hull: do not check erasing jMin
                    bool check = last != 0;
                    while (check)
                    {
                        const uint8_t j1 = hull[last];
                        const uint8_t j0 = hull[last - 1];
                        v1 = points[j1]    - points[j0];
                        v2 = points[jNext] - points[j0];
                        // if points[jNext] below or on line(points[j0], points[j1]): erase j1 from hull
                        check = (v1[0] * v2[1] <= v1[1] * v2[0]) && (--last != 0);
                    }
                    hull[++last] = jNext;
                    // if points[jNext] below line(points[jMin], points[jMax]), skip check whether above:
                    sorted[next] = UINT8_MAX;
                }
            }

            // Upper half: points above line(points[jMin], points[jMax])
            hull[++last] = jMax;
            const uint8_t limit = last;
            for (uint8_t next = numPoints - 2; next > 0; --next)
            {
                const uint8_t jNext = sorted[next];
                if (jNext != UINT8_MAX)
                {
                    const uint8_t jLast = hull[last];
                    Vector<2, Real> v1 = points[jMin]  - points[jLast];
                    Vector<2, Real> v2 = points[jNext] - points[jLast];
                    // if points[jNext] above line(points[jLast], points[jMin]): insert jNext in hull
                    if (v1[0] * v2[1] < v1[1] * v2[0])
                    {
                        // points[jMax] must belong to hull: do not check erasing jMax
                        bool check = last != limit;
                        while (check)
                        {
                            const uint8_t j1 = hull[last];
                            const uint8_t j0 = hull[last - 1];
                            v1 = points[j1]    - points[j0];
                            v2 = points[jNext] - points[j0];
                            // if points[jNext] above or on line(points[j0], points[j1]): erase j1 from hull
                            check = (v1[0] * v2[1] <= v1[1] * v2[0]) && (--last != limit);
                        }
                        hull[++last] = jNext;
                    }
                }
            }

            return ++last;
        }

    public:

        Result operator()(AlignedBox3<Real> const& box, Cylinder3<Real> const& cyl) const
        {
            constexpr Real tolerance = static_cast<Real>(1e-10);

            struct Vertex
            {
                const Vector<3, Real> coord;
                const Real            projAxis;
                Vertex(const std::array<Real, 3>& in, const Vector<3, Real>& axis)
                    : coord(in)
                    , projAxis(Dot(coord, axis)) {}
            };

            /*
                z
                ^  6-----7
                | /|    /|
                |/ |   / |
                4--+--5  |
                |  2--+--3
                | /   | /
                |/    |/
                0-----1----> x
            */

            const std::array<Vertex, 8> vertices = {
                Vertex({ box.min[0], box.min[1], box.min[2] }, cyl.axis.direction),
                Vertex({ box.max[0], box.min[1], box.min[2] }, cyl.axis.direction),
                Vertex({ box.min[0], box.max[1], box.min[2] }, cyl.axis.direction),
                Vertex({ box.max[0], box.max[1], box.min[2] }, cyl.axis.direction),
                Vertex({ box.min[0], box.min[1], box.max[2] }, cyl.axis.direction),
                Vertex({ box.max[0], box.min[1], box.max[2] }, cyl.axis.direction),
                Vertex({ box.min[0], box.max[1], box.max[2] }, cyl.axis.direction),
                Vertex({ box.max[0], box.max[1], box.max[2] }, cyl.axis.direction) };

            const Real midAxisCyl = Dot(cyl.axis.direction, cyl.axis.origin);
            const Real halfHeight = static_cast<Real>(0.5) * cyl.height;
            const Real maxAxisCyl = midAxisCyl + halfHeight;
            const Real minAxisCyl = midAxisCyl - halfHeight;

            constexpr uint8_t binInside = 0;
            constexpr uint8_t binAbove = 1;
            constexpr uint8_t binBelow = 2;
            std::array<uint8_t, 3> histogram = { 0, 0, 0 };
            for (const Vertex& vtx : vertices)
                ++histogram[vtx.projAxis > maxAxisCyl ? binAbove : (vtx.projAxis < minAxisCyl ? binBelow : binInside)];
            if (histogram[binAbove] == vertices.size() ||
                histogram[binBelow] == vertices.size())
                return Result::disjoint; // cyl.axis.direction is a separating axis

            // Compute the vertices of a polyhedron which is the box capped by the cylinder's
            // top and bottom planes and project them into a plane perpendicular to the
            // cylinder axis direction.
            const Vector<3, Real> direction1 = fabs(cyl.axis.direction[0]) > 0.9
                ? UnitCross({ 0.0, 1.0, 0.0 }, cyl.axis.direction)
                : UnitCross({ 1.0, 0.0, 0.0 }, cyl.axis.direction);
            const Vector<3, Real> direction2 = Cross(cyl.axis.direction, direction1);

            std::array<Vector<2, Real>, 24> points;
            int numPoints = 0;

            for (uint32_t idx = 0; idx < vertices.size(); ++idx)
            {
                const Vertex& vtx = vertices[idx];
                if (vtx.projAxis > maxAxisCyl)
                {
                    for (uint32_t shift = 0; shift < 3; ++shift)
                    {
                        const Vertex& neighbor = vertices[idx ^ (1U << shift)];
                        if (neighbor.projAxis > maxAxisCyl || fabs(vtx.projAxis - neighbor.projAxis) < tolerance)
                            continue; // both cropped by maxAxisCyl
                        const Real a = (maxAxisCyl - neighbor.projAxis) / (vtx.projAxis - neighbor.projAxis);
                        const Vector<3, Real> point = a * vtx.coord + (static_cast<Real>(1.0) - a) * neighbor.coord;
                        points[numPoints++] = { Dot(direction1, point), Dot(direction2, point) };
                    }
                }
                else if (vtx.projAxis < minAxisCyl)
                {
                    for (uint32_t shift = 0; shift < 3; ++shift)
                    {
                        const Vertex& neighbor = vertices[idx ^ (1U << shift)];
                        if (neighbor.projAxis < minAxisCyl || fabs(vtx.projAxis - neighbor.projAxis) < tolerance)
                            continue; // both cropped by minAxisCyl
                        const Real a = (minAxisCyl - neighbor.projAxis) / (vtx.projAxis - neighbor.projAxis);
                        const Vector<3, Real> point = a * vtx.coord + (static_cast<Real>(1.0) - a) * neighbor.coord;
                        points[numPoints++] = { Dot(direction1, point), Dot(direction2, point) };
                    }
                }
                else
                {
                    points[numPoints++] = { Dot(direction1, vtx.coord), Dot(direction2, vtx.coord) };
                }
            }

            // The polyhedron is convex. Its projection is a convex polygon.
            // The convex polygon can be computed as the convex hull of the projected polyhedron's vertices.
            // Determine whether the polygon contains a point sufficiently close to the cylinder's center.
            const Vector<2, Real> center = { Dot(direction1, cyl.axis.origin), Dot(direction2, cyl.axis.origin) };
            const Real radiusSquared = cyl.radius * cyl.radius;
            if (numPoints == 0)
                return Result::disjoint;
            else if (numPoints == 1)
            {
                Vector<2, Real> diff = points[0] - center;
                return Dot(diff, diff) > radiusSquared ? Result::disjoint : Result::intersects;
            }

            std::array<uint8_t, points.size()> hull;
            const uint8_t end = convexHull2dSimple(hull, numPoints, points);

            // Compare center's distance w.r.t. hull and/or whether center is enclosed in the hull
            bool enclosed = true;
            for (uint8_t idx = 0; idx < end; ++idx)
            {
                const Vector<2, Real>& P = points[hull[idx]];
                const Vector<2, Real>& Q = points[hull[(idx + 1) % end]];
                const Vector<2, Real> PQ = Q - P;
                const Vector<2, Real> PC = center - P;
                if (PQ[0] * PC[1] < PQ[1] * PC[0])
                    enclosed = false;
                // point(z) = z * Q + (1 - z) * P
                // |point(z) - C|^2 = |z * PQ - PC|^2 = z*z*PQ.PQ - 2*z*PQ.PC + PC.PC
                const Real PQsq = Dot(PQ, PQ);
                const Real PCsq = Dot(PC, PC);
                Real value = static_cast<Real>(0.0);
                if (PQsq < tolerance * tolerance)
                    value = PCsq;
                else
                {
                    const Real PQdotPC = Dot(PQ, PC);
                    const Real z = PQdotPC / PQsq;
                    if (z <= static_cast<Real>(0.0))
                        value = PCsq;
                    else if (z >= static_cast<Real>(1.0))
                    {
                        const Vector<2, Real> QC = center - Q;
                        value = Dot(QC, QC);
                    }
                    else
                        value = z * (z * PQsq - static_cast<Real>(2.0) * PQdotPC) + PCsq;
                }
                if (value <= radiusSquared)
                    return Result::intersects;
            }
            if (enclosed)
                return Result::intersects;
            else
                return Result::disjoint;
        }
    };
}


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
    std::cout << "result1.intersect = " << result1.intersect << (result1.intersect == 0 ? " (ok)" : " (not ok)") << std::endl;
    std::cout << "result1.numLCPIterations = " << result1.numLCPIterations << std::endl;
    std::cout << "result2.intersect = " << result2.intersect << (result2.intersect == 1 ? " (ok)" : " (not ok)") << std::endl;
    std::cout << "result2.numLCPIterations = " << result2.numLCPIterations << std::endl;
    std::cout << "result3.intersect = " << result3.intersect << (result3.intersect == 0 ? " (ok)" : " (not ok)") << std::endl;
    std::cout << "result3.numLCPIterations = " << result3.numLCPIterations << std::endl;
    std::cout << "result4.intersect = " << result4.intersect << (result4.intersect == 1 ? " (ok)" : " (not ok)") << std::endl;
    std::cout << "result4.numLCPIterations = " << result4.numLCPIterations << std::endl;
    std::cout << "result5.intersect = " << result5.intersect << (result5.intersect == 0 ? " (ok)" : " (not ok)") << std::endl;
    std::cout << "result5.numLCPIterations = " << result5.numLCPIterations << std::endl;
    std::cout << "result6.intersect = " << result6.intersect << (result6.intersect == 1 ? " (ok)" : " (not ok)") << std::endl;
    std::cout << "result6.numLCPIterations = " << result6.numLCPIterations << std::endl;
    std::cout << "time [seconds] = " << timeS << std::endl;

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
    std::cout << "c1 = " << solver2.str(c1) << std::endl;
    std::cout << "c2 = " << solver2.str(c2) << std::endl;
    std::cout << "c3 = " << solver2.str(c3) << std::endl;
    std::cout << "c4 = " << solver2.str(c4) << std::endl;
    std::cout << "c5 = " << solver2.str(c5) << std::endl;
    std::cout << "c6 = " << solver2.str(c6) << std::endl;
    std::cout << "time [seconds] = " << timeS << std::endl;

    return 0;
}

