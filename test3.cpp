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

    const Vector<3, Real> boxMinPoint = { 0.168216, 0.389476,   0.0822103 };
    const Vector<3, Real> boxMaxPoint = { 0.291133, 0.452832,   0.194713  };
    const Vector<3, Real> cylAxisOrig = { 0.490062, 0.516011,   0.705511  };
    const Vector<3, Real> cylAxisDir  = { 0.357778, 0.00360225, 0.9338    };
    const Real            cylRadius   = 0.130149;
    const Real            cylHeight   = 0.203217;

    const AlignedBox3<Real> box(boxMinPoint, boxMaxPoint);
    const Cylinder3<Real>   cyl(Line<3, Real>(cylAxisOrig, cylAxisDir), cylRadius, cylHeight);

    const auto resultLcp  = TIQuery<Real, AlignedBox3<Real>, Cylinder3<Real>>()(box, cyl);
    const auto resultProj = Solver<Real>()(box, cyl);

    std::cout << "box.min            = { " << box.min[0] << ", " << box.min[1] << ", " << box.min[2] << " }" << std::endl;
    std::cout << "box.max            = { " << box.max[0] << ", " << box.max[1] << ", " << box.max[2] << " }" << std::endl;
    std::cout << "cyl.axis.origin    = { " << cyl.axis.origin[0]    << ", " << cyl.axis.origin[1]    << ", " << cyl.axis.origin[2]    << " }" << std::endl;
    std::cout << "cyl.axis.direction = { " << cyl.axis.direction[0] << ", " << cyl.axis.direction[1] << ", " << cyl.axis.direction[2] << " }" << std::endl;
    std::cout << "cyl.radius         = " << cyl.radius << std::endl;
    std::cout << "cyl.height         = " << cyl.height << std::endl;

    std::cout << "LCP:  intersect    = " << resultLcp.intersect  << std::endl;
    std::cout << "Proj: intersect    = " << resultProj.intersect << std::endl;
    const Real term1 = box.max[2];
    const Real term2 = cyl.axis.origin[2] - static_cast<Real>(0.5) * cyl.height - cyl.radius;
    std::cout << "box:  max[2]                                 = " << term1 << std::endl;
    std::cout << "cyl:  axis.origin[2] - 0.5 * height - radius = " << term2 << std::endl;

    return 0;
}


