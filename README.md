# AabbCyl

Small command line comparison of axis aligned bounding box (aabb) - cylinder (cyl) collision detection.
The LCP solver from Geometric Tools (https://github.com/davideberly/GeometricTools) and a custom method are compared.
The custom method follows the following steps:
 - checking whether the cyl's axis is a separating axis;
 - cropping the aabb with the cyl's top and bottom planes = convex polyhedron;
 - projecting the polyhedron's vertices to a plane perpendicular to the cylinder's axis;
 - creating a convex hull based on the projected polyhedron's vertices;
 - comparing the convex hull with the cyl's projection.

# Bugs, remarks & questions

--> sebastianwouters [at] gmail [dot] com

# Copyright

Sebastian Wouters

Copyright (c) 2021

Distributed under the Boost Software License, Version 1.0.

https://www.boost.org/LICENSE_1_0.txt
