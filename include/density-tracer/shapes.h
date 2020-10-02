#ifndef DENSITY_TRACER_SHAPES_H_GUARD
#define DENSITY_TRACER_SHAPES_H_GUARD

#include "density-tracer/quaternion.h"

real tetrahedron(const quaternion& loc);

real merkaba(const quaternion& loc);

real tesseract(const quaternion& loc);

real sphere(const quaternion& loc);

real torus(const quaternion& loc, const real& major_radius);

#endif
