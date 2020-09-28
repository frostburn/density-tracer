#ifndef DENSITY_TRACER_JULIA_H_GUARD
#define DENSITY_TRACER_JULIA_H_GUARD

#include <cmath>
#include <utility>

#include "density-tracer/quaternion.h"

std::pair<real, quaternion> mandelbrot(quaternion q, const quaternion& c, const quaternion& dc, const int& exponent, const int& num_iter);

#endif
