#ifndef DENSITY_TRACER_RAY_H_GUARD
#define DENSITY_TRACER_RAY_H_GUARD

#include <memory>
#include <utility>
#include <vector>

#include "density-tracer/quaternion.h"

const real FREQ_INFRARED = 6;
const real FREQ_NEAR_INFRARED = 7;
const real FREQ_RED = 8;
const real FREQ_YELLOW = 9;
const real FREQ_GREEN = 10;
const real FREQ_CYAN = 11;
const real FREQ_BLUE = 12;
const real FREQ_VIOLET = 13;
const real FREQ_NEAR_ULTRAVIOLET = 14;
const real FREQ_ULTRAVIOLET = 15;

// Linear approximation from black to red to green to blue to magenta to black
// This is only physically approximate. See the reference if you want to get fancy.
// https://en.wikipedia.org/wiki/CIE_1931_color_space
color frequency_to_rgb(const real& frequency);

class Pigment {
public:
    virtual real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& frequency) const = 0;
};

class Density {
public:
    virtual std::pair<real, real> eval(const quaternion& location, const quaternion& direction, const real& frequency) const = 0;
};

class SkySphere {
public:
    virtual real eval(const quaternion& direction, const real& frequency) const = 0;
};

class Traceable {
public:
    quaternion location = Q_ZERO;
    quaternion left_transform = Q_ONE;
    quaternion right_transform = Q_ONE;
    const Pigment *pigment;
    const Pigment *reflectivity;
    virtual std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const = 0;
};

class Plane : public Traceable {
public:
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class Tetrahedron : public Traceable {
public:
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class Ball : public Traceable {
public:
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class RayPath {
    real frequency;
    quaternion start;
    quaternion direction;
    real length;
    real total_length;
    real end_amplitude;
    real end_alpha;
    RayPath *child;
public:
    RayPath(const quaternion& origin, const quaternion& direction_, const std::vector<std::shared_ptr<Traceable>>& objects, const SkySphere& sky_sphere, const real& max_length, const int& max_depth, const real& frequency_);

    real eval(const Density& density, const int& num_samples) const;
};
#endif
