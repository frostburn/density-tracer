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

class Black : public Pigment {
public:
    real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& frequency) const {
        return 0;
    }
};

class White : public Pigment {
public:
    real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& frequency) const {
        return 1;
    }
};

class Reflectivity {
public:
    virtual real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& ior, const real& frequency) const = 0;
};

class Dull : public Reflectivity {
    real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& ior, const real& frequency) const {
        return 0;
    }
};

class PerfectMirror : public Reflectivity {
    real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& ior, const real& frequency) const {
        return 1;
    }
};

class Fresnel : public Reflectivity {
public:
    real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& ior, const real& frequency) const;
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
    const Pigment *finish;  // Additive surface highlights
    const Pigment *ior;  // Index of refraction
    const Reflectivity *reflectivity;  // Reflectivity of the surface. Can depend on IOR, usually a Fresnel instance.
    const Pigment *transparency;  // Relative transparency applied after reflections
    const Pigment *pigment;  // Base ambient pigment
    virtual bool inside(const quaternion& location) const = 0;
    virtual std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const = 0;
};

class RayPath {
    real frequency;
    quaternion start;
    quaternion direction;
    real length;
    real path_length;
    std::shared_ptr<Traceable> end_object;
    quaternion end_normal;
    real reflection_weight;
    real refraction_weight;
    RayPath *reflected_path;
    RayPath *refracted_path;
public:
    RayPath(
        const quaternion& origin,
        const quaternion& direction_,
        const std::vector<std::shared_ptr<Traceable>>& objects,
        const real& max_length,
        const int& max_depth,
        const real& frequency_);

    ~RayPath();

    real eval(const Density& density, const SkySphere& sky_sphere, const int& num_samples) const;
    quaternion end() const;
    quaternion get_direction() const;
    const RayPath* get_reflected_path() const;
    const RayPath* get_refracted_path() const;
};
#endif
