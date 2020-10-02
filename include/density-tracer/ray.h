#ifndef DENSITY_TRACER_RAY_H_GUARD
#define DENSITY_TRACER_RAY_H_GUARD

#include <memory>
#include <utility>
#include <vector>

#include "density-tracer/quaternion.h"

class Pigment {
public:
    virtual color eval(const quaternion& location, const quaternion& direction, const quaternion& normal) const = 0;
};

class Density {
public:
    virtual std::pair<color, color> eval(const quaternion& location, const quaternion& direction) const = 0;
};

class SkySphere {
public:
    virtual color eval(const quaternion& direction) const = 0;
};

class Traceable {
public:
    quaternion location = Q_ZERO;
    quaternion left_transform = Q_ONE;
    quaternion right_transform = Q_ONE;
    real reflectivity = 0.0;
    const Pigment *pigment;
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
    quaternion start;
    quaternion direction;
    real length;
    real total_length;
    color end_color;
    real end_alpha;
    RayPath *child;
public:
    RayPath(const quaternion& origin, const quaternion& direction_, const std::vector<std::shared_ptr<Traceable>>& objects, const SkySphere& sky_sphere, const real& max_length, const int& max_depth);

    color eval(const Density& density, const int& num_samples) const;
};
#endif
