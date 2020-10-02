#ifndef DENSITY_TRACER_RAY_H_GUARD
#define DENSITY_TRACER_RAY_H_GUARD

#include <memory>
#include <utility>
#include <vector>

#include "density-tracer/quaternion.h"

class Traceable {
public:
    quaternion location;
    quaternion left_transform;
    quaternion right_transform;
    virtual std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const = 0;
};

class Plane : public Traceable {
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
    RayPath *child;
public:
    RayPath(const quaternion& origin, const quaternion& direction_, const std::vector<std::shared_ptr<Traceable>>& objects, const int& max_depth);

    quaternion get_location(const real& t) const;
};
#endif
