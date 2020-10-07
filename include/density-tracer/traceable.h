#ifndef DENSITY_TRACER_TRACEABLE_H_GUARD
#define DENSITY_TRACER_TRACEABLE_H_GUARD

#include "density-tracer/ray.h"

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

#endif