#ifndef DENSITY_TRACER_TRACEABLE_H_GUARD
#define DENSITY_TRACER_TRACEABLE_H_GUARD

#include <vector>

#include "density-tracer/ray.h"

class Plane : public Traceable {
public:
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class Tetrahedron : public Traceable {
public:
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class ConvexPolyhedron : public Traceable {
public:
    std::vector<quaternion> planes;
    ConvexPolyhedron(const std::vector<quaternion> planes);
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class TrianglePrism : public ConvexPolyhedron {
public:
    TrianglePrism();
};

class HemiConvex : public Traceable {
public:
    std::vector<quaternion> planes;
    HemiConvex(const std::vector<quaternion> planes);
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class Hexahedron : public HemiConvex {
public:
    Hexahedron();
};

class Octahedron : public HemiConvex {
public:
    Octahedron();
};

class Dodecahedron : public HemiConvex {
public:
    Dodecahedron();
};

class Icosahedron : public HemiConvex {
public:
    Icosahedron();
};

class Ball : public Traceable {
public:
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

#endif