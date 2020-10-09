#ifndef DENSITY_TRACER_TRACEABLE_H_GUARD
#define DENSITY_TRACER_TRACEABLE_H_GUARD

#include <vector>

#include "density-tracer/ray.h"

class Transformation : public Traceable {
public:
    std::shared_ptr<Traceable> object;
    Transformation(std::shared_ptr<Traceable> obj) : object(obj) {}
};

class Translate : public Transformation {
public:
    quaternion location = Q_ZERO;
    Translate(std::shared_ptr<Traceable> obj, quaternion loc) : Transformation(obj), location(loc) {}
    bool inside(const quaternion& location) const;
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class Rotate : public Transformation {
    quaternion left_rotor = Q_ONE;
    quaternion right_rotor = Q_ONE;
public:
    Rotate(std::shared_ptr<Traceable> obj, quaternion axis, real amount);
    bool inside(const quaternion& location) const;
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class Scale : public Transformation {
public:
    quaternion amount = {0, 1, 1, 1};
    Scale(std::shared_ptr<Traceable> obj, quaternion a) : Transformation(obj), amount(a) {}
    bool inside(const quaternion& location) const;
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class CSG : public Traceable {
public:
    std::shared_ptr<Traceable> first;
    std::shared_ptr<Traceable> second;
    CSG(std::shared_ptr<Traceable> a, std::shared_ptr<Traceable> b) : first(a), second(b) {}
};

class Merge : public CSG {
public:
    Merge(std::shared_ptr<Traceable> a, std::shared_ptr<Traceable> b) : CSG(a, b) {}
    bool inside(const quaternion& location) const;
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class Plane : public Traceable {
public:
    bool inside(const quaternion& location) const;
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class Tetrahedron : public Traceable {
public:
    bool inside(const quaternion& location) const;
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

class ConvexPolyhedron : public Traceable {
public:
    std::vector<quaternion> planes;
    ConvexPolyhedron(const std::vector<quaternion> planes);

    bool inside(const quaternion& location) const;
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

    bool inside(const quaternion& location) const;
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
    bool inside(const quaternion& location) const;
    std::pair<real, quaternion> trace(const quaternion& origin, const quaternion& direction) const;
};

#endif