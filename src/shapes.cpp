#include <algorithm>
#include <cmath>

#include "density-tracer/shapes.h"

real tetrahedron(const quaternion& loc) {
    real u = std::max(loc.x + loc.y + loc.z, loc.x - loc.y - loc.z);
    u = std::max(u, loc.y - loc.x - loc.z);
    u = std::max(u, loc.z - loc.x - loc.y);
    return u * 0.5773502691896258;
}

real merkaba(const quaternion& loc) {
    return std::min(tetrahedron(loc), tetrahedron(-loc));
}

real tesseract(const quaternion& loc) {
    return std::max(std::max(std::max(fabs(loc.w), fabs(loc.x)), fabs(loc.y)), fabs(loc.z));
}
