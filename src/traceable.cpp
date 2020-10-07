#include <cmath>
#include <limits>

#include "density-tracer/traceable.h"

const real EPSILON = 1e-8;
const real ISQRT_3 = 0.5773502691896258;
const real PHI = 1.618033988749895;
const real IPHI = 0.6180339887498948;

std::pair<real, quaternion> Plane::trace(const quaternion& origin, const quaternion& direction) const {
    if ((origin.y < 0 && direction.y > EPSILON) || (origin.y > 0 && direction.y < -EPSILON)) {
        return std::make_pair(-origin.y / direction.y, Q_J);
    }
    return std::make_pair(std::numeric_limits<real>::infinity(), Q_ZERO);
}

std::pair<real, quaternion> Tetrahedron::trace(const quaternion& origin, const quaternion& direction) const {
    real o = origin.x + origin.y + origin.z;
    real d = direction.x + direction.y + direction.z;
    real t_a;
    if (fabs(d) > EPSILON) {
        t_a = (1-o) / d;
    } else {
        t_a = -1;
    }
    if (t_a > 0) {
        quaternion i = origin + direction*t_a;
        if (i.x - i.y - i.z > 1 || i.y - i.x - i.z > 1 || i.z - i.x - i.y > 1) {
            t_a = -1;
        }
    }
    o = origin.x - origin.y - origin.z;
    d = direction.x - direction.y - direction.z;
    real t_b;
    if (fabs(d) > EPSILON) {
        t_b = (1-o) / d;
    } else {
        t_b = -1;
    }
    if (t_b > 0) {
        quaternion i = origin + direction*t_b;
        if (i.x + i.y + i.z > 1 || i.y - i.x - i.z > 1 || i.z - i.x - i.y > 1) {
            t_b = -1;
        }
    }
    o = origin.y - origin.x - origin.z;
    d = direction.y - direction.x - direction.z;
    real t_c;
    if (fabs(d) > EPSILON) {
        t_c = (1-o) / d;
    } else {
        t_c = -1;
    }
    if (t_c > 0) {
        quaternion i = origin + direction*t_c;
        if (i.x + i.y + i.z > 1 || i.x - i.y - i.z > 1 || i.z - i.x - i.y > 1) {
            t_c = -1;
        }
    }
    o = origin.z - origin.y - origin.x;
    d = direction.z - direction.y - direction.x;
    real t_d;
    if (fabs(d) > EPSILON) {
        t_d = (1-o) / d;
    } else {
        t_d = -1;
    }
    if (t_d > 0) {
        quaternion i = origin + direction*t_d;
        if (i.x + i.y + i.z > 1 || i.y - i.x - i.z > 1 || i.x - i.z - i.y > 1) {
            t_d = -1;
        }
    }
    if (t_a > 0 && (t_a < t_b || t_b < 0)) {
        if (t_a < t_c || t_c < 0) {
            if (t_a < t_d || t_d < 0) {
                return std::make_pair(t_a, (quaternion){0, ISQRT_3, ISQRT_3, ISQRT_3});
            }
        }
    }
    if (t_b > 0 && (t_b < t_c || t_c < 0)) {
        if (t_b < t_d || t_d < 0) {
            return std::make_pair(t_b, (quaternion){0, ISQRT_3, -ISQRT_3, -ISQRT_3});
        }
    }
    if (t_c > 0 && (t_c < t_d || t_d < 0)) {
        return std::make_pair(t_c, (quaternion){0, -ISQRT_3, ISQRT_3, -ISQRT_3});
    }
    if (t_d > 0) {
        return std::make_pair(t_d, (quaternion){0, -ISQRT_3, -ISQRT_3, ISQRT_3});
    }
    return std::make_pair(std::numeric_limits<real>::infinity(), Q_ZERO);
}

HemiConvex::HemiConvex(const std::vector<quaternion> planes) {
    this->planes = planes;
}

std::pair<real, quaternion> HemiConvex::trace(const quaternion& origin, const quaternion& direction) const {
    real closest = std::numeric_limits<real>::infinity();
    quaternion closest_plane = Q_ONE;

    std::vector<quaternion>::const_iterator plane;
    for (plane = this->planes.begin(); plane != this->planes.end(); ++plane) {
        quaternion pplane = *plane;
        real o = dot(origin, pplane);
        real d = dot(direction, pplane);
        if (o < -1 || (o < 1 && d < 0)) {
            pplane = -pplane;
            o = -o;
            d = -d;
        }
        real t;
        if (fabs(d) > EPSILON) {
            t = (1-o) / d;
        } else {
            continue;
        }
        quaternion intersection = origin + t*direction;
        std::vector<quaternion>::const_iterator edge_plane;
        for (edge_plane = this->planes.begin(); edge_plane != this->planes.end(); ++edge_plane) {
            if (edge_plane == plane) {
                continue;
            }
            if (fabs(dot(intersection, *edge_plane)) > 1) {
                t = -1;
                break;
            };
        }
        if (t > 0 && t < closest) {
            closest = t;
            closest_plane = pplane;
        }
    }

    return std::make_pair(closest, normalize(closest_plane));
}

Hexahedron::Hexahedron() : HemiConvex({}) {
    this->planes = {Q_I, Q_J, Q_K};
}

Octahedron::Octahedron() : HemiConvex({}) {
    this->planes = {{0, 1, 1, 1}, {0, 1, 1, -1}, {0, 1, -1, 1}, {0, 1, -1, -1}};
}

Dodecahedron::Dodecahedron() : HemiConvex({}) {
    this->planes = {
        {0, 0, 1, PHI}, {0, 0, 1, -PHI},
        {0, 1, PHI, 0}, {0, 1, -PHI, 0},
        {0, PHI, 0, 1}, {0, -PHI, 0, 1},
    };
}

Icosahedron::Icosahedron() : HemiConvex({}) {
    this->planes = {
        {0, 1, 1, 1}, {0, 1, 1, -1}, {0, 1, -1, 1}, {0, 1, -1, -1},
        {0, PHI, IPHI, 0}, {0, PHI, -IPHI, 0},
        {0, 0, PHI, IPHI}, {0, 0, PHI, -IPHI},
        {0, IPHI, 0, PHI}, {0, -IPHI, 0, PHI},
    };
}

std::pair<real, quaternion> Ball::trace(const quaternion& origin, const quaternion& direction) const {
    real o2 = origin.x*origin.x + origin.y*origin.y + origin.z*origin.z;
    real od = origin.x*direction.x + origin.y*direction.y + origin.z*direction.z;
    real d2 = direction.x*direction.x + direction.y*direction.y + direction.z*direction.z;

    real discriminant = od*od - (o2-1)*d2;
    if (discriminant < 0) {
        return std::make_pair(std::numeric_limits<real>::infinity(), Q_ZERO);
    }
    discriminant = sqrt(discriminant);
    real distance = -od - discriminant;
    if (distance < 0) {
        distance = discriminant - od;
    }
    if (distance > 0) {
        // On a unit ball the point of intersection is the same as the vector of outgoing normal
        const quaternion normal = origin + distance * direction;
        return std::make_pair(distance, normal);
    }
    return std::make_pair(std::numeric_limits<real>::infinity(), Q_ZERO);
}