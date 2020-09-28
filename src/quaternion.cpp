#include <iostream>
#include <limits>
#include <cmath>

#include "density-tracer/quaternion.h"

const real EPSILON = 1e-12;

std::ostream& operator<<(std::ostream& os, const quaternion& q)
{
    return os << "(quaternion){" <<  q.w << ", " << q.x << ", " << q.y << ", " << q.z << "}"; 
}

real norm(const quaternion& q) {
    return sqrt(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
}

real norm2(const quaternion& q) {
    return q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z;
}

quaternion conjugate(const quaternion& q) {
    return (quaternion){q.w, -q.x, -q.y, -q.z};
}

quaternion inverse(const quaternion& q) {
    real ir = -1.0 / (q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
    return (quaternion){-q.w*ir, q.x*ir, q.y*ir, q.z*ir};
}

quaternion operator -(const quaternion& q) {
    return (quaternion){-q.w, -q.x, -q.y, -q.z};
}

quaternion operator *(const quaternion& a, const real& b) {
    return (quaternion){a.w * b, a.x * b, a.y * b, a.z * b};
}

quaternion operator *(const real& b, const quaternion& a) {
    return (quaternion){a.w * b, a.x * b, a.y * b, a.z * b};
}

quaternion operator /(const quaternion& a, const real& b) {
    return a * (1.0 / b);
}

quaternion operator /(const real& a, const quaternion& b) {
    real r = b.w*b.w + b.x*b.x + b.y*b.y + b.z*b.z;
    return (a / r) * conjugate(b);
}

quaternion operator +(const quaternion& a, const quaternion& b) {
    return (quaternion){a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z};
}

quaternion operator +(const real& a, const quaternion& b) {
    return (quaternion){a + b.w, b.x, b.y, b.z};
}

quaternion operator -(const quaternion& a, const quaternion& b) {
    return (quaternion){a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z};
}

quaternion operator *(const quaternion& a, const quaternion& b) {
    return (quaternion){
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w
    };
}

quaternion operator /(const quaternion& a, const quaternion& b) {
    real r = b.w*b.w + b.x*b.x + b.y*b.y + b.z*b.z;
    return a * conjugate(b) / r;
}

bool operator ==(const quaternion& a, const quaternion& b) {
    return a.w == b.w && a.x == b.x && a.y == b.y && a.z == b.z;
}

bool operator ==(const quaternion& a, const real& b) {
    return a.w == b && a.x == 0 && a.y == 0 && a.z == 0;
}

bool operator ==(const real& a, const quaternion& b) {
    return a == b.w && b.x == 0 && b.y == 0 && b.z == 0;
}

bool isclose(const quaternion& a, const quaternion& b) {
    return fabs(a.w - b.w) < EPSILON && fabs(a.x - b.x) < EPSILON && fabs(a.y - b.y) < EPSILON && fabs(a.z - b.z) < EPSILON;
}

quaternion infinity () {
    real inf = std::numeric_limits<real>::infinity();
    return (quaternion){inf, inf, inf, inf};
}

quaternion imag(const quaternion& q) {
    return (quaternion){0, q.x, q.y, q.z};
}

inline quaternion square(const quaternion& q) {
    return (quaternion){
        q.w*q.w - q.x*q.x - q.y*q.y - q.z*q.z,
        2*q.w*q.x,
        2*q.w*q.y,
        2*q.w*q.z
    };
}

// Quaternions in the same imaginary plane have zero cross terms.
inline quaternion aligned_mul(const quaternion& a, const quaternion& b) {
    return (quaternion){
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
        a.w*b.x + a.x*b.w,
        a.w*b.y + a.y*b.w,
        a.w*b.z + a.z*b.w
    };
}

quaternion pow(const quaternion& a, int b) {
    quaternion res = Q_ONE;
    quaternion m = a;
    int exponent = abs(b);
    while (exponent > 0) {
        if (exponent & 1) {
            res = aligned_mul(res, m);
        }
        m = square(m);
        exponent >>= 1;
    }
    if (b < 0) {
        return inverse(res);
    }
    return res;
}

real dot(const quaternion& a, const quaternion& b) {
    return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
}

quaternion project(const quaternion& a, const quaternion& b) {
    return dot(a, b) * a;
}

quaternion cross_align(const quaternion& a, const quaternion& b) {
    // Assumes |a| = 1
    quaternion c = b;
    real angle = dot(a, b);
    if (angle > EPSILON) {
        c = b / angle - a;
    } else if (angle < -EPSILON) {
        c = a - b / angle;
    }
    return c / norm(c);
}

quaternion exp(const quaternion& q) {
    quaternion im = imag(q);
    real r = norm(im);
    return exp(q.w) * (cos(r) + sin(r) * im);
}

color clip_color(const color& a) {
    color pixel = a;
    if (pixel.w > 1) {
        pixel.w = 1;
    }
    if (pixel.w < 0) {
        pixel.w = 0;
    }
    if (std::isnan(pixel.w)) {
        pixel.w = 0.5;
    }
    if (pixel.x > 1) {
        pixel.x = 1;
    }
    if (pixel.x < 0) {
        pixel.x = 0;
    }
    if (std::isnan(pixel.x)) {
        pixel.x = 0.5;
    }
    if (pixel.y > 1) {
        pixel.y = 1;
    }
    if (pixel.y < 0) {
        pixel.y = 0;
    }
    if (std::isnan(pixel.y)) {
        pixel.y = 0.5;
    }
    if (pixel.z > 1) {
        pixel.z = 1;
    }
    if (pixel.z < 0) {
        pixel.z = 0;
    }
    if (std::isnan(pixel.z)) {
        pixel.z = 0.5;
    }
    return pixel;
}