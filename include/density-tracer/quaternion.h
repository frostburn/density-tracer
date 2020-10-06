#ifndef DENSITY_TRACER_QUATERNION_H_GUARD
#define DENSITY_TRACER_QUATERNION_H_GUARD

#include <iostream>

typedef double real;

struct quaternion
{
    real w;
    real x;
    real y;
    real z;
};

typedef quaternion color;

const color C_BLACK = {0, 0, 0, 0};
const color C_RED = {0, 1, 0, 0};
const color C_YELLOW = {0, 1, 1, 0};
const color C_GREEN = {0, 0, 1, 0};
const color C_CYAN = {0, 0, 1, 1};
const color C_BLUE = {0, 0, 0, 1};
const color C_MAGENTA = {0, 1, 0, 1};
const color C_WHITE = {0, 1, 1, 1};

const quaternion Q_ONE = {1, 0, 0, 0};
const quaternion Q_ZERO = {0, 0, 0, 0};

const quaternion Q_W = {1, 0, 0, 0};
const quaternion Q_I = {0, 1, 0, 0};
const quaternion Q_J = {0, 0, 1, 0};
const quaternion Q_K = {0, 0, 0, 1};

std::ostream& operator<<(std::ostream& os, const quaternion& q);

real norm(const quaternion& q);

real norm2(const quaternion& q);

quaternion conjugate(const quaternion& q);

quaternion inverse(const quaternion& q);

quaternion operator -(const quaternion& q);

quaternion operator *(const quaternion& a, const real& b);

quaternion operator *(const real& b, const quaternion& a);

quaternion operator /(const quaternion& a, const real& b);

quaternion operator /(const real& a, const quaternion& b);

quaternion operator +(const quaternion& a, const quaternion& b);

quaternion operator +(const real& a, const quaternion& b);

quaternion operator -(const quaternion& a, const quaternion& b);

quaternion operator *(const quaternion& a, const quaternion& b);

quaternion operator /(const quaternion& a, const quaternion& b);

bool operator ==(const quaternion& a, const quaternion& b);

bool operator ==(const quaternion& a, const real& b);

bool operator ==(const real& a, const quaternion& b);

bool isclose(const quaternion& a, const quaternion& b);
bool isclose(const quaternion& a, const quaternion& b, const real& tolerance);

quaternion infinity ();

quaternion imag(const quaternion& q);

quaternion normalize(const quaternion& q);

quaternion pow(const quaternion& a, const int& b);
quaternion pow(const quaternion& q, const real& b);

real dot(const quaternion& a, const quaternion& b);

quaternion project(const quaternion& a, const quaternion& b);

quaternion cross(const quaternion& a, const quaternion& b);

quaternion rotor(quaternion axis, const real& angle);

quaternion rotate(const quaternion& q, const quaternion& axis, const real& angle);

quaternion exp(const quaternion& q);

color clip_color(const color& c);

inline quaternion square(const quaternion& q) {
    return (quaternion){
        q.w*q.w - q.x*q.x - q.y*q.y - q.z*q.z,
        2*q.w*q.x,
        2*q.w*q.y,
        2*q.w*q.z
    };
}

inline quaternion cube(const quaternion& q) {
    const real w2 = q.w*q.w;
    const real xyz2 = q.x*q.x + q.y*q.y + q.z*q.z;
    const real wxyz = 3*w2 - xyz2;
    return (quaternion){
        q.w*(w2 - 3*xyz2),
        q.x*wxyz,
        q.y*wxyz,
        q.z*wxyz
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

#endif
