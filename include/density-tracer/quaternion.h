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

#define Q_ONE ((quaternion){1, 0, 0, 0})
#define Q_ZERO ((quaternion){0, 0, 0, 0})

#define Q_W ((quaternion){1, 0, 0, 0})
#define Q_I ((quaternion){0, 1, 0, 0})
#define Q_J ((quaternion){0, 0, 1, 0})
#define Q_K ((quaternion){0, 0, 0, 1})

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

real dot(const quaternion& a, const quaternion& b);

quaternion project(const quaternion& a, const quaternion& b);

quaternion cross(const quaternion& a, const quaternion& b);

quaternion rotor(quaternion axis, const real& angle);

quaternion rotate(const quaternion& q, const quaternion& axis, const real& angle);

quaternion exp(const quaternion& q);

color clip_color(const color& c);

#endif
