#include <cmath>
#include <utility>

#include "density-tracer/julia.h"

const real BAILOUT = 256;
const real EPSILON = 1e-12;

std::pair<real, quaternion> mandelbrot(quaternion q, const quaternion& c, const quaternion& dc, const int& exponent, const int& num_iter) {
    int i = 0;
    quaternion dq = Q_ONE;
    real r = 0;
    for (; i < num_iter; ++i) {
        r = norm2(q);
        if (r > BAILOUT) {
            break;
        }
        dq = dq*exponent*pow(q, exponent-1) + dc;
        q = pow(q, exponent) + c;
    }
    if (r < M_E) {
        return std::make_pair(-sqrt(r), dq);
    }
    r = log(log(r)*0.5) / log(exponent) + (num_iter - 1 - i);
    dq = dq * conjugate(q);
    real dr = norm(dq);
    dq = dq / (dr + (dr < EPSILON));
    return std::make_pair(r, dq);
}


real classic_mandelbrot(quaternion q, const quaternion& c, const int& num_iter) {
    int i = 0;
    real w = q.w;
    real x = q.x;
    real y = q.y;
    real z = q.z;
    real r = 0;
    for (; i < num_iter; ++i) {
        r = w*w + x*x + y*y + z*z;
        if (r > BAILOUT) {
            break;
        }
        real t = w;
        w = w*w - x*x - y*y - z*z + c.w;
        x = 2*t*x + c.x;
        y = 2*t*y + c.y;
        z = 2*t*z + c.z;
    }
    if (r < M_E*M_E) {
        return -sqrt(r);
    }
    return log(log(r)) * 1.4426950408889634 - 1.0 + (num_iter - 1 - i);
}


real abc_julia(quaternion q, const quaternion& a, const quaternion& b, const quaternion& c, const int& num_iter) {
    const quaternion d = a + b;
    const quaternion e = a - b;
    int i = 0;
    real r = 0;
    for (; i < num_iter; ++i) {
        r = norm2(q);
        if (r > BAILOUT) {
            break;
        }
        // Optimized version of: q = q*q + q*a + b*q + c
        const quaternion f = q + d;
        q = (quaternion) {
            q.w*f.w - q.x*f.x - q.y*f.y - q.z*f.z + c.w,
            q.w*f.x + q.x*f.w + q.y*e.z - q.z*e.y + c.x,
            q.w*f.y + q.y*f.w - q.x*e.z + q.z*e.x + c.y,
            q.w*f.z + q.z*f.w + q.x*e.y - q.y*e.x + c.z,
        };
    }
    if (r < M_E*M_E) {
        return -sqrt(r);
    }
    return log(log(r)) * 1.4426950408889634 - 2.0 + num_iter - i;
}
