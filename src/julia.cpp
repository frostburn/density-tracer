#include <cmath>
#include <utility>

#include "density-tracer/julia.h"

const real BAILOUT = 128;
const real EPSILON = 1e-12;

std::pair<real, quaternion> mandelbrot(quaternion q, const quaternion& c, const quaternion& dc, const int& exponent, const int& num_iter) {
    int i = 0;
    quaternion dq = Q_ONE;
    for (; i < num_iter; ++i) {
        real r = norm2(q);
        if (r > BAILOUT) {
            break;
        }
        dq = dq*exponent*pow(q, exponent-1) + dc;
        q = pow(q, exponent) + c;
    }
    real r = norm(q);
    if (r < M_E) {
        return std::make_pair(-r, dq);
    }
    r = log(log(r)) / log(exponent) + (num_iter - i);
    dq = dq * conjugate(q);
    real dr = norm(dq);
    dq = dq / (dr + (dr < EPSILON));
    return std::make_pair(r, dq);
}
