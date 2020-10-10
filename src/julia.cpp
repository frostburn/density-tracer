#include <algorithm>
#include <cmath>
#include <limits>

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
    real r2 = 0;
    for (; i < num_iter; ++i) {
        r2 = norm2(q);
        if (r2 > BAILOUT) {
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
    if (r2 < M_E*M_E) {
        return -sqrt(r2);
    }
    return log(log(r2)) * 1.4426950408889634 - 2.0 + num_iter - i;
}

real orthoplex(quaternion q, const quaternion& c, const int& exponent, const int& num_iter) {
    int i = 0;
    real r2 = 0;
    for (; i < num_iter; ++i) {
        r2 = norm2(q);
        if (r2 > BAILOUT) {
            break;
        }
        q = 0.25*(pow(q, exponent) - Q_I * pow(Q_I*q, exponent) - Q_J * pow(Q_J*q, exponent) - Q_K * pow(Q_K*q, exponent)) + c;
    }
    if (r2 < M_E*M_E) {
        return -sqrt(r2);
    }
    return log(log(r2)*0.5) / log(exponent) + (num_iter - 1 - i);
}

const quaternion P1 = {-0.25, 0.5590169943749475, 0.5590169943749475, 0.5590169943749475};
const quaternion P2 = {-0.25, 0.5590169943749475, -0.5590169943749475, -0.5590169943749475};
const quaternion P3 = {-0.25, -0.5590169943749475, 0.5590169943749475, -0.5590169943749475};
const quaternion P4 = {-0.25, -0.5590169943749475, -0.5590169943749475, 0.5590169943749475};

const quaternion IP1 = {-0.25, -0.5590169943749475, -0.5590169943749475, -0.5590169943749475};
const quaternion IP2 = {-0.25, -0.5590169943749475, 0.5590169943749475, 0.5590169943749475};
const quaternion IP3 = {-0.25, 0.5590169943749475, -0.5590169943749475, 0.5590169943749475};
const quaternion IP4 = {-0.25, 0.5590169943749475, 0.5590169943749475, -0.5590169943749475};

void multibranch_pentatope_accumulate(unsigned long long *inside_counter, real *outside_accumulator, quaternion q, const quaternion& c, const int& num_iter, const unsigned long long int& inside_cutoff, const bool& clip_outside, const real& bailout, R2Mapper r2_mapper, Reducer reducer) {
    if (inside_cutoff && *inside_counter >= inside_cutoff) {
        return;
    }
    const real w2 = q.w*q.w;
    const real im2 = q.x*q.x + q.y*q.y + q.z*q.z;
    const real r2 = w2 + im2;
    if (r2 >= bailout || num_iter == 0) {
        if (clip_outside || *inside_counter == 0) {
            const real v = r2_mapper(r2, num_iter, 1.4426950408889634);
            *outside_accumulator = reducer(*outside_accumulator, v);
        }
        return;
    }
    if (num_iter == 0) {
        *inside_counter = *inside_counter + 1;
        return;
    }
    multibranch_pentatope_accumulate(inside_counter, outside_accumulator, {w2 - im2 + c.w, 2*q.w*q.x + c.x, 2*q.w*q.y + c.y, 2*q.w*q.z + c.z}, c, num_iter - 1, inside_cutoff, clip_outside, bailout, r2_mapper, reducer);
    multibranch_pentatope_accumulate(inside_counter, outside_accumulator, q*IP1*q + c, c, num_iter - 1, inside_cutoff, clip_outside, bailout, r2_mapper, reducer);
    multibranch_pentatope_accumulate(inside_counter, outside_accumulator, q*IP2*q + c, c, num_iter - 1, inside_cutoff, clip_outside, bailout, r2_mapper, reducer);
    multibranch_pentatope_accumulate(inside_counter, outside_accumulator, q*IP3*q + c, c, num_iter - 1, inside_cutoff, clip_outside, bailout, r2_mapper, reducer);
    multibranch_pentatope_accumulate(inside_counter, outside_accumulator, q*IP4*q + c, c, num_iter - 1, inside_cutoff, clip_outside, bailout, r2_mapper, reducer);
}

std::pair<unsigned long long, real> multibranch_pentatope(quaternion q, const quaternion& c, const int& num_iter, const unsigned long long& inside_cutoff, const bool& clip_outside, const real& bailout, R2Mapper r2_mapper, Reducer reducer, const real& empty) {
    unsigned long long int inside_counter = 0;
    real outside_accumulator = empty;

    multibranch_pentatope_accumulate(&inside_counter, &outside_accumulator, q, c, num_iter, inside_cutoff, clip_outside, bailout, r2_mapper, reducer);

    return std::make_pair(inside_counter, outside_accumulator);
}

real pentatope(quaternion q, const quaternion& c, const int& num_iter) {
    int i = 0;
    real r2 = norm2(q);
    for (; i < num_iter; ++i) {
        if (r2 > BAILOUT) {
            break;
        }
        const quaternion q0 = square(q) + c;
        const quaternion q1 = q*IP1*q + c;
        const quaternion q2 = q*IP2*q + c;
        const quaternion q3 = q*IP3*q + c;
        const quaternion q4 = q*IP4*q + c;

        q = q0;
        r2 = norm2(q);

        real r2_n = norm2(q1);
        if (r2_n < r2) {
            q = q1;
            r2 = r2_n;
        }

        r2_n = norm2(q2);
        if (r2_n < r2) {
            q = q2;
            r2 = r2_n;
        }

        r2_n = norm2(q3);
        if (r2_n < r2) {
            q = q3;
            r2 = r2_n;
        }

        r2_n = norm2(q4);
        if (r2_n < r2) {
            q = q4;
            r2 = r2_n;
        }
    }
    if (r2 < M_E*M_E) {
        return -sqrt(r2);
    }
    return log(log(r2)) * 1.4426950408889634 - 2.0 + num_iter - i;
}

real multi_c_julia(quaternion q, const std::vector<quaternion>& cs, const int& exponent, const int& num_iter) {
    int i = 0;
    real r2 = 0;
    for (; i < num_iter; ++i) {
        r2 = norm2(q);
        if (r2 > BAILOUT) {
            break;
        }
        const quaternion qn = pow(q, exponent);
        real min_r2 = std::numeric_limits<real>::infinity();
        for (std::vector<quaternion>::const_iterator c = cs.begin(); c != cs.end(); ++c) {
            const quaternion qnc = qn + (*c);
            const real r2_n = norm2(qnc);
            if (r2_n < min_r2) {
                min_r2 = r2_n;
                q = qnc;
            }
        }
    }
    if (r2 < M_E*M_E) {
        return -sqrt(r2);
    }
    return log(log(r2)*0.5) / log(exponent) + (num_iter - 1 - i);
}

real min_r(const real& a, const real& b) {
    return std::min(a, b);
}

real shifted_potential(const real& r2, const int& remaining_iterations, const real& inverse_log_exponent) {
    return log(log(r2 + 1.0 + EPSILON)*0.5) * inverse_log_exponent + remaining_iterations;
}


MultiBranchMandelbrot::MultiBranchMandelbrot(const int& numerator, const int& denominator, const int& num_iter, const unsigned long long& inside_cutoff, const bool& clip_outside, const real& bailout, R2Mapper r2_mapper, Reducer reducer, const real& empty) {
    this->roots_of_unity = new real[2*denominator];
    for (int i = 0; i < denominator; ++i) {
        this->roots_of_unity[2*i+0] = cos(i*2*M_PI/(real)denominator);
        this->roots_of_unity[2*i+1] = sin(i*2*M_PI/(real)denominator);
    }
    this->denominator = denominator;
    this->exponent = numerator / (real)denominator;
    this->inverse_log_exponent = 1.0 / log(this->exponent);

    this->num_iter = num_iter;
    this->inside_cutoff = inside_cutoff;
    this->clip_outside = clip_outside;
    this->bailout = bailout;
    this->r2_mapper = r2_mapper;
    this->reducer = reducer;
    this->empty = empty;
}

MultiBranchMandelbrot::~MultiBranchMandelbrot() {
    delete[] this->roots_of_unity;
}

void MultiBranchMandelbrot::accumulate(unsigned long long *inside_counter, real *outside_accumulator, quaternion q, const quaternion& c, const int& num_iter) const {
    if (this->inside_cutoff && *inside_counter >= this->inside_cutoff) {
        return;
    }
    const real w2 = q.w*q.w;
    const real im2 = q.x*q.x + q.y*q.y + q.z*q.z;
    const real r2 = w2 + im2;
    if (r2 >= this->bailout || num_iter == 0) {
        if (!this->clip_outside || *inside_counter == 0) {
            const real v = this->r2_mapper(r2, num_iter, this->inverse_log_exponent);
            *outside_accumulator = this->reducer(*outside_accumulator, v);
        }
        return;
    }
    if (num_iter == 0) {
        *inside_counter = *inside_counter + 1;
        return;
    }
    const real rim = sqrt(im2);
    if (fabs(rim) < EPSILON) {
        const real theta = atan2(q.x, q.w) * this->exponent;
        const real r = pow(r2, this->exponent*0.5);
        const real w = cos(theta) * r;
        const real im = sin(theta) * r;
        for (int i = 0; i < this->denominator; ++i) {
            real w_ = w * this->roots_of_unity[2*i] + im * this->roots_of_unity[2*i + 1];
            real im_ = im * this->roots_of_unity[2*i] - w *this->roots_of_unity[2*i + 1];
            this->accumulate(
                inside_counter, outside_accumulator,
                (quaternion){w_ + c.w, im_ + c.x, c.y, c.z},
                c, num_iter - 1
            );
        }
    } else {
        const real theta = atan2(rim, q.w) * this->exponent;
        const real rimn = 1.0 / rim;
        const real r = pow(r2, this->exponent*0.5);
        const real w = cos(theta) * r;
        const real im = sin(theta) * r;
        const real wn = w * rimn;
        const real imn = im * rimn;
        for (int i = 0; i < this->denominator; ++i) {
            real w_ = w * this->roots_of_unity[2*i] + im * this->roots_of_unity[2*i + 1];
            real imn_ = imn * this->roots_of_unity[2*i] - wn * this->roots_of_unity[2*i + 1];
            this->accumulate(
                inside_counter, outside_accumulator,
                (quaternion){w_ + c.w, q.x*imn_ + c.x, q.y*imn_ + c.y, q.z*imn_ + c.z},
                c, num_iter - 1
            );
        }
    }
}

std::pair<unsigned long long int, real> MultiBranchMandelbrot::eval(const quaternion& q, const quaternion& c) const {
    unsigned long long int inside_counter = 0;
    real outside_accumulator = this->empty;

    this->accumulate(&inside_counter, &outside_accumulator, q, c, this->num_iter);

    return std::make_pair(inside_counter, outside_accumulator);
}
