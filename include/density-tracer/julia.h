#ifndef DENSITY_TRACER_JULIA_H_GUARD
#define DENSITY_TRACER_JULIA_H_GUARD

#include <cmath>
#include <utility>
#include <functional>
#include <vector>

#include "density-tracer/quaternion.h"

std::pair<real, quaternion> mandelbrot(quaternion q, const quaternion& c, const quaternion& dc, const int& exponent, const int& num_iter);

real classic_mandelbrot(quaternion q, const quaternion& c, const int& num_iter);

real abc_julia(quaternion q, const quaternion& a, const quaternion& b, const quaternion& c, const int& num_iter);

real orthoplex(quaternion q, const quaternion& c, const int& exponent, const int& num_iter);

real pentatope(quaternion q, const quaternion& c, const int& num_iter);

real min_axis_mandelbrot(quaternion q, const quaternion c, const std::vector<quaternion>& ax, const int& exponent, const int& num_iter);

real multi_c_julia(quaternion q, const std::vector<quaternion>& cs, const int& exponent, const int& num_iter);

typedef std::function<real(const real&, const real&)> Reducer;
typedef std::function<real(const real&, const int&, const real&)> R2Mapper;

real min_r(const real& a, const real& b);

real shifted_potential(const real& r2, const int& remaining_iterations, const real& inverse_log_exponent);

class MultiBranchMandelbrot {
    // Parameters
    int denominator;
    int num_iter;
    real bailout;
    unsigned long long inside_cutoff;
    bool clip_outside;
    Reducer reducer;
    R2Mapper r2_mapper;
    real empty;

    // Temporary variables
    real exponent;
    quaternion c;
    real inverse_log_exponent;
    real *roots_of_unity;

public:
    MultiBranchMandelbrot(const int& numerator, const int& denominator, const int& num_iter, const unsigned long long& inside_cutoff, const bool& clip_outside, const real& bailout, R2Mapper r2_mapper, Reducer reducer, const real& empty);
    ~MultiBranchMandelbrot();
    void accumulate(unsigned long long *inside_counter, real *outside_accumulator, quaternion q, const quaternion& c, const int& num_iter) const;
    std::pair<unsigned long long int, real> eval(const quaternion& q, const quaternion& c) const;
};

std::pair<unsigned long long, real> multibranch_pentatope(quaternion q, const quaternion& c, const int& num_iter, const unsigned long long& inside_cutoff, const bool& clip_outside, const real& bailout, R2Mapper r2_mapper, Reducer reducer, const real& empty);

#endif
