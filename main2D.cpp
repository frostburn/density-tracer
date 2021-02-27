#include <cmath>
#include <iostream>
#include <random>

#include "density-tracer/julia.h"
#include "density-tracer/ppm.h"
#include "density-tracer/image.h"

using namespace std;

real potential(const real& r2, const int& remaining_iterations, const real& inverse_log_exponent) {
    const real v = shifted_potential(r2, remaining_iterations, inverse_log_exponent);
    return (v - 1)*(v-1);
}

real rcos(const real& t) {
    return 0.5  - 0.5*cos(t);
}

int main(int argc, char *argv[])
{
    const int scale = 10;
    const int anti_alias = 5;
    // const int image_width = 108 * scale;
    const int image_width = 192 * scale;
    const int image_height = 108 * scale;
    const real rotation = -M_PI*0.5;

    int width = image_width * anti_alias;
    int height = image_height * anti_alias;
    const real cos_rot = cos(rotation);
    const real sin_rot = sin(rotation);
    color *pixels;
    pixels = new color [width*height];
    int progress = 0;

    // NonEscapingMultiBranch brot(-5, 2, 7, max_q);

    // vector<quaternion> ax = {
    //      {0.3, -0.6, 0, 0}, {0.5, 0.7, 0, 0}, {-0.2, 0.61, 0, 0}
    // };

    MultiBranchMandelbrot brot(5, 2, 14, 64, true, 128, shifted_potential, min_r, 1e100);

    #pragma omp parallel for
    for (int j = 0; j < height; ++j) {
        #pragma omp atomic update
            progress++;
        if (progress % 10 == 0) {
            cerr << (100 * progress) / height << "%" << endl;
        }
        for (int i = 0; i < width; ++i) {
            int idx = i + j*width;
            real view_x = (2*i - width) / (real) height * 1.0;
            real view_y = (2*j - height) / (real) height * 1.0;

            const quaternion q = {1.0*(view_x*cos_rot + view_y*sin_rot) + 1.8, 1.0*(view_y*cos_rot - view_x*sin_rot)+0.001, 0, 0};
            //const quaternion c = {0.123, -0.2, 0.4, -0.3211};
            // const quaternion res = multi_c_nonescaping2(q, ax, 5, 160);
            auto res = brot.eval(q, q*q);


            const unsigned long long int inside = res.first;
            const real r = res.second;

            // cerr << inside << endl;

            // pixels[idx] = {0, r, 0, 0};
            pixels[idx] = {0, 0.8*exp(-r*r*0.2) - inside * 0.007, 0.8*exp(-r*0.15) + exp(-1.0*(r-10)*(r-10))*0.5 - inside*0.01 - inside*inside*0.00001, 0.9*exp(-r*0.15) - inside*0.001};

            // const real r = multi_c_julia(q, ax, 2, 6);

            // pixels[idx] = {0, res.x, res.y + 0.2*r, res.z - 0.05*r};
            // if (isnan(res.w)) {
            //     const real vr = view_x*view_x + view_y*view_y;
            //     pixels[idx] = {0, exp(-vr), exp(-2*vr), exp(-3*vr)};
            // }
            // pixels[idx] = {0, 0.1*r , 0.3*r, 0.9*r};
        }
    }

    color *aa_pixels = downscale(pixels, image_width, image_height, anti_alias);
    delete[] pixels;
    cout_ppm(aa_pixels, image_width, image_height);
    delete[] aa_pixels;
    return EXIT_SUCCESS;
}
