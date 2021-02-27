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

int main(int argc, char *argv[])
{
    const int scale = 10;
    const int anti_alias = 2;
    const int depth = 1 << 11;
    const real clip_start = -2.5;
    const real clip_end = 2.5;
    const int image_width = 108 * scale;
    const int image_height = 108 * scale;

    int width = image_width * anti_alias;
    int height = image_height * anti_alias;
    color *pixels;
    pixels = new color [width*height];
    int progress = 0;

    default_random_engine generator;
    uniform_real_distribution<real> distribution(0, 1);
    const real du = (clip_end - clip_start) / (real) depth;

    // NonEscapingMultiBranch brot(-5, 2, 5, max_q);

    const real PHI = 1.618033988749895;
    vector<quaternion> ax = {
        {0, 1, PHI, 0}, {0, 1, -PHI, 0},
        {1, PHI, 0, 0}, {1, -PHI, 0, 0},
        {PHI, 0, 1, 0}, {-PHI, 0, 1, 0}
    };
    // vector<quaternion> ax = {
    //     Q_ONE,
    //     {-pow(5, -0.5), 1, 1, 1},
    //     {-pow(5, -0.5), 1, -1, -1},
    //     {-pow(5, -0.5), -1, 1, -1},
    //     {-pow(5, -0.5), -1, -1, 1}
    // };

    // vector<quaternion> ax = {
    //     {1, 1, 1, 0},
    //     {1, -1, -1, 0},
    //     {-1, 1, -1, 0},
    //     {-1, -1, 1, 0 }
    // };


    // for (auto i = ax.begin(); i != ax.end(); ++i) {
    //     *i = normalize(*i);
    // }

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

            color pixel = C_WHITE*0.015;
            for (int k = 0; k < depth; ++k) {
                real t = (depth - k - distribution(generator)) * du + clip_start;
                quaternion q = {2*view_x, 2*view_y+0.001, t, 0};

                // quaternion rot = rotor({0, 4, -5, 6}, 1.15);
                quaternion rot = normalize({1, 0.4, 2.25, 0});
                q = rot*q*rot;

                // const quaternion c = {0.4, 0.2, -0.3, -0.213};
                // const quaternion res = brot.eval(q, c);

                const quaternion res = max_axis_nonescaping(q, q, ax, -3, 42);

                const real r = norm2(res)-2.9;

                color illumination, absorption;
                if (r < 0) {
                    illumination = {0, 0, 5*exp(5*r), 0};
                    absorption = {0, 10-r, 10-r, 10-r};
                } else {
                    illumination = {0, exp(-r*0.5), exp(-r), exp(-r*2) + exp(-(r-1)*(r-1))};
                    absorption = {0, 0.1 + 5*exp(-(r-2)*(r-2)), 0.1, 0.1};
                }

                pixel = pixel + du * illumination;
                pixel = {0, pixel.x * exp(-absorption.x*du), pixel.y * exp(-absorption.y*du), pixel.z * exp(-absorption.z*du)};
            }
            pixels[idx] = pixel;
        }
    }

    color *aa_pixels = downscale(pixels, image_width, image_height, anti_alias);
    delete[] pixels;
    cout_ppm(aa_pixels, image_width, image_height);
    delete[] aa_pixels;
    return EXIT_SUCCESS;
}
