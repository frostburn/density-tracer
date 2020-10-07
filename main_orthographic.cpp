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
    const int anti_alias = 3;
    const int depth = 1 << 12;
    const real clip_start = -2.0;
    const real clip_end = 2.0;
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

    // MultiBranchMandelbrot brot(5, 2, 11, 0, false, 1<<12, potential, min_r, numeric_limits<real>::infinity());

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
                quaternion q = {1.8*view_x+0.2, 1.8*view_y+0.05, t, 0};

                q = rotate(q, {0, 1, 2, -3}, -2.5);
                quaternion rot = rotor({0, 4, -5, 6}, 1.15);
                q = rot*q*rot;

                // real r = orthoplex(q, {0.4, 0.5, 0.6, -0.7}, 4, 22);
                real r = pentatope(q, q, 12);

                // auto [inside, r] = multibranch_pentatope(q, q, 6, 1, true, 1<<6, shifted_potential, min_r, numeric_limits<real>::infinity());


                color illumination, absorption;
                if (r < 0) {
                    illumination = C_BLACK;
                    absorption = {0, 15, 15, 15};
                } else {
                    illumination = {0, 10*exp(-0.05*r*r), 12*exp(-0.3*r*r), 22*exp(-0.7*r*r)};
                    absorption = {0, 0.1, 0.1, 0.1};
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
