#include <cmath>
#include <iostream>
#include <random>

#include "density-tracer/julia.h"
#include "density-tracer/ppm.h"
#include "density-tracer/image.h"

using namespace std;

int main(int argc, char *argv[])
{
    int scale = 10;
    int anti_alias = 3;
    int depth = 1<<11;
    int image_width = 108 * scale;
    int image_height = 108 * scale;

    const quaternion c = {0.5, 0.2, -0.3, 0.35};

    default_random_engine generator;
    uniform_real_distribution<real> distribution(-0.5, 0.5);
    real du = 2.0 / (real) depth;

    int width = image_width * anti_alias;
    int height = image_height * anti_alias;
    color *pixels;
    pixels = new color [width*height];
    #pragma omp parallel for
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            int idx = i + j*width;
            real x = (2*i - width) / (real) height * 0.9;
            real y = (2*j - height) / (real) height * 1.3;

            color pixel = Q_ZERO;
            for (int k = 0; k < depth; ++k) {
                real z = (2*(k + distribution(generator)) - depth) / (real) depth;

                quaternion q = {x, y, z, 0};
                real r = classic_mandelbrot(q, c, 13);
                color illumination = Q_ZERO;
                color absorption = Q_ZERO;
                if (r > 0) {
                    illumination = {0, exp(-0.5*r), exp(-0.4*r), exp(-0.3*r)};
                    absorption = {0, 0.1, 0.1, 0.1};
                } else {
                    illumination = Q_ZERO;
                    absorption = {0, 1.0, 1.0, 1.0};
                }
                pixel = pixel + du * illumination * 2;
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
