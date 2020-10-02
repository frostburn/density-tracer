#include <cmath>
#include <iostream>
#include <random>

#include "density-tracer/julia.h"
#include "density-tracer/ppm.h"
#include "density-tracer/image.h"
#include "density-tracer/shapes.h"

using namespace std;

int main(int argc, char *argv[])
{
    const int scale = 10;
    const int anti_alias = 2;
    const int depth = 1<<11;
    const real view_depth = 15.0;
    const int image_width = 108 * scale;
    const int image_height = 108 * scale;

    const quaternion origin = {0, 0, 0, -3};

    default_random_engine generator;
    uniform_real_distribution<real> distribution(-0.5, 0.5);
    const real du = view_depth / (real) depth;

    int width = image_width * anti_alias;
    int height = image_height * anti_alias;
    color *pixels;
    pixels = new color [width*height];
    #pragma omp parallel for
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            int idx = i + j*width;
            real view_x = (2*i - width) / (real) height * 1.0;
            real view_y = (2*j - height) / (real) height * 1.0;

            color pixel = Q_ZERO;
            for (int k = 0; k < depth; ++k) {
                real t = (depth - (k + distribution(generator))) / (real) depth * view_depth;

                quaternion ray_direction = (quaternion){0, view_x, view_y, 0} - origin;
                ray_direction = ray_direction / norm(ray_direction);

                const quaternion loc = origin + t*ray_direction;

                real r = 1e9;
                for (int n = 0; n < 3; ++n) {
                    for (int m = 0; m < 3; ++m) {
                        quaternion local_loc = loc - (quaternion){0, 0.5*n - 0.5, 0.6, m + m*m};
                        local_loc = rotate(local_loc, {0, n, m, 2}, n+3*m);
                        r = min(r, merkaba(local_loc) - 0.1);
                    }
                }
                color illumination = Q_ZERO;
                color absorption = Q_ZERO;
                if (r > 0) {
                    illumination = {0, exp(-13*r), exp(-16*r), exp(-12*r)};
                    illumination = illumination * exp(-loc.z*0.1);
                    absorption = {0, 0.01, 0.01, 0.01};
                } else {
                    illumination = Q_ZERO;
                    absorption = {0, 20.0, 30.0, 40.0};
                }

                r = norm(loc - (quaternion){0, 0, -0.6, 2}) - 0.4;
                if (r > 0) {
                    illumination = illumination + (color){0, exp(-10*r), exp(-9*r), exp(-11*r*r)};
                } else {
                    absorption = absorption + (color){0, 9, 8, 7};
                }

                pixel = pixel + du * illumination * 1.8;
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
