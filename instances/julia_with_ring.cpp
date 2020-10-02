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
    const int anti_alias = 3;
    const int depth = 1<<11;
    const real clip_start = 2.0;
    const real clip_end = 5.0;
    const int image_width = 108 * scale;
    const int image_height = 108 * scale;

    const quaternion origin = {0, 0, 0, -3.1};

    default_random_engine generator;
    uniform_real_distribution<real> distribution(-0.5, 0.5);
    const real du = (clip_end - clip_start) / (real) depth;

    int width = image_width * anti_alias;
    int height = image_height * anti_alias;
    color *pixels;
    pixels = new color [width*height];
    int progress = 0;
    #pragma omp parallel for
    for (int j = 0; j < height; ++j) {
        #pragma omp atomic update
            progress++;
        if (progress % 10 == 0) {
            cerr << (100 * progress) / width << "%" << endl;
        }
        for (int i = 0; i < width; ++i) {
            int idx = i + j*width;
            real view_x = (2*i - width) / (real) height * 1.0;
            real view_y = (2*j - height) / (real) height * 1.0;

            color pixel = Q_ZERO;
            for (int k = 0; k < depth; ++k) {
                real t = (depth - (k + distribution(generator))) / (real) depth * (clip_end - clip_start) + clip_start;

                quaternion ray_direction = (quaternion){0, view_x, view_y, 0} - origin;
                ray_direction = ray_direction / norm(ray_direction);

                const quaternion loc = origin + t*ray_direction;
                color illumination = Q_ZERO;
                color absorption = Q_ZERO;


                quaternion local_loc = loc - (quaternion){0, 0, 0.2, 0.2};
                local_loc = rotate(local_loc, {0, -2, -1, 0.5}, 0.82);
                quaternion q = {local_loc.x, local_loc.y, local_loc.z, 0};
                real r = abc_julia(q, {0.42, 0.31, -0.3, -0.4}, {-0.71, 0.4, 0.39, 0.22}, {0.52, -0.3, 0.14, -0.1}, 15);

                if (r > 0) {
                    illumination = {0, 2*exp(-0.4*r), exp(-0.5*r) + 1.2*exp(-0.5*r*r), exp(-0.33*r)};
                    illumination = illumination * 2;
                    absorption = {0, 0.01, 0.015, 0.01};
                } else {
                    illumination = Q_ZERO;
                    absorption = {0, 200, 100, 300};
                }

                local_loc = loc - (quaternion){0, 0.3, -0.08, 0};
                local_loc = rotate(local_loc, {0, 1, 2, 3}, 1.5);
                r = torus(local_loc, 0.5) - 0.1;
                if (r > 0) {
                    illumination = illumination + (color){0, exp(-10*r), 1.1*exp(-18*r*r), 1.3*exp(-16*r*r)};
                } else {
                    absorption = absorption + (color){0, 100, 90, 80};
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
