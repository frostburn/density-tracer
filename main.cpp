#include <cmath>
#include <iostream>
#include <random>

#include "density-tracer/julia.h"
#include "density-tracer/ppm.h"
#include "density-tracer/image.h"
#include "density-tracer/shapes.h"
#include "density-tracer/ray.h"

using namespace std;

int main(int argc, char *argv[])
{
    const int scale = 1;
    const int anti_alias = 3;
    const int depth = 1<<12;
    const real clip_start = 0.05;
    const real clip_end = 30.0;
    const int image_width = 108 * scale;
    const int image_height = 108 * scale;
    const int max_reflections = 128;

    const quaternion origin = {0, 0, 0, -5};

    vector<shared_ptr<Traceable>> objects;

    shared_ptr<Ball> ball = make_shared<Ball>();
    ball->location = {0, 0, -0.45, 1};
    objects.push_back(ball);

    shared_ptr<Plane> plane = make_shared<Plane>();
    plane->location = {0, 0, 0.56, 0};
    objects.push_back(plane);

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
            view_y += 0.5;

            const quaternion ray_direction = normalize((quaternion){0, view_x, view_y, 0} - origin);

            const RayPath path(origin + clip_start*ray_direction, ray_direction, objects, max_reflections);

            color pixel = Q_ZERO;
            for (int k = 0; k < depth; ++k) {
                real t = (depth - (k + distribution(generator))) / (real) depth * (clip_end - clip_start);

                const quaternion loc = path.get_location(t);
                color illumination = Q_ZERO;
                color absorption = Q_ZERO;

                real x = loc.x + 0.5;
                x = (x - floor(x + 0.5));
                real y = (loc.y - 0.5);
                real r = x*x + y*y;
                illumination = {0, exp(-150*r), 1.2*exp(-162*r), 1.2*exp(-174*r)};
                illumination = illumination * 1.5;

                real z = loc.z + 0.5;
                z = (z - floor(z + 0.5));
                r = z*z + y*y;
                illumination = illumination + ((quaternion){0, exp(-250*r), 1.2*exp(-262*r), 1.2*exp(-274*r)}) * 2.5;

                absorption = {0, 0.2*exp(loc.y*loc.y), 0.2*exp(loc.y), 0.2*exp(loc.y)};

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
