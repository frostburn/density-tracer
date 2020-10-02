#include <cmath>
#include <iostream>
#include <random>

#include "density-tracer/julia.h"
#include "density-tracer/ppm.h"
#include "density-tracer/image.h"
#include "density-tracer/shapes.h"
#include "density-tracer/ray.h"

using namespace std;

class Cloud : public Density {
public:
    std::pair<color, color> eval(const quaternion& location, const quaternion& direction) const {
        real x = location.x + 0.5;
        x = (x - floor(x + 0.5));
        real y = (location.y - 0.5);
        real r = x*x + y*y;
        color illumination = {0, exp(-150*r), 1.2*exp(-162*r), 1.2*exp(-174*r)};
        illumination = illumination * 1.5;

        real z = location.z + 0.5;
        z = (z - floor(z + 0.5));
        r = z*z + y*y;
        illumination = illumination + ((quaternion){0, exp(-250*r), 1.2*exp(-262*r), 1.2*exp(-274*r)}) * 2.5;

        color absorption = {0, 0.2*exp(location.y*location.y), 0.2*exp(location.y), 0.2*exp(location.y)};
        return std::make_pair(illumination, absorption);
    }
};

class Sky : public SkySphere {
public:
    color eval(const quaternion& direction) const {
        return {0, fabs(direction.x), fabs(direction.y), fabs(direction.z)};
    }
};

class Phong : public Pigment {
public:
    color eval(const quaternion& location, const quaternion& direction, const quaternion& normal) const {
        quaternion light_location = {0, 1, -3, -3};
        quaternion light = normalize(light_location - location);
        quaternion reflection = 2 * dot(light, normal) * normal - light;
        real intensity = dot(light, normal) + pow(dot(-reflection, -direction), 2);
        return {0, intensity, intensity, intensity};
    }
};

int main(int argc, char *argv[])
{
    const int scale = 10;
    const int anti_alias = 2;
    const int num_samples = 1<<10;
    const real clip_start = 0.05;
    const real clip_end = 20.0;
    const int image_width = 108 * scale;
    const int image_height = 108 * scale;
    const int max_reflections = 128;

    const quaternion origin = {0, 0, 0, -3};

    const Sky sky_sphere;
    const Cloud density;
    const Phong phong;

    vector<shared_ptr<Traceable>> objects;

    shared_ptr<Ball> ball = make_shared<Ball>();
    ball->location = {0, 1.25, -0.45, 4};
    ball->pigment = &phong;
    ball->reflectivity = 0.2;
    objects.push_back(ball);

    shared_ptr<Tetrahedron> tetra = make_shared<Tetrahedron>();
    tetra->location = {0, -1.25, 0, 4};
    tetra->right_transform = rotor({0, 1, 2, 3}, 1.5);
    tetra->left_transform = inverse(tetra->right_transform);
    tetra->pigment = &phong;
    tetra->reflectivity = 0.1;
    objects.push_back(tetra);

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

            const quaternion ray_direction = normalize((quaternion){0, view_x, view_y, 0} - origin);

            const RayPath path(origin + clip_start*ray_direction, ray_direction, objects, sky_sphere, clip_end - clip_start, max_reflections);

            pixels[idx] = path.eval(density, num_samples);
        }
    }

    color *aa_pixels = downscale(pixels, image_width, image_height, anti_alias);
    delete[] pixels;
    cout_ppm(aa_pixels, image_width, image_height);
    delete[] aa_pixels;
    return EXIT_SUCCESS;
}
