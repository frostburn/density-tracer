#include <cmath>
#include <iostream>
#include <random>

#include "density-tracer/julia.h"
#include "density-tracer/ppm.h"
#include "density-tracer/image.h"
#include "density-tracer/shapes.h"
#include "density-tracer/ray.h"
#include "density-tracer/traceable.h"

using namespace std;

enum ColorMode { RGB, RYGCBV, WIDE, FULL_LINEAR, FULL_RANDOM };

class Cloud : public Density {
public:
    std::pair<real, real> eval(const quaternion& location, const quaternion& direction, const real& frequency) const {
        return std::make_pair(0, 0);
        real x = location.x + 0.5;
        x = (x - floor(x + 0.5));
        real y = (location.y - 0.5);
        real r = x*x + y*y;
        real illumination = exp(-50*r*frequency);

        real z = location.z + 0.5;
        z = (z - floor(z + 0.5));
        r = z*z + y*y;
        illumination += exp(-50*r*frequency);

        real absorption = exp(-0.01*frequency) * 0.01;
        return std::make_pair(illumination, absorption);
    }
};

class Sky : public SkySphere {
public:
    real eval(const quaternion& direction, const real& frequency) const {
        return max(0.0, min(1.0, 2*sin(7*direction.x + frequency)*cos(6*direction.y + frequency*frequency)*cos(5*direction.z)));
    }
};

class Checkers : public Pigment {
public:
    real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& frequency) const {
        real c = floor(2*location.x + 0.1) + floor(2*location.y + 0.1) + floor(2*location.z + 0.1);
        return c - 2*floor(c*0.5);
    }
};

class Phong : public Pigment {
public:
    real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& frequency) const {
        quaternion light_location = {0, 1, -3, -3};
        quaternion light = normalize(light_location - location);
        quaternion reflection = 2 * dot(light, normal) * normal - light;
        real intensity = dot(light, normal) + pow(fabs(dot(reflection, direction)), 2);
        return intensity * 0.5;
    }
};

class Mirror : public Pigment {
public:
    real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& frequency) const {
        return 0.05;
    }
};

class GlassTransparency : public Pigment {
public:
    real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& frequency) const {
        return 0.9;
    }
};

class GlassIOR : public Pigment {
public:
    real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& frequency) const {
        const real mu = (frequency - FREQ_RED) / (FREQ_VIOLET - FREQ_RED);
        return 1.03 + mu * 0.02;
    }
};

int main(int argc, char *argv[])
{
    const int scale = 1;
    const int anti_alias = 1;
    const int num_samples = 1<<10;
    const int color_samples = 1<<8;
    const ColorMode color_mode = WIDE;
    const real clip_start = 0.05;
    const real clip_end = 100.0;
    const int image_width = 108 * scale;
    const int image_height = 108 * scale;
    const int max_reflections = 16;

    const quaternion origin = {0, 0, 0, -3};
    const real air_ior = 1.0;

    const Sky sky_sphere;
    const Cloud density;
    const Phong phong;
    const Mirror mirror;
    const GlassTransparency transparency;
    const GlassIOR ior;
    const Checkers checkers;
    const Black black;
    const White white;

    vector<shared_ptr<Traceable>> objects;

    shared_ptr<Tetrahedron> tetra_a = make_shared<Tetrahedron>();
    shared_ptr<Tetrahedron> tetra_b = make_shared<Tetrahedron>();
    shared_ptr<Rotate> tetra_c = make_shared<Rotate>(tetra_b, Q_I, 0.5*M_PI);
    shared_ptr<Merge> tetra_d = make_shared<Merge>(tetra_a, tetra_c);
    shared_ptr<Rotate> tetra_e = make_shared<Rotate>(tetra_d, Q_J, 0.25*M_PI);
    shared_ptr<Rotate> tetra_f = make_shared<Rotate>(tetra_e, Q_I, 0.25*M_PI);
    shared_ptr<Translate> tetra = make_shared<Translate>(tetra_f, (quaternion){0, 0, -0.75, 5});
    tetra->pigment = &phong;
    tetra->reflectivity = &mirror;
    tetra->transparency = &transparency;
    tetra->ior = &ior;
    objects.push_back(tetra);

    shared_ptr<Plane> plane = make_shared<Plane>();
    shared_ptr<Translate> plane_t = make_shared<Translate>(plane, (quaternion){0, 0, 1, 0});
    plane_t->pigment = &checkers;
    plane_t->reflectivity = &black;
    plane_t->transparency = &black;
    plane_t->ior = &white;
    objects.push_back(plane_t);

    real dc;
    int num_color_samples;
    if (color_mode == RGB) {
        dc = 1.0;
        num_color_samples = 3;
    } else if (color_mode == RYGCBV) {
        dc = 1.0 / 3.0;
        num_color_samples = 6;
    } else if (color_mode == WIDE) {
        dc = 0.22;
        num_color_samples = 13;
    } else {
        dc = 2.0 / (real) color_samples;
        num_color_samples = color_samples;
    }

    std::default_random_engine generator;
    std::uniform_real_distribution<real> distribution(FREQ_NEAR_INFRARED, FREQ_NEAR_ULTRAVIOLET);

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
            cerr << (100 * progress) / height << "%" << endl;
        }
        for (int i = 0; i < width; ++i) {
            int idx = i + j*width;
            real view_x = (2*i - width) / (real) height * 1.0;
            real view_y = (2*j - height) / (real) height * 1.0;

            const quaternion ray_direction = normalize((quaternion){0, view_x + 0.001, view_y - 0.000321, 0} - origin);

            pixels[idx] = C_BLACK;
            for (int k = 0; k < num_color_samples; ++k) {
                real frequency;
                switch(color_mode) {
                    case RGB:
                        switch(k) {
                            case 0: frequency = FREQ_RED; break;
                            case 1: frequency = FREQ_GREEN; break;
                            case 2: frequency = FREQ_BLUE; break;
                        } break;
                    case RYGCBV:
                        switch(k) {
                            case 0: frequency = FREQ_RED; break;
                            case 1: frequency = FREQ_GREEN; break;
                            case 2: frequency = FREQ_BLUE; break;
                            case 3: frequency = FREQ_YELLOW; break;
                            case 4: frequency = FREQ_CYAN; break;
                            case 5: frequency = FREQ_VIOLET; break;
                        } break;
                    case WIDE: frequency = FREQ_RED + (k - 1) * 0.5; break;
                    case FULL_LINEAR: frequency = FREQ_NEAR_INFRARED + (FREQ_NEAR_ULTRAVIOLET - FREQ_NEAR_INFRARED) * (k + 0.5) / (real) num_color_samples; break;
                    default: frequency = distribution(generator); break;
                }
                const RayPath path(origin + clip_start*ray_direction, ray_direction, objects, clip_end - clip_start, max_reflections, frequency, air_ior);
                const real amplitude = path.eval(density, sky_sphere, num_samples);
                pixels[idx] = pixels[idx] + frequency_to_rgb(frequency) * amplitude * dc;
            }
        }
    }

    color *aa_pixels = downscale(pixels, image_width, image_height, anti_alias);
    delete[] pixels;
    cout_ppm(aa_pixels, image_width, image_height);
    delete[] aa_pixels;
    return EXIT_SUCCESS;
}
