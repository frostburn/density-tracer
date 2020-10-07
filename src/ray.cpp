#include <cmath>
#include <limits>
#include <random>

#include "density-tracer/ray.h"

const real EPSILON = 1e-8;

color frequency_to_rgb(const real& frequency) {
    if (frequency < FREQ_NEAR_INFRARED) {
        return C_BLACK;
    }
    if (frequency < FREQ_RED) {
        return (color){0, frequency - FREQ_NEAR_INFRARED, 0, 0};
    }
    if (frequency < FREQ_YELLOW) {
        return (color){0, 1, frequency - FREQ_RED, 0};
    }
    if (frequency < FREQ_GREEN) {
        return (color){0, FREQ_GREEN - frequency, 1, 0};
    }
    if (frequency < FREQ_CYAN) {
        return (color){0, 0, 1, frequency - FREQ_GREEN};
    }
    if (frequency < FREQ_BLUE) {
        return (color){0, 0, FREQ_BLUE - frequency, 1};
    }
    if (frequency < FREQ_VIOLET) {
        return (color){0, 0.5*(frequency - FREQ_BLUE), 0, 1};
    }
    if (frequency < FREQ_NEAR_ULTRAVIOLET) {
        return (color){0, 0.5*(FREQ_NEAR_ULTRAVIOLET - frequency), 0, FREQ_NEAR_ULTRAVIOLET - frequency};
    }
    return C_BLACK;
}

RayPath::RayPath(
    const quaternion& origin, const quaternion& direction,
    const std::vector<std::shared_ptr<Traceable>>& objects,
    const SkySphere& sky_sphere,
    const real& max_length, const int& max_depth,
    const real& frequency, const real& index) {
    this->start = origin;
    this->direction = direction;
    this->total_length = max_length;
    this->frequency = frequency;
    this->index = index;
    this->length = std::numeric_limits<real>::infinity();
    if (max_depth > 0) {
        quaternion closest_normal;
        std::shared_ptr<Traceable> closest_object;
        std::vector<std::shared_ptr<Traceable>>::const_iterator obj;
        for (obj = objects.begin(); obj != objects.end(); ++obj) {
            const quaternion local_origin = (*obj)->left_transform*(origin - (*obj)->location)*(*obj)->right_transform;
            const quaternion local_direction = (*obj)->left_transform*direction*(*obj)->right_transform;
            auto [distance, normal] = (*obj)->trace(local_origin, local_direction);
            if (distance < this->length) {
                this->length = distance;
                closest_normal = normal;
                closest_object = (*obj);
            }
        }
        if (this->length < max_length) {
            closest_normal = inverse(closest_object->left_transform)*closest_normal*inverse(closest_object->right_transform);
            quaternion end = origin + direction * this->length;
            this->end_amplitude = closest_object->pigment->eval(end, direction, closest_normal, frequency);
            this->reflection_weight = closest_object->reflectivity->eval(end, direction, closest_normal, frequency);
            this->refraction_weight = closest_object->transparency->eval(end, direction, closest_normal, frequency);
            if (this->reflection_weight == 0) {
                this->reflected_path = nullptr;
            } else {
                const quaternion reflected_direction = direction - 2 * dot(direction, closest_normal) * closest_normal;
                this->reflected_path = new RayPath(
                    end + reflected_direction * EPSILON,
                    reflected_direction,
                    objects, sky_sphere,
                    max_length - this->length, max_depth - 1,
                    frequency, index);
            }
            if (this->refraction_weight == 0 || this->reflection_weight == 1) {
                this->refracted_path = nullptr;
            } else {
                const real refractive_index = closest_object->ior->eval(end, direction, closest_normal, frequency);
                // Adapted from:
                // https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
                real cosi = dot(direction, closest_normal);
                real etai = index;
                real etat = refractive_index;
                quaternion n = closest_normal;
                if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -n; }
                const real eta = etai / etat;
                const real k = 1 - eta*eta * (1 - cosi*cosi);
                if (k > 0) {
                    const quaternion refracted_direction = normalize(eta * direction + (eta * cosi - sqrt(k)) * n);
                    this->refracted_path = new RayPath(
                        end + refracted_direction * EPSILON,
                        refracted_direction,
                        objects, sky_sphere,
                        max_length - this->length, max_depth - 1,
                        frequency, refractive_index / index);
                } else {
                    this->refracted_path = nullptr;
                }
            }
        } else {
            this->length = max_length;
            this->end_amplitude = sky_sphere.eval(direction, frequency);
            this->reflection_weight = 0;
            this->refraction_weight = 0;
            this->reflected_path = nullptr;
            this->refracted_path = nullptr;
        }
    } else {
        this->length = EPSILON;
        this->end_amplitude = 0;
        this->reflection_weight = 0;
        this->refraction_weight = 0;
        this->reflected_path = nullptr;
        this->refracted_path = nullptr;
    }
}

real RayPath::eval(const Density& density, const int& num_samples) const {
    std::default_random_engine generator;
    std::uniform_real_distribution<real> distribution(0.0, 1.0);

    const int trunk_samples = (int) (num_samples * this->length / this->total_length);
    real result = this->end_amplitude;
    if (this->refracted_path != nullptr) {
        result *= (1 - this->refraction_weight);
        result += this->refracted_path->eval(density, num_samples - trunk_samples) * this->refraction_weight;
    }
    if (this->reflected_path != nullptr) {
        result *= (1 - this->reflection_weight);
        result += this->reflected_path->eval(density, num_samples - trunk_samples) * this->reflection_weight;
    }
    real dt = 1.0 / (real)(trunk_samples);
    real du = this->length * dt;
    for (int i = 0; i < trunk_samples; ++i) {
        real t = (trunk_samples - i - distribution(generator)) * du;
        auto [illumination, absorption] = density.eval(this->start + t*this->direction, this->direction, this->frequency);
        result += du*illumination;
        result *= exp(-absorption*du);
    }
    return result;
}
