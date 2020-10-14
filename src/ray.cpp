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

// Adapted from:
// https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
// XXX: This only works to and from air.
// XXX: RayPaths could use properties for polarization and wavelength.
real Fresnel::eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& ior, const real& frequency) const {
    real cosi = dot(direction, normal);
    real etai = 1.0;
    real etat = ior;
    if (cosi > 0) {
        std::swap(etai, etat);
    }
    // Compute sin using Snell's law
    real sint = etai / etat * sqrt(std::max(0.0, 1.0 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1) {
        return 1.0;
    }
    else {
        real cost = sqrt(std::max(0.0, 1 - sint * sint));
        cosi = fabs(cosi);
        real Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        real Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        return (Rs * Rs + Rp * Rp) * 0.5;
    }
}

RayPath::RayPath(
    const quaternion& origin, const quaternion& direction,
    const std::vector<std::shared_ptr<Traceable>>& objects,
    const real& max_length, const int& max_depth,
    const real& frequency) {
    this->start = origin;
    this->direction = direction;
    this->path_length = max_length;
    this->frequency = frequency;
    this->length = std::numeric_limits<real>::infinity();
    if (max_depth > 0) {
        quaternion closest_normal;
        std::shared_ptr<Traceable> closest_object;
        std::vector<std::shared_ptr<Traceable>>::const_iterator obj;
        for (obj = objects.begin(); obj != objects.end(); ++obj) {
            auto [distance, normal] = (*obj)->trace(origin, direction);
            if (distance < this->length) {
                this->length = distance;
                closest_normal = normal;
                closest_object = (*obj);
            }
        }
        if (this->length < max_length) {
            quaternion end = origin + direction * this->length;
            this->end_object = closest_object;
            this->end_normal = closest_normal;
            const real ior = closest_object->ior->eval(end, direction, closest_normal, frequency);
            this->reflection_weight = closest_object->reflectivity->eval(end, direction, closest_normal, ior, frequency);
            this->refraction_weight = closest_object->transparency->eval(end, direction, closest_normal, frequency);
            if (this->reflection_weight == 0 || max_depth == 1) {
                this->reflected_path = nullptr;
            } else {
                const quaternion reflected_direction = direction - 2 * dot(direction, closest_normal) * closest_normal;
                this->reflected_path = new RayPath(
                    end + reflected_direction * EPSILON,
                    reflected_direction,
                    objects,
                    max_length - this->length, max_depth - 1,
                    frequency);
            }
            if (this->refraction_weight == 0 || this->reflection_weight == 1 || max_depth == 1) {
                this->refracted_path = nullptr;
            } else {
                // Adapted from:
                // https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
                // XXX: Doesn't work if objects contain other objects.
                real cosi = dot(direction, closest_normal);
                real etai = 1.0;
                real etat = ior;
                quaternion n = closest_normal;
                if (cosi < 0) {
                    cosi = -cosi;
                } else {
                    std::swap(etai, etat);
                    n = -n;
                }
                const real eta = etai / etat;
                const real k = 1 - eta*eta * (1 - cosi*cosi);
                if (k > 0) {
                    const quaternion refracted_direction = normalize(eta * direction + (eta * cosi - sqrt(k)) * n);
                    this->refracted_path = new RayPath(
                        end + refracted_direction * EPSILON,
                        refracted_direction,
                        objects,
                        max_length - this->length, max_depth - 1,
                        frequency);
                } else {
                    this->refracted_path = nullptr;
                }
            }
        } else {
            this->length = max_length;
            this->end_object = nullptr;
            this->reflection_weight = 0;
            this->refraction_weight = 0;
            this->reflected_path = nullptr;
            this->refracted_path = nullptr;
        }
    } else {
        this->length = EPSILON;
        this->end_object = nullptr;
        this->reflection_weight = 0;
        this->refraction_weight = 0;
        this->reflected_path = nullptr;
        this->refracted_path = nullptr;
    }
}

RayPath::~RayPath() {
    delete this->reflected_path;
    delete this->refracted_path;
}

real RayPath::eval(const Density& density, const SkySphere& sky_sphere, const int& num_samples) const {
    std::default_random_engine generator;
    std::uniform_real_distribution<real> distribution(0.0, 1.0);

    const int trunk_samples = (int) (num_samples * this->length / this->path_length);
    const quaternion end = this->start + this->direction * this->length;
    real result;
    if (this->end_object == nullptr) {
        result = sky_sphere.eval(this->direction, this->frequency);
    } else {
        result = this->end_object->pigment->eval(end, this->direction, this->end_normal, this->frequency);
    }
    if (this->refracted_path != nullptr) {
        result *= (1 - this->refraction_weight);
        result += this->refracted_path->eval(density, sky_sphere, num_samples - trunk_samples) * this->refraction_weight;
    }
    if (this->reflected_path != nullptr) {
        result *= (1 - this->reflection_weight);
        result += this->reflected_path->eval(density, sky_sphere, num_samples - trunk_samples) * this->reflection_weight;
    }
    if (this->end_object != nullptr) {
        result += this->end_object->finish->eval(end, this->direction, this->end_normal, this->frequency);
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

quaternion RayPath::end() const {
    return this->start + this->direction * this->length;
}

quaternion RayPath::get_direction() const {
    return this->direction;
}

const RayPath* RayPath::get_reflected_path() const {
    return this->reflected_path;
}

const RayPath* RayPath::get_refracted_path() const {
    return this->refracted_path;
}
