#include <cmath>
#include <limits>
#include <random>

#include "density-tracer/ray.h"

const real EPSILON = 1e-8;
const real ISQRT_3 = 0.5773502691896258;

std::pair<real, quaternion> Plane::trace(const quaternion& origin, const quaternion& direction) const {
    if (origin.y > 0) {
        if (direction.y < -EPSILON) {
            return std::make_pair(-origin.y / direction.y, Q_J);
        }
        return std::make_pair(std::numeric_limits<real>::infinity(), Q_ZERO);
    }
    if (direction.y > EPSILON) {
        return std::make_pair(-origin.y / direction.y, -Q_J);
    }
    return std::make_pair(std::numeric_limits<real>::infinity(), Q_ZERO);
}

std::pair<real, quaternion> Tetrahedron::trace(const quaternion& origin, const quaternion& direction) const {
    real s = ISQRT_3;
    quaternion l = origin;
    if (l.x + l.y + l.z < 1 && l.x - l.y -l.z < 1 && l.y - l.x - l.z < 1 && l.z - l.x - l.y < 1) {
        s = -s;
    }
    real o = origin.x + origin.y + origin.z;
    real d = direction.x + direction.y + direction.z;
    real t_a;
    if (fabs(d) > EPSILON) {
        t_a = (1-o) / d;
    } else {
        t_a = -1;
    }
    if (t_a > 0) {
        quaternion i = origin + direction*t_a;
        if (i.x - i.y - i.z > 1 || i.y - i.x - i.z > 1 || i.z - i.x - i.y > 1) {
            t_a = -1;
        }
    }
    o = origin.x - origin.y - origin.z;
    d = direction.x - direction.y - direction.z;
    real t_b;
    if (fabs(d) > EPSILON) {
        t_b = (1-o) / d;
    } else {
        t_b = -1;
    }
    if (t_b > 0) {
        quaternion i = origin + direction*t_b;
        if (i.x + i.y + i.z > 1 || i.y - i.x - i.z > 1 || i.z - i.x - i.y > 1) {
            t_b = -1;
        }
    }
    o = origin.y - origin.x - origin.z;
    d = direction.y - direction.x - direction.z;
    real t_c;
    if (fabs(d) > EPSILON) {
        t_c = (1-o) / d;
    } else {
        t_c = -1;
    }
    if (t_c > 0) {
        quaternion i = origin + direction*t_c;
        if (i.x + i.y + i.z > 1 || i.x - i.y - i.z > 1 || i.z - i.x - i.y > 1) {
            t_c = -1;
        }
    }
    o = origin.z - origin.y - origin.x;
    d = direction.z - direction.y - direction.x;
    real t_d;
    if (fabs(d) > EPSILON) {
        t_d = (1-o) / d;
    } else {
        t_d = -1;
    }
    if (t_d > 0) {
        quaternion i = origin + direction*t_d;
        if (i.x + i.y + i.z > 1 || i.y - i.x - i.z > 1 || i.x - i.z - i.y > 1) {
            t_d = -1;
        }
    }
    if (t_a > 0 && (t_a < t_b || t_b < 0)) {
        if (t_a < t_c || t_c < 0) {
            if (t_a < t_d || t_d < 0) {
                return std::make_pair(t_a, (quaternion){0, s, s, s});
            }
        }
    }
    if (t_b > 0 && (t_b < t_c || t_c < 0)) {
        if (t_b < t_d || t_d < 0) {
            return std::make_pair(t_b, (quaternion){0, s, -s, -s});
        }
    }
    if (t_c > 0 && (t_c < t_d || t_d < 0)) {
        return std::make_pair(t_c, (quaternion){0, -s, s, -s});
    }
    if (t_d > 0) {
        return std::make_pair(t_d, (quaternion){0, -s, -s, s});
    }
    return std::make_pair(std::numeric_limits<real>::infinity(), Q_ZERO);
}

std::pair<real, quaternion> Ball::trace(const quaternion& origin, const quaternion& direction) const {
    real o2 = origin.x*origin.x + origin.y*origin.y + origin.z*origin.z;
    real od = origin.x*direction.x + origin.y*direction.y + origin.z*direction.z;
    real d2 = direction.x*direction.x + direction.y*direction.y + direction.z*direction.z;

    real discriminant = od*od - (o2-1)*d2;
    if (discriminant < 0) {
        return std::make_pair(std::numeric_limits<real>::infinity(), Q_ZERO);
    }
    discriminant = sqrt(discriminant);
    real distance = -od - discriminant;
    if (distance < 0) {
        distance = discriminant - od;
    }
    if (distance > 0) {
        quaternion intersection = origin + distance * direction;
        return std::make_pair(distance, intersection);
    }
    return std::make_pair(std::numeric_limits<real>::infinity(), Q_ZERO);
}

RayPath::RayPath(const quaternion& origin, const quaternion& direction, const std::vector<std::shared_ptr<Traceable>>& objects, const SkySphere& sky_sphere, const real& max_length, const int& max_depth) {
    this->start = origin;
    this->direction = direction;
    this->length = std::numeric_limits<real>::infinity();
    this->total_length = max_length;
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
            this->end_color = closest_object->pigment->eval(end, direction, closest_normal);
            this->end_alpha = closest_object->reflectivity;
            if (closest_object->reflectivity != 0) {
                const quaternion reflected_direction = direction - 2 * dot(direction, closest_normal) * closest_normal;
                this->child = new RayPath(end + reflected_direction * EPSILON, reflected_direction, objects, sky_sphere, max_length - this->length, max_depth - 1);
            } else {
                this->child = nullptr;
            }
        } else {
            this->length = max_length;
            this->end_color = sky_sphere.eval(direction);
            this->end_alpha = 0;
            this->child = nullptr;
        }
    } else {
        this->length = EPSILON;
        this->end_color = Q_ZERO;
        this->end_alpha = 0;
        this->child = nullptr;
    }
}

quaternion RayPath::eval(const Density& density, const int& num_samples) const {
    std::default_random_engine generator;
    std::uniform_real_distribution<real> distribution(0.0, 1.0);

    const int trunk_samples = (int) (num_samples * this->length / this->total_length);
    color result;
    if (this->child) {
        result = child->eval(density, num_samples - trunk_samples) * this->end_alpha;
        result = result + this->end_color * (1 - this->end_alpha);
    } else {
        result = this->end_color;
    }
    real dt = 1.0 / (real)(trunk_samples);
    real du = this->length * dt;
    for (int i = 0; i < trunk_samples; ++i) {
        real t = (trunk_samples - i - distribution(generator)) * du;
        auto [illumination, absorption] = density.eval(this->start + t*this->direction, this->direction);
        result = result + du*illumination;
        result = {
            0,
            result.x * exp(-absorption.x*du),
            result.y * exp(-absorption.y*du),
            result.z * exp(-absorption.z*du)
        };
    }
    return result;
}
