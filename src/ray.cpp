#include <cmath>
#include <limits>

#include "density-tracer/ray.h"

const real EPSILON = 1e-8;

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

RayPath::RayPath(const quaternion& origin, const quaternion& direction, const std::vector<std::shared_ptr<Traceable>>& objects, const int& max_depth) {
    this->start = origin;
    this->direction = direction;
    this->length = std::numeric_limits<real>::infinity();
    if (max_depth > 0) {
        quaternion closest_normal;
        std::vector<std::shared_ptr<Traceable>>::const_iterator obj;
        for (obj = objects.begin(); obj != objects.end(); ++obj) {
            const quaternion local_origin = origin - (*obj)->location;
            const quaternion local_direction = direction;
            auto [distance, normal] = (*obj)->trace(local_origin, local_direction);
            if (distance < this->length) {
                this->length = distance;
                closest_normal = normal;
            }
        }
        if (this->length < std::numeric_limits<real>::infinity()) {
            const quaternion reflected_direction = direction - 2 * dot(direction, closest_normal) * closest_normal;
            this->child = new RayPath(origin + direction * this->length + reflected_direction * EPSILON, reflected_direction, objects, max_depth - 1);
        }
    }
}

quaternion RayPath::get_location(const real& t) const {
    if (t < this->length) {
        return this->start + t * this->direction;
    }
    return this->child->get_location(t - this->length);
}
