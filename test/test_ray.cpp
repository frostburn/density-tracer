#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>

#include "density-tracer/ray.h"
#include "density-tracer/traceable.h"

const real EPSILON = 1e-8;

class Glass : public Pigment {
public:
    real eval(const quaternion& location, const quaternion& direction, const quaternion& normal, const real& frequency) const {
        return 1.9;
    }
};

void test_color() {
    assert(frequency_to_rgb(FREQ_NEAR_INFRARED) == C_BLACK);
    assert(frequency_to_rgb(FREQ_RED) == C_RED);
    assert(frequency_to_rgb(0.5*FREQ_RED + 0.5*FREQ_YELLOW) == 0.5*C_RED + 0.5*C_YELLOW);
    assert(frequency_to_rgb(FREQ_YELLOW) == C_YELLOW);
    assert(frequency_to_rgb(FREQ_GREEN) == C_GREEN);
    assert(frequency_to_rgb(FREQ_CYAN) == C_CYAN);
    assert(frequency_to_rgb(FREQ_BLUE) == C_BLUE);
    assert((frequency_to_rgb(0.5*FREQ_BLUE + 0.5*FREQ_VIOLET) == 0.5*C_BLUE + (color){0, 0.25, 0, 0.5}));
    assert((frequency_to_rgb(FREQ_VIOLET) == (color){0, 0.5, 0, 1}));
    assert(frequency_to_rgb(FREQ_NEAR_ULTRAVIOLET) == C_BLACK);
}

void test_mirror() {
    std::vector<std::shared_ptr<Traceable>> objects;
    const PerfectMirror mirror;
    const White air;
    const Black opaque;

    std::shared_ptr<Plane> plane = std::make_shared<Plane>();
    plane->reflectivity = &mirror;
    plane->transparency = &opaque;
    plane->ior = &air;
    objects.push_back(plane);

    const quaternion origin = {0, -1, 1, 0};
    const quaternion ray_direction = {0, sqrt(0.5), -sqrt(0.5), 0};

    const RayPath path(origin, ray_direction, objects, 2*sqrt(2), 2, FREQ_BLUE);

    assert(isclose(path.end(), Q_ZERO));
    assert(isclose(path.get_reflected_path()->end(), {0, 1, 1, 0}, EPSILON));
}

void test_refraction() {
    std::vector<std::shared_ptr<Traceable>> objects;
    const Dull dull;
    const Glass glass;
    const White transparent;

    std::shared_ptr<Hexahedron> hexahedron = std::make_shared<Hexahedron>();
    hexahedron->reflectivity = &dull;
    hexahedron->transparency = &transparent;
    hexahedron->ior = &glass;
    objects.push_back(hexahedron);

    const quaternion origin = {0, -1, 2, 0};
    const quaternion ray_direction = {0, sqrt(0.5), -sqrt(0.5), 0};

    const RayPath path(origin, ray_direction, objects, 4, 3, FREQ_GREEN);

    assert(isclose(path.end(), Q_J));
    const quaternion internal_direction = path.get_refracted_path()->get_direction();
    assert(fabs(internal_direction.x) < fabs(internal_direction.y));
    assert(fabs(internal_direction.z) < EPSILON);

    assert(fabs(path.get_refracted_path()->end().y + 1) < EPSILON);

    const quaternion external_direction = path.get_refracted_path()->get_refracted_path()->get_direction();
    assert(isclose(external_direction, ray_direction));
}

int main() {
    test_color();
    test_mirror();
    test_refraction();
}
