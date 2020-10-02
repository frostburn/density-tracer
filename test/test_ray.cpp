#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>

#include "density-tracer/ray.h"

const real EPSILON = 1e-7;


void test_ball() {
    quaternion origin = {0, -2, 0, 0};
    quaternion direction = {0, 1, 0, 0};
    Ball ball;
    auto [distance, normal] = ball.trace(origin, direction);
    assert(fabs(distance - 1) < EPSILON);
    assert(isclose(normal, -Q_I));

    origin = {0, 1, 1, 0.5};
    direction = normalize((quaternion){0, -0.9, -0.95, -0.49});
    auto [d, n] = ball.trace(origin, direction);
    assert(d < 2);
    assert(dot(direction, n) < 0);

    direction = normalize((quaternion){0, 1, 2, 3});
    auto [d2, n2] = ball.trace(origin, direction);
    assert(d2 == std::numeric_limits<real>::infinity());

    origin = {0, 0, 2, 0};
    direction = {0, 0, 0, 1};
    auto [d3, n3] = ball.trace(origin, direction);
    assert(d3 == std::numeric_limits<real>::infinity());
}

void test_reflection() {
    std::vector<std::shared_ptr<Traceable>> objects;

    std::shared_ptr<Plane> plane = std::make_shared<Plane>();
    plane->location = {0, 0, -1, 0};
    objects.push_back(plane);

    const RayPath path({0, -1, 0, 0}, normalize({0, 1, -1, 0}), objects, 5);
    assert((path.get_location(0) == (quaternion){0, -1, 0, 0}));
    assert(isclose(path.get_location(2*sqrt(2)), {0, 1, 0, 0}, EPSILON));
}


int main() {
    test_ball();
    test_reflection();
}
