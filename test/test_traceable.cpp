#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>

#include "density-tracer/traceable.h"

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

void test_tetrahedron() {
    quaternion origin = {0, 2, 2, 2};
    quaternion direction = normalize(-origin);
    Tetrahedron tetra;
    auto [distance, normal] = tetra.trace(origin, direction);
    assert(distance < 3);
    assert(isclose(normal, normalize(origin)));

    direction = {0, 1, 0, 0};
    auto [d, n] = tetra.trace(origin, direction);
    assert(d == std::numeric_limits<real>::infinity());

    origin = {0, 3, 3, 3};
    direction = normalize(-origin);
    auto [d2, n2] = tetra.trace(origin, direction);
    assert(d2 < 5);
    assert(isclose(n2, normalize(origin)));
}

void test_plane() {
    quaternion origin = {0, 1, 2, 3};
    quaternion direction = normalize(-origin);
    Plane plane;
    auto [distance, normal] = plane.trace(origin, direction);
    assert(fabs(distance - norm(origin)) < EPSILON);
    assert(normal == Q_J);
}

int main() {
    test_ball();
    test_tetrahedron();
    test_plane();
}
