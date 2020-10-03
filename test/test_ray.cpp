#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>

#include "density-tracer/ray.h"

const real EPSILON = 1e-7;

void test_color() {
    for (real x = FREQ_INFRARED; x < FREQ_ULTRAVIOLET; x += 0.25) {
        std::cout << x << frequency_to_rgb(x) << std::endl;
    }
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

int main() {
    test_color();
    test_ball();
    test_tetrahedron();
}
