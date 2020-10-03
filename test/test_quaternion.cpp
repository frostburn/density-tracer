#include <cassert>
#include <cmath>

#include "density-tracer/quaternion.h"

void test_sir_william_rowan_hamilton_on_16th_of_october_1843() {
    assert(Q_I*Q_I == Q_J*Q_J);
    assert(Q_J*Q_J == Q_K*Q_K);
    assert(Q_I*Q_J*Q_K == -1);
}

void test_multiplication() {
    quaternion a = {3, 5, 7, 11};
    quaternion b = {13, 17, 19, 23};

    quaternion c = a * b;
    assert(c.w == -432);
    assert(c.x == 68);
    assert(c.y == 220);
    assert(c.z == 188);

    c = b * a;
    assert(c.w == -432);
    assert(c.x == 164);
    assert(c.y == 76);
    assert(c.z == 236);
}

void test_hurwitz() {
    assert(pow((quaternion){0.5, 0.5, 0.5, 0.5}, 6) == 1);
    assert(pow((quaternion){-0.5, 0.5, 0.5, 0.5}, 3) == 1);
    assert(pow(Q_I, -4) == 1);
}

void test_eulers_identity() {
    assert(isclose(Q_ONE + exp(M_PI * Q_I), Q_ZERO));
    assert(isclose(Q_ONE + exp(M_PI * Q_J), Q_ZERO));
    assert(isclose(Q_W + exp(M_PI * Q_K), Q_ZERO));
}

void test_cross() {
    assert(cross(Q_I, Q_J) == Q_K);
}

void test_rotate() {
    assert(isclose(rotate(2*Q_I, 3*Q_K, 0.5*M_PI), 2*Q_J));
}

int main() {
    test_sir_william_rowan_hamilton_on_16th_of_october_1843();
    test_multiplication();
    test_hurwitz();
    test_eulers_identity();
    test_cross();
    test_rotate();
}
