#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>

#include "density-tracer/ray.h"

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

int main() {
    test_color();
}
