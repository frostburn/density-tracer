#include <cmath>
#include <iostream>
#include <random>

#include "density-tracer/julia.h"
#include "density-tracer/ppm.h"
#include "density-tracer/image.h"

using namespace std;

real potential(const real& r2, const int& remaining_iterations, const real& inverse_log_exponent) {
    const real v = shifted_potential(r2, remaining_iterations, inverse_log_exponent);
    return (v - 1)*(v-1);
}

real rcos(const real& t) {
    return 0.5  - 0.5*cos(t);
}

int main(int argc, char *argv[])
{
    const int scale = 10;
    const int anti_alias = 5;
    const int image_width = 108 * scale;
    const int image_height = 108 * scale;

    int width = image_width * anti_alias;
    int height = image_height * anti_alias;
    color *pixels;
    pixels = new color [width*height];
    int progress = 0;

    NonEscapingMultiBranch brot(-5, 2, 7, max_q);

    #pragma omp parallel for
    for (int j = 0; j < height; ++j) {
        #pragma omp atomic update
            progress++;
        if (progress % 10 == 0) {
            cerr << (100 * progress) / height << "%" << endl;
        }
        for (int i = 0; i < width; ++i) {
            int idx = i + j*width;
            real view_x = (2*i - width) / (real) height * 1.0;
            real view_y = (2*j - height) / (real) height * 1.0;

            const quaternion q = {-2*view_y, 2*view_x+0.001, 0, 0};
            const quaternion c = {0.123, -0.2, 0.4, -0.3211};
            const quaternion res = brot.eval(q, c);

            const real r = pow(norm2(res-c), 0.2)-2;

            pixels[idx] = {0, 0.6*r, 0.1*r, 0.01*r};
        }
    }

    color *aa_pixels = downscale(pixels, image_width, image_height, anti_alias);
    delete[] pixels;
    cout_ppm(aa_pixels, image_width, image_height);
    delete[] aa_pixels;
    return EXIT_SUCCESS;
}
