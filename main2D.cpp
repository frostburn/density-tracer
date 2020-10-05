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
    const int anti_alias = 9;
    const int image_width = 108 * scale;
    const int image_height = 108 * scale;

    int width = image_width * anti_alias;
    int height = image_height * anti_alias;
    color *pixels;
    pixels = new color [width*height];
    int progress = 0;

    MultiBranchMandelbrot brot(14, 3, 6, 0, false, 1<<12, potential, min_r, numeric_limits<real>::infinity());

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

            // quaternion q = {-0.05*view_y-0.5, 0.05*view_x, 0, 0};
            quaternion q = {-0.65+0.5*view_y, 0.5*view_x+0.001, 0, 0};

            auto [inside, outside] = brot.eval(q, {0.52, -0.4, 0, 0});

            // const real residue = -min(0.0, outside);
            // outside = max(0.0, outside);

            // pixels[idx] = {0, rcos(outside) + rcos(residue*0.1 + 0.1*residue*residue), rcos(outside*2.2 + 0.2) + rcos(residue*0.45 - residue*residue*0.011)*(0.9 + 0.1*sin(residue)), rcos(outside*2.9) + rcos(residue*0.393 + residue*residue*0.011)*(0.3 + 0.7*sin(residue+2))};
            pixels[idx] = {0, 15*outside, 24*outside, 9*outside};
        }
    }

    color *aa_pixels = downscale(pixels, image_width, image_height, anti_alias);
    delete[] pixels;
    cout_ppm(aa_pixels, image_width, image_height);
    delete[] aa_pixels;
    return EXIT_SUCCESS;
}
