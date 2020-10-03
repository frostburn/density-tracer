#include "density-tracer/image.h"

color* downscale(color *pixels, int result_height, int result_width, int factor) {
    color *result;
    result = new color [result_width*result_height];

    real scale = 1.0 / (real)(factor*factor);

    #pragma omp parallel for
    for (int j = 0; j < result_height; ++j) {
        for (int i = 0; i < result_width; ++i) {
            color pixel = C_BLACK;
            for (int sub_j = 0; sub_j < factor; ++sub_j) {
                for (int sub_i = 0; sub_i < factor; ++sub_i) {
                    int idx = i*factor + sub_i + (j*factor + sub_j) * (result_width*factor);
                    pixel = pixel + clip_color(pixels[idx]);
                }
            }
            result[i + j*result_width] = pixel * scale;
        }
    }
    return result;
}
