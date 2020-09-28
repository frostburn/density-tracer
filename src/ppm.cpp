#include <iostream>
#include <iomanip>
#include <random>

#include "density-tracer/ppm.h"

void cout_ppm(const color* pixels, const int& width, const int& height) {
    int color_depth = 255;
    std::cout << "P3" << std::endl;
    std::cout << width << " " << height << std::endl;
    std::cout << color_depth << std::endl;

    std::default_random_engine generator;
    std::uniform_real_distribution<real> distribution(0.0, 1.0/(real)color_depth);

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            color pixel = pixels[i + j*width];
            color dither = {0, distribution(generator), distribution(generator), distribution(generator)};
            pixel = clip_color(pixel + dither);
            std::cout << std::setw(4) << std::left << (int) (pixel.x * color_depth);
            std::cout << std::setw(4) << std::left << (int) (pixel.y * color_depth);
            std::cout << std::setw(4) << std::left << (int) (pixel.z * color_depth);
        }
        std::cout << std::endl;
    }
}
