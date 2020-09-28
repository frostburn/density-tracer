#include <cmath>
#include <iostream>

#include "density-tracer/julia.h"
#include "density-tracer/ppm.h"

using namespace std;

int main(int argc, char *argv[])
{
    int image_width = 1080;
    int image_height = 1080;
    int anti_alias = 4;

    int width = image_width * anti_alias;
    int height = image_height * anti_alias;
    color *pixels;
    pixels = new color [width*height];
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            real x = (2*i - width) / (real) height * 1.2 - 0.7;
            real y = (2*j - height) / (real) height * 1.2 + 0.001;
            int idx = i + j*width;

            quaternion q = {x, y, 0, 0};
            auto res = mandelbrot(q, q, Q_ONE, 2, 256);
            real r = res.first;
            quaternion grad = res.second;
            if (r > 0) {
                pixels[idx] = {0, sin(0.1*r)*0.5+0.5, cos(0.01*r)*0.5+0.5, sin(0.001*r+1)*0.5+0.5};
                pixels[idx] = pixels[idx] * (0.5*grad.w+0.5);
            } else {
                pixels[idx] = {0, grad.w*0.1+0.1, grad.x*0.1+0.1, (grad.w+grad.x)*0.1+0.2};
            }
        }
    }
    color *aa_pixels;
    aa_pixels = new color [image_width*image_height];
    for (int j = 0; j < image_height; ++j) {
        for (int i = 0; i < image_width; ++i) {
            color pixel = Q_ZERO;
            for (int sub_j = 0; sub_j < anti_alias; ++sub_j) {
                for (int sub_i = 0; sub_i < anti_alias; ++sub_i) {
                    int idx = i*anti_alias + sub_i + (j*anti_alias + sub_j) * width;
                    pixel = pixel + pixels[idx];
                }
            }
            aa_pixels[i + j*image_width] = pixel / (anti_alias*anti_alias);
        }
    }
    delete[] pixels;
    cout_ppm(aa_pixels, image_width, image_height);
    delete[] aa_pixels;
    return EXIT_SUCCESS;
}
