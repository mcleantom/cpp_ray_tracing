#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"

struct Sphere {
    Vec3f center;
    float radius;

    Sphere(const Vec3f &c, const float &r) : center(c), radius(r) {}

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
        // https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
        Vec3f L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0 = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N) {
    float spheres_dist = std::numeric_limits<float>::max();

    for (size_t i=0; i<spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - spheres[i].center).normalize();
        }
    }

    return spheres_dist < 1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres) {
    Vec3f point, N;
    if (!scene_intersect(orig, dir, spheres, point, N)) {
        return Vec3f(0.2, 0.7, 0.8);  // Background Colour
    }
    return Vec3f(0.4, 0.4, 0.3);  // Sphere colour
}

void render(const std::vector<Sphere> &spheres) {
    const int width    = 1024;
    const int height   = 768;
    const int fov      = M_PI/2.;

    std::vector<Vec3f> framebuffer(width*height); // One dimensional array of Vec3f values (r, g, b)

    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x = (2*(i + 0.5)/(float)width - 1)*tan(fov/2.)*width/(float)height;
            float y = -(2*(j+0.5)/(float)height - 1)*tan(fov/2.);
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i+j*width] = cast_ray(Vec3f(0, 0, 0), dir, spheres);
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("../out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
    std::vector<Sphere> spheres = {
        Sphere(Vec3f(-3, 0, -16), 2),
        Sphere(Vec3f(-1.0, -1.5, -12), 2.5),
        Sphere(Vec3f(1.5, -0.5, -18), 3),
        Sphere(Vec3f(7, 5, -18), 4)
    };
    render(spheres);
    return 0;
}