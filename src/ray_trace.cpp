#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"

struct Material {
    Material(const Vec2f &a, const Vec3f &color, const float &spec) :albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : albedo(1,0), diffuse_color(), specular_exponent() {}
    Vec2f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
};


struct Sphere {
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f &c, const float &r, const Material &m) : center(c), radius(r), material(m) {}

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

struct Light {
    Light(const Vec3f &p, const float &i): position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N, Material &material) {
    // Makes sure that the closest sphere gets drawn
    float spheres_dist = std::numeric_limits<float>::max();

    for (size_t i=0; i<spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }

    return spheres_dist < 1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    /*
        Phone reflection model
        Intensity_point = Ambiant + Sum(each light) { Diffusive + Reflective (specular) }

        Kd - Diffusive constant
        Ks - Specular constent
        Lm - Direction vector from the point to the light source
        N - Normal vector of the surface
        i_m,(d/s) - The diffusive and specular intensity component of the  light source
        alpha - Shininess constant
        Rm - The reflection of Lm characterized by N using:
            Rm = 2(Lm dot N)*N - Lm
        V - The direction pointing towards the viewer/camera

        Diffusive = Kd * (Lm dot N)*i_m,d
        Specular = Ks * (Rm * V)^alpha * i_m,s

        Note: Diffuse component not affected by the viewer direction. The specular term is
        only large when the viewer direction is alligned with the reflection direction.
    */
    
    Vec3f point, N;
    Material material;
    float diffuse_light_intensity = 0,
          specular_light_intensity = 0;

    if (!scene_intersect(orig, dir, spheres, point, N, material)) {
        return Vec3f(0.2, 0.7, 0.8);  // Background Colour
    }

    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir = (lights[i].position - point).normalize();
        diffuse_light_intensity += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(
            std::max(0.f, reflect(light_dir, N)*dir), 
            material.specular_exponent*lights[i].intensity
        );
    }

    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity*material.albedo[0];  // Sphere colour
}

void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    const int width    = 1024;
    const int height   = 768;
    const int fov      = M_PI/2.;

    std::vector<Vec3f> framebuffer(width*height); // One dimensional array of Vec3f values (r, g, b)

    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x = (2*(i + 0.5)/(float)width - 1)*tan(fov/2.)*width/(float)height;
            float y = -(2*(j+0.5)/(float)height - 1)*tan(fov/2.);
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i+j*width] = cast_ray(Vec3f(0, 0, 0), dir, spheres, lights);
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("../out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max); // Clip color brightness to 255
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
    Material ivory(Vec2f(0.6, 0.3), Vec3f(0.4, 0.4, 0.3), 50.);
    Material red_rubber(Vec2f(0.9, 0.1), Vec3f(0.3, 0.1, 0.1), 10.);
    
    std::vector<Sphere> spheres = {
        Sphere(Vec3f(-3, 0, -16), 2, ivory),
        Sphere(Vec3f(-1.0, -1.5, -12), 2.5, red_rubber),
        Sphere(Vec3f(1.5, -0.5, -18), 3, red_rubber),
        Sphere(Vec3f(7, 5, -18), 4, ivory)
    };

    std::vector<Light> lights = {
        Light(Vec3f(-20, 20, 20), 1.5)
    };
    render(spheres, lights);
    return 0;
}