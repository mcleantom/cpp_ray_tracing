#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

struct Material {
    Material(const float &r, const Vec4f &a, const Vec3f &color, const float &spec) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : refractive_index(1), albedo(1,0,0,0), diffuse_color(), specular_exponent() {}
    float refractive_index;
    Vec4f albedo;
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

struct SkyBox {
    SkyBox(char const *f) : filename(f), req_comp(0) {
        unsigned char *pixmap = stbi_load(filename, &width, &height, &comp, req_comp);
        if (!pixmap || 3!=comp) {
            std::cerr << "Cannot load the skybox image." << std::endl;
            std::abort();
        }

        envmap = std::vector<Vec3f>(width*height);
        for (int j = height-1; j>=0; j--) {
            for (int i = 0; i<width; i++) {
                envmap[i+j*width] = Vec3f(pixmap[(i+j*width)*3 + 0], pixmap[(i+j*width)*3 + 1], pixmap[(i+j*width)*3 + 2])*(1/255.);
            }
        }

        stbi_image_free(pixmap);
    }
    char const *filename;
    std::vector<Vec3f> envmap;
    int width;
    int height;
    int comp;
    int req_comp;
};

struct Light {
    Light(const Vec3f &p, const float &i): position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

Vec3f refract(const Vec3f &I, const Vec3f &N, const float &refractive_index) {
    // Snells law - https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
    float cosi = - std::max(-1.f, std::min(1.f, I*N));
    float etai = 1, etat = refractive_index;
    Vec3f n = N;
    if (cosi < 0) {
        /* theta_1 must be positive, which it is if N is the normal vector that
        points towards the point source, so swap the direction.
        */
        cosi = -cosi;
        std::swap(etai, etat);
        n = -N;
    }
    float eta = etai/etat; // r = (n1/n2)
    float k = 1 - eta*eta*(1 - cosi*cosi); // cos(theta2)^2

    return k < 0 ? Vec3f(0, 0, 0) : I*eta + n*(eta * cosi - sqrtf(k));
}

Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    // Reflection vector
    // Vf = I + 2*cos(01)*N = I * 2*(I*N)*N
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

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, const SkyBox &skybox, size_t depth=0) {
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

    if (depth >= 4 || !scene_intersect(orig, dir, spheres, point, N, material)) {
        int x = (int) ((atan2(dir.z, dir.x) / (2 * M_PI) + 0.5) * skybox.width);
        int y = (int) (acos(dir.y) / M_PI * skybox.height);
        x = std::max(0, std::min(x, skybox.width-1));
        y = std::max(0, std::min(y, skybox.height-1));
        return skybox.envmap[x + y*skybox.width];
    }

    Vec3f reflect_dir = reflect(dir, N).normalize();
    Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, skybox, depth+1);

    Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();
    Vec3f refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f refract_color = cast_ray(refract_orig, refract_dir, spheres, lights, skybox, depth+1);

    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position-point).norm();

        Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;

        /*
            If a ray going from the point towards the light intersects another object, then it is in the shadow of that 
            object, and skip that light source.
        */
        if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(
            std::max(0.f, -reflect(-light_dir, N)*dir), 
            material.specular_exponent*lights[i].intensity
        );
    }

    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity*material.albedo[1] + reflect_color*material.albedo[2] + refract_color*material.albedo[3];  // Sphere colour
}

void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights, const SkyBox &skybox) {
    const int width    = 1024*4;
    const int height   = 768*4;
    const int fov      = M_PI/2.;

    std::vector<Vec3f> framebuffer(width*height); // One dimensional array of Vec3f values (r, g, b)
    
    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        #pragma omp parallel for
        for (size_t i = 0; i<width; i++) {
            float x = (2*(i + 0.5)/(float)width - 1)*tan(fov/2.)*width/(float)height;
            float y = -(2*(j+0.5)/(float)height - 1)*tan(fov/2.);
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i+j*width] = cast_ray(Vec3f(0, 0, 0), dir, spheres, lights, skybox);
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
    SkyBox    coastline("../envmap.jpg");
    Material      ivory(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3),   50.);
    Material      glass(1.5, Vec4f( 0.,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8),  125.);
    Material red_rubber(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1),   10.);
    Material     mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.);
    
    std::vector<Sphere> spheres = {
        Sphere(Vec3f(-3,    0,   -16),  2,      ivory),
        Sphere(Vec3f(-1.0, -1.5, -12),  2,      glass),
        Sphere(Vec3f( 1.5, -0.5, -18),  3, red_rubber)
        //Sphere(Vec3f( 7,    5,   -18),  4,     mirror)
    };

    std::vector<Light> lights = {
        Light(Vec3f(-20, 20,  20), 1.5),
        Light(Vec3f( 30, 50, -25), 1.8),
        Light(Vec3f( 30, 20,  30), 1.)
    };
    render(spheres, lights, coastline);
    return 0;
}