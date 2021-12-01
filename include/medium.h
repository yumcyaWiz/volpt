#ifndef _MEDIUM_H
#define _MEDIUM_H
#include <memory>

#include "core.h"
#include "sampler.h"

class PhaseFunction {
 public:
  // evaluate phase function
  virtual float evaluate(const Vec3f& wo, const Vec3f& wi) const = 0;

  // sample incident direction
  // its pdf is propotional to the shape of phase function
  virtual float sampleDirection(const Vec3f& wo, Sampler& sampler,
                                Vec3f& wi) const = 0;
};

// https://pbr-book.org/3ed-2018/Volume_Scattering/Phase_Functions#PhaseHG
class HenyeyGreenstein : public PhaseFunction {
 private:
  const float g;

 public:
  HenyeyGreenstein(float g) : g(g) {}

  float evaluate(const Vec3f& wo, const Vec3f& wi) const override {
    const float cosTheta = dot(wo, wi);
    const float denom = 1 + g * g + 2 * g * cosTheta;
    return PI_MUL_4_INV * (1 - g * g) / (denom * std::sqrt(denom));
  }

  // https://pbr-book.org/3ed-2018/Light_Transport_II_Volume_Rendering/Sampling_Volume_Scattering#SamplingPhaseFunctions
  float sampleDirection(const Vec3f& wo, Sampler& sampler,
                        Vec3f& wi) const override {
    const Vec2 u = sampler.getNext2D();

    // sample cosTheta
    float cosTheta;
    if (std::abs(g) < 1e-3) {
      cosTheta = 1 - 2 * u[0];
    } else {
      const float sqrTerm = (1 - g * g) / (1 - g + 2 * g * u[0]);
      cosTheta = (1 + g * g - sqrTerm * sqrTerm) / (2 * g);
    }

    // compute direction
    const float sinTheta =
        std::sqrt(std::max(1.0f - cosTheta * cosTheta, 0.0f));
    const float phi = 2 * PI * u[1];
    wi = Vec3f(std::cos(phi) * sinTheta, cosTheta, std::sin(phi) * sinTheta);

    return evaluate(wo, wi);
  }
};

class Medium {
 protected:
  const std::shared_ptr<PhaseFunction> phaseFunction;

 public:
  Medium(float g) : phaseFunction(std::make_shared<HenyeyGreenstein>(g)) {}

  // false means there is no collision
  virtual bool sampleMedium(Ray& ray, float distToSurface, Sampler& sampler,
                            Vec3f& Le, bool& terminate) const = 0;
};

class HomogeneousMedium : public Medium {
 private:
  const float sigma_a;  // absorption coefficient
  const float sigma_s;  // scattering coefficient
  const float sigma_t;  // extinction coefficient

 public:
  HomogeneousMedium(float g, float sigma_a, float sigma_s)
      : Medium(g),
        sigma_a(sigma_a),
        sigma_s(sigma_s),
        sigma_t(sigma_a + sigma_s) {}

  bool sampleMedium(Ray& ray, float distToSurface, Sampler& sampler, Vec3f& Le,
                    bool& terminate) const override {
    // sample collision-free distance
    const float t =
        -std::log(std::max(1.0f - sampler.getNext1D(), 0.0f)) / sigma_t;

    // hit volume boundary, no collision
    if (t > distToSurface) {
      return false;
    }

    const float p_a = sigma_a / sigma_t;
    // absorption/emission
    if (sampler.getNext1D() < p_a) {
      Le = Vec3f(0);
      terminate = true;
    }
    // scattering
    else {
      // sample direction
      Vec3f wi;
      phaseFunction->sampleDirection(-ray.direction, sampler, wi);

      // advance ray, and set new direction
      ray.origin = ray(t);
      ray.direction = wi;
    }

    return true;
  }
};

#endif