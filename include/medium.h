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
      // when g is small, sample direction uniformly
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
                            Vec3f& throughput, Vec3f& Le,
                            bool& terminate) const = 0;

  virtual Vec3f transmittance(const Vec3f& p1, const Vec3f& p2) const = 0;
};

class HomogeneousMedium : public Medium {
 private:
  const Vec3f sigma_a;  // absorption coefficient
  const Vec3f sigma_s;  // scattering coefficient
  const Vec3f sigma_t;  // extinction coefficient

 public:
  HomogeneousMedium(float g, Vec3f sigma_a, Vec3f sigma_s)
      : Medium(g),
        sigma_a(sigma_a),
        sigma_s(sigma_s),
        sigma_t(sigma_a + sigma_s) {}

  bool sampleMedium(Ray& ray, float distToSurface, Sampler& sampler,
                    Vec3f& throughput, Vec3f& Le,
                    bool& terminate) const override {
    // NOTE: Hero wavelength sampling with balance heuristics
    // sample wavelength
    int channel = 3 * sampler.getNext1D();
    if (channel == 3) channel--;

    // sample collision-free distance
    const float t = -std::log(std::max(1.0f - sampler.getNext1D(), 0.0f)) /
                    sigma_t[channel];

    // hit volume boundary, no collision
    if (t > distToSurface) {
      const Vec3f tr = transmittance(ray.origin, ray(distToSurface));
      throughput = tr[channel] / ((tr[0] + tr[1] + tr[2]) / 3.0f);
      return false;
    }

    const Vec3f tr = transmittance(ray.origin, ray(t));

    // sample collision event
    const float p_a = sigma_a[channel] / sigma_t[channel];
    // emission
    if (sampler.getNext1D() < p_a) {
      Le = Vec3f(0);
      const Vec3f tr_sigma_a = tr * sigma_a;
      throughput = tr_sigma_a[channel] /
                   ((tr_sigma_a[0] + tr_sigma_a[1] + tr_sigma_a[2]) / 3.0f);
      terminate = true;
    }
    // in-scattering
    else {
      // sample direction
      Vec3f wi;
      phaseFunction->sampleDirection(-ray.direction, sampler, wi);

      // advance ray, and set new direction
      ray.origin = ray(t);
      ray.direction = wi;

      const Vec3f tr_sigma_s = tr * sigma_s;
      throughput = tr_sigma_s[channel] /
                   ((tr_sigma_s[0] + tr_sigma_s[1] + tr_sigma_s[2]) / 3.0f);
      terminate = false;
    }

    return true;
  }

  Vec3f transmittance(const Vec3f& p1, const Vec3f& p2) const override {
    const float dist = length(p1 - p2);
    Vec3f ret;
    for (int i = 0; i < 3; ++i) {
      ret[i] = std::exp(-sigma_t[i] * dist);
    }
    return ret;
  }
};

#endif