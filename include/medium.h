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
  virtual bool sampleMedium(const Ray& ray, float distToSurface,
                            Sampler& sampler, Vec3f& pos, Vec3f& dir,
                            Vec3f& throughput) const = 0;

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

  // NOTE: ignore emission
  bool sampleMedium(const Ray& ray, float distToSurface, Sampler& sampler,
                    Vec3f& pos, Vec3f& dir, Vec3f& throughput) const override {
    // sample wavelength by throughput * single scattering albedo

    // Wrenninge, Magnus, Ryusuke Villemin, and Christophe Hery. Path traced
    // subsurface scattering using anisotropic phase functions and
    // non-exponential free flights. Tech. Rep. 17-07, Pixar. https://graphics.
    // pixar. com/library/PathTracedSubsurface, 2017.
    const Vec3f throughput_albedo = ray.throughput * (sigma_s / sigma_t);
    DiscreteEmpiricalDistribution1D distribution(throughput_albedo.getPtr(), 3);
    const Vec3f pdf_wavelength(distribution.getPDF(0), distribution.getPDF(1),
                               distribution.getPDF(2));
    float _pdf;
    const uint32_t channel = distribution.sample(sampler.getNext1D(), _pdf);

    // sample collision-free distance
    const float t = -std::log(std::max(1.0f - sampler.getNext1D(), 0.0f)) /
                    sigma_t[channel];

    // hit volume boundary, no collision
    if (t > distToSurface - RAY_EPS) {
      pos = ray(distToSurface);
      dir = ray.direction;

      const Vec3f tr = transmittance(ray.origin, pos);
      const Vec3f p_surface = tr;
      const Vec3f pdf = pdf_wavelength * p_surface;
      throughput = tr / (pdf[0] + pdf[1] + pdf[2]);
      return false;
    }

    // in-scattering
    // sample direction
    phaseFunction->sampleDirection(-ray.direction, sampler, dir);

    pos = ray(t);
    const Vec3f tr = transmittance(ray.origin, pos);
    const Vec3f pdf_distance = sigma_t * tr;
    const Vec3f pdf = pdf_wavelength * pdf_distance;
    throughput = (tr * sigma_s) / (pdf[0] + pdf[1] + pdf[2]);

    return true;
  }

  Vec3f transmittance(const Vec3f& p1, const Vec3f& p2) const override {
    const float dist = length(p1 - p2);
    return exp(-sigma_t * dist);
  }
};

class HomogeneousMediumAchromatic : public Medium {
 private:
  const float sigma_a;  // absorption coefficient
  const float sigma_s;  // scattering coefficient
  const float sigma_t;  // extinction coefficient

 public:
  HomogeneousMediumAchromatic(float g, float sigma_a, float sigma_s)
      : Medium(g),
        sigma_a(sigma_a),
        sigma_s(sigma_s),
        sigma_t(sigma_a + sigma_s) {}

  // NOTE: ignore emission
  bool sampleMedium(const Ray& ray, float distToSurface, Sampler& sampler,
                    Vec3f& pos, Vec3f& dir, Vec3f& throughput) const override {
    // sample collision-free distance
    const float t =
        -std::log(std::max(1.0f - sampler.getNext1D(), 0.0f)) / sigma_t;

    // hit volume boundary, no collision
    if (t > distToSurface - RAY_EPS) {
      pos = ray(distToSurface);
      dir = ray.direction;
      throughput = Vec3f(1);
      return false;
    }

    // in-scattering
    // sample direction
    phaseFunction->sampleDirection(-ray.direction, sampler, dir);

    pos = ray(t);
    throughput = Vec3f(sigma_s / sigma_t);

    return true;
  }

  Vec3f transmittance(const Vec3f& p1, const Vec3f& p2) const override {
    const float dist = length(p1 - p2);
    return exp(-Vec3f(sigma_t) * dist);
  }
};

class HomogeneousMediumMIS : public Medium {
 private:
  const Vec3f sigma_a;  // absorption coefficient
  const Vec3f sigma_s;  // scattering coefficient
  const Vec3f sigma_t;  // extinction coefficient

 public:
  HomogeneousMediumMIS(float g, Vec3f sigma_a, Vec3f sigma_s)
      : Medium(g),
        sigma_a(sigma_a),
        sigma_s(sigma_s),
        sigma_t(sigma_a + sigma_s) {}

  // NOTE: ignore emission
  bool sampleMedium(const Ray& ray, float distToSurface, Sampler& sampler,
                    Vec3f& pos, Vec3f& dir, Vec3f& throughput) const override {
    // NOTE: Hero wavelength sampling with balance heuristics
    // sample wavelength
    int channel = 3 * sampler.getNext1D();
    if (channel == 3) channel--;

    // sample collision-free distance
    const float t = -std::log(std::max(1.0f - sampler.getNext1D(), 0.0f)) /
                    sigma_t[channel];

    // hit volume boundary, no collision
    if (t > distToSurface - RAY_EPS) {
      pos = ray(distToSurface);
      dir = ray.direction;
      const Vec3f tr = transmittance(ray.origin, pos);
      throughput = tr / ((tr[0] + tr[1] + tr[2]) / 3.0f);
      return false;
    }

    // in-scattering
    // sample direction
    phaseFunction->sampleDirection(-ray.direction, sampler, dir);

    pos = ray(t);
    const Vec3f tr = transmittance(ray.origin, pos);
    const Vec3f tr_sigma_t = tr * sigma_t;
    throughput = (tr * sigma_s) /
                 ((tr_sigma_t[0] + tr_sigma_t[1] + tr_sigma_t[2]) / 3.0f);

    return true;
  }

  Vec3f transmittance(const Vec3f& p1, const Vec3f& p2) const override {
    const float dist = length(p1 - p2);
    return exp(-sigma_t * dist);
  }
};

// no MIS version
class HomogeneousMediumNoMIS : public Medium {
 private:
  const Vec3f sigma_a;  // absorption coefficient
  const Vec3f sigma_s;  // scattering coefficient
  const Vec3f sigma_t;  // extinction coefficient

 public:
  HomogeneousMediumNoMIS(float g, Vec3f sigma_a, Vec3f sigma_s)
      : Medium(g),
        sigma_a(sigma_a),
        sigma_s(sigma_s),
        sigma_t(sigma_a + sigma_s) {}

  // NOTE: ignore emission
  bool sampleMedium(const Ray& ray, float distToSurface, Sampler& sampler,
                    Vec3f& pos, Vec3f& dir, Vec3f& throughput) const override {
    // sample wavelength
    int channel = 3 * sampler.getNext1D();
    if (channel == 3) channel--;
    const float pdf_wavelength = 1.0f / 3.0f;

    // sample collision-free distance
    const float t = -std::log(std::max(1.0f - sampler.getNext1D(), 0.0f)) /
                    sigma_t[channel];
    const float pdf_distance =
        sigma_t[channel] * std::exp(-sigma_t[channel] * t);

    // hit volume boundary, no collision
    if (t > distToSurface - RAY_EPS) {
      pos = ray(distToSurface);
      dir = ray.direction;
      const Vec3f tr = transmittance(ray.origin, pos);
      const float p_surface = std::exp(-sigma_t[channel] * distToSurface);
      throughput = 1.0f / 3.0f * tr / (pdf_wavelength * p_surface);
      return false;
    }

    // in-scattering
    // sample direction
    phaseFunction->sampleDirection(-ray.direction, sampler, dir);

    pos = ray(t);
    throughput = 1.0f / 3.0f * transmittance(ray.origin, pos) * sigma_s /
                 (pdf_wavelength * pdf_distance);

    return true;
  }

  Vec3f transmittance(const Vec3f& p1, const Vec3f& p2) const override {
    const float dist = length(p1 - p2);
    return exp(-sigma_t * dist);
  }
};

#endif