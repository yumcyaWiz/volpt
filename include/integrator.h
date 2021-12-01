#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H

#include <omp.h>

#include <optional>

#include "camera.h"
#include "core.h"
#include "image.h"
#include "scene.h"

class Integrator {
 protected:
  const std::shared_ptr<Camera> camera;

 public:
  Integrator(const std::shared_ptr<Camera>& camera) : camera(camera) {}

  // render scene
  virtual void render(const Scene& scene, Sampler& sampler, Image& image) = 0;

  // compute cosine term
  // NOTE: need to account for the asymmetry of BSDF when photon tracing
  // https://pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/The_Path-Space_Measurement_Equation#x3-Non-symmetryDuetoShadingNormals
  // Veach, Eric. Robust Monte Carlo methods for light transport simulation.
  // Stanford University, 1998. Section 5.3
  static float cosTerm(const Vec3f& wo, const Vec3f& wi,
                       const SurfaceInfo& surfaceInfo,
                       const TransportDirection& transport_dir) {
    const float wi_ns = dot(wi, surfaceInfo.shadingNormal);
    const float wi_ng = dot(wi, surfaceInfo.geometricNormal);
    const float wo_ns = dot(wo, surfaceInfo.shadingNormal);
    const float wo_ng = dot(wo, surfaceInfo.geometricNormal);

    // prevent light leaks
    if (wi_ng * wi_ns <= 0 || wo_ng * wo_ns <= 0) {
      return 0;
    }

    if (transport_dir == TransportDirection::FROM_CAMERA) {
      return std::abs(wi_ns);
    } else if (transport_dir == TransportDirection::FROM_LIGHT) {
      return std::abs(wo_ns) * std::abs(wi_ng) / std::abs(wo_ng);
    } else {
      spdlog::error("[Integrator] invalid transport direction");
      std::exit(EXIT_FAILURE);
    }
  }
};

// abstraction of path based integrator
class PathIntegrator : public Integrator {
 private:
  // number of samples in each pixel
  const uint32_t n_samples;

 public:
  // compute radiance coming from the given ray
  virtual Vec3f integrate(const Ray& ray, const Scene& scene,
                          Sampler& sampler) const = 0;

  PathIntegrator(const std::shared_ptr<Camera>& camera, uint32_t n_samples)
      : Integrator(camera), n_samples(n_samples) {}

  void render(const Scene& scene, Sampler& sampler,
              Image& image) override final {
    const uint32_t width = image.getWidth();
    const uint32_t height = image.getHeight();

    spdlog::info("[PathIntegrator] rendering...");
#pragma omp parallel for collapse(2) schedule(dynamic, 1)
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        // init sampler for each pixel
        const std::unique_ptr<Sampler> sampler_per_pixel = sampler.clone();
        sampler_per_pixel->setSeed((sampler.getSeed() + 1) * (j + width * i));

        // warmup sampler
        for (uint32_t k = 0; k < 100; ++k) {
          sampler_per_pixel->getNext1D();
        }

        // iteration
        for (uint32_t k = 0; k < n_samples; ++k) {
          // SSAA
          const float u =
              (2.0f * (j + sampler_per_pixel->getNext1D()) - width) / height;
          const float v =
              (2.0f * (i + sampler_per_pixel->getNext1D()) - height) / height;

          Ray ray;
          float pdf;
          if (camera->sampleRay(Vec2f(u, v), *sampler_per_pixel, ray, pdf)) {
            // compute incoming radiance
            const Vec3f radiance =
                integrate(ray, scene, *sampler_per_pixel) / pdf;

            // invalid radiance check
            if (std::isnan(radiance[0]) || std::isnan(radiance[1]) ||
                std::isnan(radiance[2])) {
              spdlog::error("[PathIntegrator] radiance is NaN");
              continue;
            } else if (std::isinf(radiance[0]) || std::isinf(radiance[1]) ||
                       std::isinf(radiance[2])) {
              spdlog::error("[PathIntegrator] radiance is inf");
              continue;
            } else if (radiance[0] < 0 || radiance[1] < 0 || radiance[2] < 0) {
              spdlog::error("[PathIntegrator] radiance is minus");
              continue;
            }

            image.addPixel(i, j, radiance);
          } else {
            image.setPixel(i, j, Vec3f(0));
          }
        }
      }
    }
    spdlog::info("[PathIntegrator] done");

    // take average
    image /= Vec3f(n_samples);
  }
};

// implementation of path tracing
// NOTE: for reference purpose
class PathTracing : public PathIntegrator {
 private:
  const uint32_t maxDepth;

 public:
  PathTracing(const std::shared_ptr<Camera>& camera, uint32_t n_samples,
              uint32_t maxDepth = 100)
      : PathIntegrator(camera, n_samples), maxDepth(maxDepth) {}

  Vec3f integrate(const Ray& ray_in, const Scene& scene,
                  Sampler& sampler) const override {
    Vec3f radiance(0);
    Ray ray = ray_in;
    Vec3f throughput(1, 1, 1);

    for (uint32_t k = 0; k < maxDepth; ++k) {
      IntersectInfo info;
      if (scene.intersect(ray, info)) {
        // russian roulette
        if (k > 0) {
          const float russian_roulette_prob = std::min(
              std::max(throughput[0], std::max(throughput[1], throughput[2])),
              1.0f);
          if (sampler.getNext1D() >= russian_roulette_prob) {
            break;
          }
          throughput /= russian_roulette_prob;
        }

        // Le
        if (info.hitPrimitive->hasAreaLight()) {
          radiance += throughput *
                      info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
        }

        // sample direction by BxDF
        Vec3f dir;
        float pdf_dir;
        Vec3f f = info.hitPrimitive->sampleBxDF(
            -ray.direction, info.surfaceInfo, TransportDirection::FROM_CAMERA,
            sampler, dir, pdf_dir);

        // update throughput and ray
        throughput *= f *
                      cosTerm(-ray.direction, dir, info.surfaceInfo,
                              TransportDirection::FROM_CAMERA) /
                      pdf_dir;
        ray = Ray(info.surfaceInfo.position, dir);
      } else {
        break;
      }
    }

    return radiance;
  }
};

#endif