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

  static bool isTransmitted(const Vec3f& wo, const Vec3f& wi, const Vec3f& n) {
    return (dot(wo, n) < 0) != (dot(wi, n) < 0);
  }

  static bool isEntered(const Vec3f& wi, const Vec3f& n) {
    return dot(wi, n) < 0;
  }

  // push or pop medium
  static void updateMedium(Ray& ray, const Vec3f& wi,
                           const IntersectInfo& info) {
    if (isTransmitted(-ray.direction, wi, info.surfaceInfo.shadingNormal)) {
      if (isEntered(wi, info.surfaceInfo.shadingNormal)) {
        if (info.hitPrimitive->hasMedium()) {
          ray.pushMedium(info.hitPrimitive->getMedium());
        }
      } else {
        ray.popMedium();
      }
    }
  }

  static Vec3f sampleDirectionToLight(const Scene& scene, const Vec3f& pos,
                                      Sampler& sampler, Vec3f& dir,
                                      float& dist_to_light, float& pdf) {
    // sample light
    float light_choose_prob;
    const std::shared_ptr<Light> light =
        scene.sampleLight(sampler, light_choose_prob);

    // sample point on light
    float light_pos_pdf;
    const SurfaceInfo light_surf = light->samplePoint(sampler, light_pos_pdf);

    dir = normalize(light_surf.position - pos);
    dist_to_light = length(light_surf.position - pos);

    // convert point pdf to directional pdf
    const float cos = std::abs(dot(-dir, light_surf.shadingNormal));
    pdf = light_choose_prob *
          (dist_to_light * dist_to_light / cos * light_pos_pdf);

    return light->Le(light_surf, -dir);
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
              spdlog::error("[PathIntegrator] radiance is NaN: ({}, {}, {})",
                            radiance[0], radiance[1], radiance[2]);
              continue;
            } else if (std::isinf(radiance[0]) || std::isinf(radiance[1]) ||
                       std::isinf(radiance[2])) {
              spdlog::error("[PathIntegrator] radiance is inf: ({}, {}, {})",
                            radiance[0], radiance[1], radiance[2]);
              continue;
            } else if (radiance[0] < 0 || radiance[1] < 0 || radiance[2] < 0) {
              spdlog::error("[PathIntegrator] radiance is minus: ({}, {}, {})",
                            radiance[0], radiance[1], radiance[2]);
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

// render shading normal
class NormalIntegrator : public PathIntegrator {
 public:
  NormalIntegrator(const std::shared_ptr<Camera>& camera)
      : PathIntegrator(camera, 1) {}

  Vec3f integrate(const Ray& ray_in, const Scene& scene,
                  Sampler& sampler) const override {
    IntersectInfo info;
    if (scene.intersect(ray_in, info)) {
      // return 0.5f * (ray_in.direction + 1.0f);
      // spdlog::info("{}, {}, {}", info.surfaceInfo.shadingNormal[0],
      //              info.surfaceInfo.shadingNormal[1],
      //              info.surfaceInfo.shadingNormal[2]);
      return 0.5f * (info.surfaceInfo.shadingNormal + 1.0f);
    } else {
      return Vec3f(0);
    }
  }
};

// implementation of unidirectional path tracing
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
    ray.throughput = Vec3f(1, 1, 1);

    for (uint32_t k = 0; k < maxDepth; ++k) {
      IntersectInfo info;
      if (scene.intersect(ray, info)) {
        // russian roulette
        if (k > 0) {
          const float russian_roulette_prob = std::min(
              (ray.throughput[0] + ray.throughput[1] + ray.throughput[2]) /
                  3.0f,
              1.0f);
          if (sampler.getNext1D() >= russian_roulette_prob) {
            break;
          }
          ray.throughput /= russian_roulette_prob;
        }

        // sample medium
        bool on_surface = true;
        if (ray.hasMedium()) {
          const Medium* medium = ray.getCurrentMedium();

          Vec3f pos;
          Vec3f dir;
          Vec3f throughput_medium;
          on_surface = !medium->sampleMedium(ray, info.t, sampler, pos, dir,
                                             throughput_medium);

          // advance ray
          ray.origin = pos;
          ray.direction = dir;

          // update throughput
          ray.throughput *= throughput_medium;
        }

        if (on_surface) {
          // Le
          if (info.hitPrimitive->hasAreaLight()) {
            radiance += ray.throughput *
                        info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
            break;
          }

          // sample direction by BxDF
          Vec3f dir = ray.direction;
          if (info.hitPrimitive->hasSurface()) {
            float pdf_dir;
            const Vec3f f = info.hitPrimitive->sampleBxDF(
                -ray.direction, info.surfaceInfo,
                TransportDirection::FROM_CAMERA, sampler, dir, pdf_dir);

            // update throughput
            ray.throughput *= f *
                              cosTerm(-ray.direction, dir, info.surfaceInfo,
                                      TransportDirection::FROM_CAMERA) /
                              pdf_dir;
          }

          // update ray's medium
          updateMedium(ray, dir, info);

          // update ray
          ray.origin = info.surfaceInfo.position;
          ray.direction = dir;
        }
      } else {
        // ray goes out to the sky
        break;
      }
    }

    return radiance;
  }
};

// implementation of path tracing with next event estimation
class PathTracingNEE : public PathIntegrator {
 private:
  const uint32_t maxDepth;

  static bool isVisible(const Ray& ray, const Scene& scene, const Vec3f& origin,
                        const Vec3f& dir, float dist_to_light,
                        Vec3f& transmittance) {
    // create shadow ray
    // copy ray's mediums
    Ray shadow_ray = ray;
    shadow_ray.origin = origin;
    shadow_ray.direction = dir;
    shadow_ray.tmax = dist_to_light - RAY_EPS;

    // trace shadow ray
    transmittance = Vec3f(1);
    while (true) {
      IntersectInfo shadow_info;
      if (scene.intersect(shadow_ray, shadow_info)) {
        // occuluded
        if (shadow_info.hitPrimitive->hasSurface()) {
          return false;
        }

        // update transmittance
        if (shadow_ray.hasMedium()) {
          const Medium* medium = shadow_ray.getCurrentMedium();
          transmittance *= medium->transmittance(
              shadow_ray.origin, shadow_info.surfaceInfo.position);
        }

        // update shadow ray's medium
        updateMedium(shadow_ray, shadow_ray.direction, shadow_info);

        // advance shadow ray
        shadow_ray.origin = shadow_info.surfaceInfo.position;
        shadow_ray.tmax -= shadow_info.t;
      } else {
        // update transmittance
        if (shadow_ray.hasMedium()) {
          const Medium* medium = shadow_ray.getCurrentMedium();
          transmittance *= medium->transmittance(
              shadow_ray.origin, shadow_info.surfaceInfo.position);
        }
        break;
      }
    }

    return true;
  }

 public:
  PathTracingNEE(const std::shared_ptr<Camera>& camera, uint32_t n_samples,
                 uint32_t maxDepth = 100)
      : PathIntegrator(camera, n_samples), maxDepth(maxDepth) {}

  Vec3f integrate(const Ray& ray_in, const Scene& scene,
                  Sampler& sampler) const override {
    Vec3f radiance(0);
    Ray ray = ray_in;
    ray.throughput = Vec3f(1, 1, 1);

    uint32_t depth = 0;
    while (depth < maxDepth) {
      IntersectInfo info;
      if (scene.intersect(ray, info)) {
        // russian roulette
        if (depth > 0) {
          const float russian_roulette_prob = std::min(
              (ray.throughput[0] + ray.throughput[1] + ray.throughput[2]) /
                  3.0f,
              1.0f);
          if (sampler.getNext1D() >= russian_roulette_prob) {
            break;
          }
          ray.throughput /= russian_roulette_prob;
        }

        // sample medium
        bool is_scattered = false;
        if (ray.hasMedium()) {
          const Medium* medium = ray.getCurrentMedium();

          // sample collision position, in-scattering direction
          Vec3f pos;
          Vec3f dir;
          Vec3f throughput_medium;
          is_scattered = medium->sampleMedium(ray, info.t, sampler, pos, dir,
                                              throughput_medium);

          // next event estimation
          if (is_scattered) {
            // sample direction to the light
            Vec3f dir;
            float dist_to_light;
            float pdf_dir;
            const Vec3f Le = sampleDirectionToLight(scene, pos, sampler, dir,
                                                    dist_to_light, pdf_dir);

            // test visibility
            Vec3f transmittance;
            if (isVisible(ray, scene, pos, dir, dist_to_light, transmittance)) {
              // add contribution
              const Medium* medium = ray.getCurrentMedium();
              const Vec3f f = medium->evalPhaseFunction(-ray.direction, dir);
              const Vec3f Ls = transmittance * f * Le / pdf_dir;
              // if (length(Ls) > 1.0f) {
              //   spdlog::info("firefly");
              //   spdlog::info("Ls: ({}, {}, {})", Ls[0], Ls[1], Ls[2]);
              //   spdlog::info("transmittance: ({}, {}, {})", transmittance[0],
              //                transmittance[1], transmittance[2]);
              //   spdlog::info("f: ({}, {}, {})", f[0], f[1], f[2]);
              //   spdlog::info("pdf_dir: {}", pdf_dir);
              // }
              radiance += ray.throughput * throughput_medium * Ls;
            }
          }

          // advance ray
          ray.origin = pos;
          ray.direction = dir;

          // update throughput
          ray.throughput *= throughput_medium;
        }

        bool is_reflected_or_refracted = false;
        if (!is_scattered) {
          // direct hit to the light
          if (info.hitPrimitive->hasAreaLight()) {
            if (depth == 0) {
              radiance +=
                  ray.throughput *
                  info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
            }
            break;
          }

          // next event estimation
          if (info.hitPrimitive->hasSurface()) {
            // sample direction to the light
            Vec3f dir;
            float dist_to_light;
            float pdf_dir;
            const Vec3f Le =
                sampleDirectionToLight(scene, info.surfaceInfo.position,
                                       sampler, dir, dist_to_light, pdf_dir);

            // test visibility
            Vec3f transmittance;
            if (isVisible(ray, scene, info.surfaceInfo.position, dir,
                          dist_to_light, transmittance)) {
              // add contribution
              const Vec3f f = info.hitPrimitive->evaluateBxDF(
                  -ray.direction, dir, info.surfaceInfo,
                  TransportDirection::FROM_CAMERA);
              radiance += ray.throughput * transmittance * f *
                          cosTerm(-ray.direction, dir, info.surfaceInfo,
                                  TransportDirection::FROM_CAMERA) *
                          Le / pdf_dir;
            }
          }

          // sample direction by BxDF
          Vec3f dir = ray.direction;
          if (info.hitPrimitive->hasSurface()) {
            float pdf_dir;
            const Vec3f f = info.hitPrimitive->sampleBxDF(
                -ray.direction, info.surfaceInfo,
                TransportDirection::FROM_CAMERA, sampler, dir, pdf_dir);

            // update throughput
            ray.throughput *= f *
                              cosTerm(-ray.direction, dir, info.surfaceInfo,
                                      TransportDirection::FROM_CAMERA) /
                              pdf_dir;

            is_reflected_or_refracted = true;
          }

          // update ray's medium
          updateMedium(ray, dir, info);

          // update ray
          ray.origin = info.surfaceInfo.position;
          ray.direction = dir;
        }

        // update depth
        if (is_scattered || is_reflected_or_refracted) {
          depth++;
        }
      } else {
        // ray goes out to the sky
        break;
      }
    }
    return radiance;
  }
};

#endif