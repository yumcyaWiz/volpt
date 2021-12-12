#include "camera.h"
#include "integrator.h"
#include "medium.h"
#include "scene.h"

int main() {
  const uint32_t width = 512;
  const uint32_t height = 512;
  const uint32_t n_samples = 10000;
  const uint32_t max_depth = 10000;

  Image image(width, height);

  const Vec3f camera_pos = 1.5 * Vec3f(48, 24, 24);
  const Vec3f camera_lookat = Vec3f(0, 12, 0);
  // const Vec3f camera_pos = 30 * Vec3f(0, 1.5, 6);
  // const Vec3f camera_lookat = 30 * Vec3f(0, 1.5, 0);
  const Vec3f camera_forward = normalize(camera_lookat - camera_pos);
  const float FOV = 0.25 * PI;

  const auto camera =
      std::make_shared<PinholeCamera>(camera_pos, camera_forward, FOV);

  // build scene
  Scene scene;
  scene.loadObj("smoke_test.obj");
  scene.loadVDB("smoke2.vdb", 0, Vec3f(0.1f, 0.8f, 0.2f),
                Vec3f(1.0f, 1.0f, 1.0f), 40.0f);

  // scene.loadObj("cornellbox_smoke.obj");
  // scene.loadVDB("smoke.vdb", 0, Vec3f(0.1f), Vec3f(1.0f), 5.0f);
  scene.build();

  // render
  UniformSampler sampler;
  PathTracing integrator(camera, n_samples, max_depth);
  // NormalIntegrator integrator(camera);
  integrator.render(scene, sampler, image);

  // gamma correction
  image.gammaCorrection(2.2f);

  // output image
  image.writePPM("output.ppm");

  return 0;
}