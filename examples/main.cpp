#include "camera.h"
#include "integrator.h"
#include "scene.h"

int main() {
  const uint32_t width = 512;
  const uint32_t height = 512;
  const uint32_t n_samples = 100;
  const uint32_t max_depth = 100;

  Image image(width, height);

  const Vec3f camPos = Vec3f(0, 1, 6);
  const Vec3f camForward = Vec3f(0, 0, -1);
  const float FOV = 0.25 * PI;

  const auto camera = std::make_shared<PinholeCamera>(camPos, camForward, FOV);
  // const auto camera = std::make_shared<ThinLensCamera>(
  //     Vec3f(0, 1, 6), Vec3f(0, 0, -1), 0.38f * PI, 32.0f);
  // camera->focus(Vec3f(-0.2496, -0.001, 0.6));

  // build scene
  Scene scene;
  scene.loadModel("CornellBox-Original.obj");
  scene.build();

  // render
  UniformSampler sampler;
  PathTracing integrator(camera, n_samples);
  integrator.render(scene, sampler, image);

  // gamma correction
  image.gammaCorrection(2.2f);

  // output image
  image.writePPM("output.ppm");

  return 0;
}