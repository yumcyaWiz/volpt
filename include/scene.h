#ifndef _SCENE_H
#define _SCENE_H
#include <embree3/rtcore.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "core.h"
#include "grid.h"
#include "medium.h"
#include "primitive.h"
#include "tiny_obj_loader.h"

const std::shared_ptr<BxDF> createDefaultBxDF() {
  return std::make_shared<Lambert>(Vec3(0.9f));
}

// create BxDF from tinyobj material
const std::shared_ptr<BxDF> createBxDF(const tinyobj::material_t& material) {
  if (material.unknown_parameter.count("no_surface") == 1) {
    return nullptr;
  }

  const Vec3f kd =
      Vec3f(material.diffuse[0], material.diffuse[1], material.diffuse[2]);
  const Vec3f ks =
      Vec3f(material.specular[0], material.specular[1], material.specular[2]);

  switch (material.illum) {
    case 5:
      // mirror
      return std::make_shared<Mirror>(Vec3(1.0f));
    case 7:
      // glass
      return std::make_shared<Glass>(Vec3(1.0f), material.ior);
    default:
      // lambert
      return std::make_shared<Lambert>(kd);
  }
}

const std::shared_ptr<Medium> createDefaultMedium() { return nullptr; }

inline Vec3f parseVec3(const std::string& str) {
  // split string by space
  std::vector<std::string> tokens;
  std::stringstream ss(str);
  std::string buf;
  while (std::getline(ss, buf, ' ')) {
    if (!buf.empty()) {
      tokens.emplace_back(buf);
    }
  }

  if (tokens.size() != 3) {
    spdlog::error("invalid vec3 string");
    std::exit(EXIT_FAILURE);
  }

  // string to float conversion
  return Vec3f(std::stof(tokens[0]), std::stof(tokens[1]),
               std::stof(tokens[2]));
}

// create Medium from tinyobj material
const std::shared_ptr<Medium> createMedium(
    const tinyobj::material_t& material) {
  if (material.unknown_parameter.count("g") == 1 &&
      material.unknown_parameter.count("medium_color") == 1 &&
      material.unknown_parameter.count("scattering_distance") == 1) {
    const float g = std::stof(material.unknown_parameter.at("g"));

    // Chiang, Matt Jen-Yuan, Peter Kutz, and Brent Burley. "Practical and
    // controllable subsurface scattering for production path tracing." ACM
    // SIGGRAPH 2016 Talks. 2016. 1-2.
    const Vec3f A = parseVec3(material.unknown_parameter.at("medium_color"));
    const Vec3f d =
        parseVec3(material.unknown_parameter.at("scattering_distance"));
    const Vec3f alpha =
        Vec3f(1.0f) - exp(-5.09406 * A + 2.61188 * A * A - 4.31805 * A * A * A);
    const Vec3f s = Vec3f(1.9) - A + Vec3f(3.5) * (A - 0.8) * (A - 0.8);
    const Vec3f sigma_t = 1.0f / (d * s);

    const Vec3f sigma_s = alpha * sigma_t;
    const Vec3f sigma_a = sigma_t - sigma_s;
    // const Vec3f sigma_s = Vec3f(0.05, 0.99, 0.5);
    // const Vec3f sigma_a = Vec3f(0.001, 0.001, 0.001);

    // spdlog::info("sigma_a: ({}, {}, {})", sigma_a[0], sigma_a[1],
    // sigma_a[2]); spdlog::info("sigma_s: ({}, {}, {})", sigma_s[0],
    // sigma_s[1], sigma_s[2]);

    return std::make_shared<HomogeneousMedium>(g, sigma_a, sigma_s);
  } else {
    return nullptr;
  }
}

// create AreaLight from tinyobj material
const std::shared_ptr<AreaLight> createAreaLight(
    const tinyobj::material_t& material, const Triangle* tri) {
  if (material.emission[0] > 0 || material.emission[1] > 0 ||
      material.emission[2] > 0) {
    const Vec3f le =
        Vec3f(material.emission[0], material.emission[1], material.emission[2]);
    return std::make_shared<AreaLight>(le, tri);
  } else {
    return nullptr;
  }
}

class Scene {
 private:
  // mesh data
  // NOTE: assuming size of normals, texcoords == size of vertices
  std::vector<float> vertices;
  std::vector<uint32_t> indices;
  std::vector<float> normals;
  std::vector<float> texcoords;

  std::vector<std::optional<tinyobj::material_t>> materials;

  // triangles
  // NOTE: per face
  std::vector<Triangle> triangles;

  // BxDFs
  // NOTE: per face
  std::vector<std::shared_ptr<BxDF>> bxdfs;

  // mediums
  // NOTE: per face
  std::vector<std::shared_ptr<Medium>> mediums;

  // lights
  // NOTE: per face
  std::vector<std::shared_ptr<Light>> lights;

  // primitives
  // NOTE: per face
  std::vector<Primitive> primitives;

  // embree
  RTCDevice device;
  RTCScene scene;

  bool hasLight(uint32_t faceID) const { return lights[faceID] != nullptr; }

  void clear() {
    vertices.clear();
    indices.clear();
    normals.clear();
    texcoords.clear();
    bxdfs.clear();

    triangles.clear();
    bxdfs.clear();
    lights.clear();
    primitives.clear();
  }

 public:
  Scene() {}
  ~Scene() {
    clear();
    rtcReleaseScene(scene);
    rtcReleaseDevice(device);
  }

  // load obj file
  // TODO: remove vertex duplication
  void loadModel(const std::filesystem::path& filepath) {
    clear();

    spdlog::info("[Scene] loading: {}", filepath.generic_string());

    tinyobj::ObjReaderConfig reader_config;
    reader_config.mtl_search_path = "./";
    reader_config.triangulate = true;

    tinyobj::ObjReader reader;
    if (!reader.ParseFromFile(filepath.generic_string(), reader_config)) {
      if (!reader.Error().empty()) {
        spdlog::error("[Scene] failed to load {} : {}",
                      filepath.generic_string(), reader.Error());
        std::exit(EXIT_FAILURE);
      }
      return;
    }

    if (!reader.Warning().empty()) {
      spdlog::warn("[Scene] {}", reader.Warning());
    }

    const auto& attrib = reader.GetAttrib();
    const auto& shapes = reader.GetShapes();
    const auto& materials = reader.GetMaterials();

    // loop over shapes
    for (size_t s = 0; s < shapes.size(); ++s) {
      size_t index_offset = 0;
      // loop over faces
      for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); ++f) {
        const size_t fv =
            static_cast<size_t>(shapes[s].mesh.num_face_vertices[f]);

        std::vector<Vec3f> vertices;
        std::vector<Vec3f> normals;
        std::vector<Vec2f> texcoords;

        // loop over vertices
        // get vertices, normals, texcoords of a triangle
        for (size_t v = 0; v < fv; ++v) {
          const tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

          const tinyobj::real_t vx =
              attrib.vertices[3 * static_cast<size_t>(idx.vertex_index) + 0];
          const tinyobj::real_t vy =
              attrib.vertices[3 * static_cast<size_t>(idx.vertex_index) + 1];
          const tinyobj::real_t vz =
              attrib.vertices[3 * static_cast<size_t>(idx.vertex_index) + 2];
          vertices.push_back(Vec3f(vx, vy, vz));

          if (idx.normal_index >= 0) {
            const tinyobj::real_t nx =
                attrib.normals[3 * static_cast<size_t>(idx.normal_index) + 0];
            const tinyobj::real_t ny =
                attrib.normals[3 * static_cast<size_t>(idx.normal_index) + 1];
            const tinyobj::real_t nz =
                attrib.normals[3 * static_cast<size_t>(idx.normal_index) + 2];
            normals.push_back(Vec3f(nx, ny, nz));
          }

          if (idx.texcoord_index >= 0) {
            const tinyobj::real_t tx =
                attrib
                    .texcoords[2 * static_cast<size_t>(idx.texcoord_index) + 0];
            const tinyobj::real_t ty =
                attrib
                    .texcoords[2 * static_cast<size_t>(idx.texcoord_index) + 1];
            texcoords.push_back(Vec2f(tx, ty));
          }
        }

        // if normals is empty, add geometric normal
        if (normals.size() == 0) {
          const Vec3f v1 = normalize(vertices[1] - vertices[0]);
          const Vec3f v2 = normalize(vertices[2] - vertices[0]);
          const Vec3f n = normalize(cross(v1, v2));
          normals.push_back(n);
          normals.push_back(n);
          normals.push_back(n);
        }

        // if texcoords is empty, add barycentric coords
        if (texcoords.size() == 0) {
          texcoords.push_back(Vec2f(0, 0));
          texcoords.push_back(Vec2f(1, 0));
          texcoords.push_back(Vec2f(0, 1));
        }

        // populate vertices, indices, normals, texcoords
        for (int i = 0; i < 3; ++i) {
          this->vertices.push_back(vertices[i][0]);
          this->vertices.push_back(vertices[i][1]);
          this->vertices.push_back(vertices[i][2]);

          this->normals.push_back(normals[i][0]);
          this->normals.push_back(normals[i][1]);
          this->normals.push_back(normals[i][2]);

          this->texcoords.push_back(texcoords[i][0]);
          this->texcoords.push_back(texcoords[i][1]);

          this->indices.push_back(this->indices.size());
        }

        // populate materials
        const int materialID = shapes[s].mesh.material_ids[f];
        std::optional<tinyobj::material_t> material = std::nullopt;
        if (materialID != -1) {
          material = materials[materialID];
        }
        this->materials.push_back(material);

        index_offset += fv;
      }
    }

    // populate  triangles
    for (size_t faceID = 0; faceID < nFaces(); ++faceID) {
      // add triangle
      this->triangles.emplace_back(this->vertices.data(), this->indices.data(),
                                   this->normals.data(), this->texcoords.data(),
                                   faceID);
    }

    // populate bxdfs, mediums
    for (size_t faceID = 0; faceID < nFaces(); ++faceID) {
      // TODO: remove duplicate
      const auto material = this->materials[faceID];
      if (material) {
        const tinyobj::material_t& m = material.value();
        this->bxdfs.push_back(createBxDF(m));
        this->mediums.push_back(createMedium(m));
      } else {
        this->bxdfs.push_back(createDefaultBxDF());
        this->mediums.push_back(createDefaultMedium());
      }
    }

    // populate lights, primitives
    for (size_t faceID = 0; faceID < nFaces(); ++faceID) {
      // add light
      std::shared_ptr<Light> light = nullptr;
      const auto material = this->materials[faceID];
      if (material) {
        const tinyobj::material_t& m = material.value();
        light = createAreaLight(m, &this->triangles[faceID]);
        if (light != nullptr) {
          lights.push_back(light);
        }
      }

      // add primitive
      primitives.emplace_back(&this->triangles[faceID],
                              this->bxdfs[faceID].get(),
                              this->mediums[faceID].get(), light.get());
    }

    spdlog::info("[Scene] vertices: {}", nVertices());
    spdlog::info("[Scene] faces: {}", nFaces());
    spdlog::info("[Scene] lights: {}", lights.size());
  }

  uint32_t nVertices() const { return vertices.size() / 3; }
  uint32_t nFaces() const { return indices.size() / 3; }

  void build() {
    spdlog::info("[Scene] building scene...");

    // setup embree
    device = rtcNewDevice(NULL);
    scene = rtcNewScene(device);

    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

    // set vertices
    float* vb = (float*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0,
                                                RTC_FORMAT_FLOAT3,
                                                3 * sizeof(float), nVertices());
    for (size_t i = 0; i < vertices.size(); ++i) {
      vb[i] = vertices[i];
    }

    // set indices
    unsigned* ib = (unsigned*)rtcSetNewGeometryBuffer(
        geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned),
        nFaces());
    for (size_t i = 0; i < indices.size(); ++i) {
      ib[i] = indices[i];
    }

    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene);

    spdlog::info("[Scene] done");
  }

  // ray-scene intersection
  bool intersect(const Ray& ray, IntersectInfo& info) const {
    RTCRayHit rayhit;
    rayhit.ray.org_x = ray.origin[0];
    rayhit.ray.org_y = ray.origin[1];
    rayhit.ray.org_z = ray.origin[2];
    rayhit.ray.dir_x = ray.direction[0];
    rayhit.ray.dir_y = ray.direction[1];
    rayhit.ray.dir_z = ray.direction[2];
    rayhit.ray.tnear = ray.tmin;
    rayhit.ray.tfar = ray.tmax;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    rtcIntersect1(scene, &context, &rayhit);

    if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
      info.t = rayhit.ray.tfar;

      // get triangle shape
      const Triangle& tri = this->triangles[rayhit.hit.primID];

      // set surface info
      info.surfaceInfo.position = ray(info.t);
      info.surfaceInfo.barycentric = Vec2f(rayhit.hit.u, rayhit.hit.v);
      info.surfaceInfo.texcoords =
          tri.getTexcoords(info.surfaceInfo.barycentric);
      info.surfaceInfo.geometricNormal = tri.getGeometricNormal();
      info.surfaceInfo.shadingNormal =
          tri.computeShadingNormal(info.surfaceInfo.barycentric);
      orthonormalBasis(info.surfaceInfo.shadingNormal, info.surfaceInfo.dpdu,
                       info.surfaceInfo.dpdv);

      // set primitive
      info.hitPrimitive = &this->primitives[rayhit.hit.primID];

      return true;
    } else {
      return false;
    }
  }

  std::shared_ptr<Light> sampleLight(Sampler& sampler, float& pdf) const {
    unsigned int lightIdx = lights.size() * sampler.getNext1D();
    if (lightIdx == lights.size()) lightIdx--;
    pdf = 1.0f / lights.size();
    return lights[lightIdx];
  }
};

#endif