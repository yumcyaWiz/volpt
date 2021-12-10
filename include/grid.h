#ifndef _GRID_H
#define _GRID_H
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

#include <filesystem>

#include "core.h"

struct AABB {
  Vec3f pMin;
  Vec3f pMax;
};

class DensityGrid {
 public:
  virtual AABB getBounds() const = 0;
  virtual float getDensity(const Vec3f& pos) const = 0;
  virtual float getMaxDensity() const = 0;
};

class OpenVDBGrid : public DensityGrid {
 private:
  openvdb::FloatGrid::Ptr gridPtr;

 public:
  OpenVDBGrid(const std::filesystem::path& filepath) {
    spdlog::info("[OpenVDBGrid] loading: {}", filepath.generic_string());

    openvdb::initialize();

    // open vdb file
    openvdb::io::File file(filepath.generic_string());
    file.open();

    // get density grid
    this->gridPtr =
        openvdb::gridPtrCast<openvdb::FloatGrid>(file.readGrid("density"));
    if (!this->gridPtr) {
      spdlog::error("[Scene] failed to load density grid");
      return;
    }
    file.close();
  }

  AABB getBounds() const override {
    AABB ret;
    const auto bbox = gridPtr->evalActiveVoxelBoundingBox();
    const auto pMin = gridPtr->indexToWorld(bbox.getStart());
    const auto pMax = gridPtr->indexToWorld(bbox.getEnd());
    for (int i = 0; i < 3; ++i) {
      ret.pMin[i] = pMin[i];
      ret.pMax[i] = pMax[i];
    }
    return ret;
  }

  float getDensity(const Vec3f& pos) const override {
    const auto gridSampler =
        openvdb::tools::GridSampler<openvdb::FloatGrid,
                                    openvdb::tools::BoxSampler>(*this->gridPtr);
    return gridSampler.wsSample(openvdb::Vec3f(pos[0], pos[1], pos[2]));
  }

  float getMaxDensity() const override {
    float minValue, maxValue;
    gridPtr->evalMinMax(minValue, maxValue);
    return maxValue;
  }
};

#endif