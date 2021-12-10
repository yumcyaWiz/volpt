#ifndef _GRID_H
#define _GRID_H
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

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
  const openvdb::FloatGrid::Ptr gridPtr;
  const openvdb::tools::GridSampler<openvdb::FloatGrid,
                                    openvdb::tools::BoxSampler>
      gridSampler;

 public:
  OpenVDBGrid(const openvdb::FloatGrid::Ptr& gridPtr)
      : gridPtr(gridPtr), gridSampler(*gridPtr) {}

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
    return gridSampler.wsSample(openvdb::Vec3f(pos[0], pos[1], pos[2]));
  }

  float getMaxDensity() const override {
    float minValue, maxValue;
    gridPtr->evalMinMax(minValue, maxValue);
    return maxValue;
  }
};

#endif