#ifndef _RAY_H
#define _RAY_H
#include <stack>

#include "core.h"

// forward declaration
class Medium;

class Ray {
 private:
  std::stack<const Medium*> mediums;

 public:
  Vec3f origin;
  Vec3f direction;
  static constexpr float tmin = RAY_EPS;
  mutable float tmax = std::numeric_limits<float>::max();

  Ray() {}
  Ray(const Vec3f& origin, const Vec3f& direction)
      : origin(origin), direction(direction) {}

  Vec3f operator()(float t) const { return origin + t * direction; }

  bool hasMedium() const { return !mediums.empty(); }

  const Medium* getCurrentMedium() const { return mediums.top(); }

  void pushMedium(const Medium* medium) { mediums.push(medium); }

  void popMedium() {
    if (hasMedium()) {
      mediums.pop();
    }
  }
};

#endif