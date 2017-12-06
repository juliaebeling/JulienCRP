#ifndef CRPROPA_DENSITY_H
#define CRPROPA_DENSITY_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"

namespace crpropa {
/**
 @class density
 @brief Abstract base class for Targetdensity
 */
class Density: public Referenced {
 public:
  virtual ~Density() {
  }
  virtual double getDensity(const Vector3d &position) const {
    return 0;
  };
  virtual double getDensity_HE(const Vector3d &position) const {
    return getDensity(&position)*0.11;		// H-He ratio
  };

} //namespace crpropa

#endif //CRPROPA_DENSITY_H
