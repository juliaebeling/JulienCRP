#ifndef CRPROPA_DENSITY_H
#define CRPROPA_DENSITY_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"

namespace crpropa {
/**
 @class Density
 @brief Abstract base class for Targetdensity
 */
class Density: public Referenced {
 public:
  virtual ~Density() {
  }
  virtual double getDensity(const Vector3d &position) const {
    return 0;
  };
  
};

}//namespace crpropa

#endif //CRPROPA_DENSITY_H
