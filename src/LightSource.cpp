#include "LightSource.h"

float LightSource::getIntensityAt(Vec3f loc){
  return i/(pos-loc).squaredLength();
}
