#include "Vec3.h"

class LightSource {
 public :
  inline LightSource (Vec3f p, Vec3f c,float intensity){pos = p;color = c;i=intensity;c.normalize();}
  inline Vec3f getPos() {return pos;}
  inline Vec3f getColor() {return color;}
  inline float getIntensity(){return i;}
  float getIntensityAt(Vec3f loc);
  
 private :
  Vec3f pos;
  Vec3f color;
  float i;
};
