// --------------------------------------------------------------------------
// Copyright(C) 2009-2016
// Tamy Boubekeur
// 
// Permission granted to use this code only for teaching projects and 
// private practice.
//
// Do not distribute this code outside the teaching assignements.                                                                           
// All rights reserved.                                                       
// --------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <iostream>

/// Vector in 3 dimensions, with basics operators overloaded.
template <typename T>
class Vec3 {
public:
  inline Vec3 (void)  { m_p[0] = m_p[1] = m_p[2] = 0.0; }

  inline Vec3 (T p0, T p1, T p2) {
    m_p[0] = p0;
    m_p[1] = p1;
    m_p[2] = p2;
  };

  inline Vec3 (const Vec3 & v) {
    init (v[0], v[1], v[2]);
  }

  ~Vec3() {}

  inline Vec3 (T* p) {
    m_p[0] = p[0];
    m_p[1] = p[1];
    m_p[2] = p[2];
  };

  inline T& operator[] (int Index) {
    return (m_p[Index]);
  };

  inline const T& operator[] (int Index) const {
    return (m_p[Index]);
  };

  inline Vec3& operator= (const Vec3 & p) {
    m_p[0] = p[0];
    m_p[1] = p[1];
    m_p[2] = p[2];
    return (*this);
  };

  inline Vec3& operator+= (const Vec3 & p) {
    m_p[0] += p[0];
    m_p[1] += p[1];
    m_p[2] += p[2];
    return (*this);
  };

  inline Vec3& operator-= (const Vec3 & p) {
    m_p[0] -= p[0];
    m_p[1] -= p[1];
    m_p[2] -= p[2];
    return (*this);
  };

  inline Vec3& operator*= (const Vec3 & p) {
    m_p[0] *= p[0];
    m_p[1] *= p[1];
    m_p[2] *= p[2];
    return (*this);
  };

  inline Vec3& operator*= (T s) {
    m_p[0] *= s;
    m_p[1] *= s;
    m_p[2] *= s;
    return (*this);
  };

  inline Vec3& operator/= (const Vec3 & p) {
    m_p[0] /= p[0];
    m_p[1] /= p[1];
    m_p[2] /= p[2];
    return (*this);
  };

  inline Vec3& operator/= (T s) {
    m_p[0] /= s;
    m_p[1] /= s;
    m_p[2] /= s;
    return (*this);
  };

  inline Vec3 operator+ (const Vec3 & p) const {
    Vec3 res;
    res[0] = m_p[0] + p[0];
    res[1] = m_p[1] + p[1];
    res[2] = m_p[2] + p[2];
    return (res);
  };

  inline Vec3 operator- (const Vec3 & p) const {
    Vec3 res;
    res[0] = m_p[0] - p[0];
    res[1] = m_p[1] - p[1];
    res[2] = m_p[2] - p[2];
    return (res);
  };

  inline Vec3 operator- () const {
    Vec3 res;
    res[0] = -m_p[0];
    res[1] = -m_p[1];
    res[2] = -m_p[2];
    return (res);
  };

  inline Vec3 operator* (const Vec3 & p) const {
    Vec3 res;
    res[0] = m_p[0] * p[0];
    res[1] = m_p[1] * p[1];
    res[2] = m_p[2] * p[2];
    return (res);
  };

  inline Vec3 operator* (T s) const {
    Vec3 res;
    res[0] = m_p[0] * s;
    res[1] = m_p[1] * s;
    res[2] = m_p[2] * s;
    return (res);
  };

  inline Vec3 operator/ (const Vec3 & p) const {
    Vec3 res;
    res[0] = m_p[0] / p[0];
    res[1] = m_p[1] / p[1];
    res[2] = m_p[2] / p[2];
    return (res);
  };

  inline Vec3 operator/ (T s) const {
    Vec3 res;
    res[0] = m_p[0] / s;
    res[1] = m_p[1] / s;
    res[2] = m_p[2] / s;
    return (res);
  };

  inline bool operator == (const Vec3 & a) const {
    return(m_p[0] == a[0] && m_p[1] == a[1] && m_p[2] == a[2]);
  };

  inline bool operator != (const Vec3 & a) const {
    return(m_p[0] != a[0] || m_p[1] != a[1] || m_p[2] != a[2]);
  };

  inline bool operator < (const Vec3 & a) const {
    return(m_p[0] < a[0] && m_p[1] < a[1] && m_p[2] < a[2]);
  };

  inline bool operator >= (const Vec3 &a) const {
    return(m_p[0] >= a[0] && m_p[1] >= a[1] && m_p[2] >= a[2]);
  };

  inline Vec3 & init (T x, T y, T z) {
    m_p[0] = x;
    m_p[1] = y;
    m_p[2] = z;
    return (*this);
  };

  inline T squaredLength() const {
    return (dot (*this, *this));
  };

  inline T length() const {
    return (T)sqrt (squaredLength());
  };

  /// Return length after normalization
  inline T normalize () {
    T l = length ();
    if (l == 0.0f)
      return 0;
    T invL = 1.0f / l;
    m_p[0] *= invL;
    m_p[1] *= invL;
    m_p[2] *= invL;
    return l;
  };

  inline void getTwoOrthogonals (Vec3 & u, Vec3 & v) const {
    if (fabs(m_p[0]) < fabs(m_p[1])) {
      if (fabs(m_p[0]) < fabs(m_p[2])) {
                    // m_p[0] is minimum
        u = Vec3 (0, -m_p[2], m_p[1]);
      }
      else {
                    // m_p[2] is mimimum
        u = Vec3 (-m_p[1], m_p[0], 0);
      }
    }
    else {
      if (fabs(m_p[1]) < fabs(m_p[2])) {

                    // m_p[1] is minimum
        u = Vec3 (m_p[2], 0, -m_p[0]);
      }
      else {
                    // m_p[2] is mimimum
        u = Vec3(-m_p[1], m_p[0], 0);
      }
    }
    v = cross (*this, u);
  }

  inline Vec3 projectOn (const Vec3 & N, const Vec3 & P) const {
    T w = dot (((*this) - P), N);
    return (*this) - (N * w);
  }

protected:
  T m_p[3];
};

template <class T>
inline T length (const Vec3<T> & a) {
  return a.length ();
}

  template <class T>
inline T dist (const Vec3<T> & a, const Vec3<T> & b) {
  return (a-b).length ();
}

template <class T>
inline T dot (const Vec3<T> & a, const Vec3<T> & b) {
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

template <class T>
inline Vec3<T> cross (const Vec3<T> & a, const Vec3<T> & b) {
  Vec3<T> r;
  r[0] = a[1] * b[2] - a[2] * b[1];
  r[1] = a[2] * b[0] - a[0] * b[2];
  r[2] = a[0] * b[1] - a[1] * b[0];
  return r;
}

template <class T>
inline Vec3<T> normalize (const Vec3<T> & x) {
  Vec3<T> n (x);
  n.normalize ();
  return n;
}

template <class T>
inline Vec3<T> mix (const Vec3<T> & u, const Vec3<T> & v, float alpha) {
  return (u * (1.0 - alpha) + v * alpha);
}

/**
 * Cartesion to polar coordinates conversion.
 * Result:
 * [0] = length
 * [1] = angle with z-axis
 * [2] = angle of projection into x,y, plane with x-axis
 */
template <class T>
   inline Vec3<T> cartesianToPolar (const Vec3<T> & v) {
    Vec3<T> polar;
    polar[0] = length (v);

    if (v[2] > 0.0f) {
      polar[1] = (T) atan (sqrt (v[0] * v[0] + v[1] * v[1]) / v[2]);
    }
    else if (v[2] < 0.0f) {
      polar[1] = (T) atan (sqrt (v[0] * v[0] + v[1] * v[1]) / v[2]) + M_PI;
    }
    else {
      polar[1] = M_PI * 0.5f;
    }


    if (v[0] > 0.0f) {
      polar[2] = (T) atan (v[1] / v[0]);
    }
    else if (v[0] < 0.0f) {
      polar[2] = (T) atan (v[1] / v[0]) + M_PI;
    }
    else if (v[1] > 0) {
      polar[2] = M_PI * 0.5f;
    }
    else {
      polar[2] = -M_PI * 0.5;
    }
    return polar;
  }

/**
 * Polar to cartesion coordinates
 * input:
 * [0] = length
 * [1] = angle with z-axis
 * [2] = angle of projection into x,y, plane with x-axis
 */
template <class T>
inline Vec3<T> polarToCartesian (const Vec3<T> & v) {
  Vec3<T> cart;
  cart[0] = v[0] * (T) sin (v[1]) * (T) cos (v[2]);
  cart[1] = v[0] * (T) sin (v[1]) * (T) sin (v[2]);
  cart[2] = v[0] * (T) cos (v[1]);
  return cart;
}

template <class T>
inline Vec3<T> projectOntoVector (const Vec3<T> & v1, const Vec3<T> & v2) {
  return v2 * dotProduct (v1, v2);
}

template <class T>
inline Vec3<T> operator * (const T &s, const Vec3<T> &P) {
  return (P * s);
}

template <class T>
std::ostream & operator<< (std::ostream & output, const Vec3<T> & v) {
  output << v[0] << " " << v[1] << " " << v[2];
  return output;
}

template <class T>
std::istream & operator>> (std::istream & input, Vec3<T> & v) {
  input >> v[0] >> v[1] >> v[2];
  return input;
}

typedef Vec3<float> Vec3f;
typedef Vec3<double> Vec3d;
typedef Vec3<int> Vec3i;
