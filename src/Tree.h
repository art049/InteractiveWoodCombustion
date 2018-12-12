#pragma once
#include <cmath>
#include <vector>
#include "Vec3.h"
#include "Triangle.h"
#include "TreeGraph.h"


/// A Tree mesh class, storing a list of vertices and a list of triangles indexed over it.
class Tree {
public:
    Tree (TreeGraph * treegraph);
    inline virtual ~Tree () {}

    inline std::vector<Vec3f> & positions () { return m_positions; }
    inline const std::vector<Vec3f> & positions () const { return m_positions; }
    inline  std::vector<Vec3f> & normals () { return m_normals; }
    inline const std::vector<Vec3f> & normals () const { return m_normals; }
    inline std::vector<Triangle> triangles () { return m_triangles; }
    inline const std::vector<Triangle> & triangles () const { return m_triangles; }

    /// Empty the positions, normals and triangles arrays.
    void clear ();
    /// Compute smooth per-vertex normals
    void recomputeNormals ();
    /// scale to the unit cube and center at original
    void centerAndScaleToUnit ();

private:
    std::vector<Vec3f> m_positions;
    std::vector<Vec3f> m_normals;
    std::vector<Triangle> m_triangles;
};
