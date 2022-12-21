#pragma once

#include <functional>
#include <geometry.hh>

namespace Nielson {

  using Boundary = std::function<Geometry::Point3D(double)>;
  using NormalFence = std::function<Geometry::Vector3D(double)>;

  Geometry::TriMesh
  evaluate(const std::array<Boundary, 3> &boundaries,
           const std::array<NormalFence, 3> &fences,
           double fullness,
           size_t resolution);

  Geometry::TriMesh
  evaluate(const Geometry::Point3D &p1,
           const Geometry::Point3D &p2,
           const Geometry::Point3D &p3,
           const Geometry::Vector3D &n1,
           const Geometry::Vector3D &n2,
           const Geometry::Vector3D &n3,
           double fullness,
           size_t resolution);

  Geometry::TriMesh
  evaluate(const Geometry::TriMesh &mesh,
           const Geometry::PointVector &normals,
           double fullness,
           size_t resolution);

}
