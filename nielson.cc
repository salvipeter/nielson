#include "nielson.hh"

#include <cassert>
#include <cmath>

using namespace Geometry;

// Transfinite patch evaluation

static double
hermite(size_t type, double t) {
  switch (type) {
  case 0: return 1 - 3 * t * t + 2 * t * t * t;
  case 1: return 3 * t * t - 2 * t * t * t;
  case 2: return t - 2 * t * t + t * t * t; // <- Note typo in the paper
  case 3: return t * t * t - t * t;
  default: assert("invalid Hermite blend type" && false);
  }
}

static Point3D
hermiteCurve(const Point3D &p0, const Point3D &p1,
             const Vector3D &n0, const Vector3D &n1,
             double fullness, double t) {
  // This is not written in the paper, but the tangents obviously
  // need to be scaled, e.g. by the distance of the endpoints.
  double d = (p1 - p0).norm();
  return p0 * hermite(0, t) + p1 * hermite(1, t) +
    ((n0 ^ (p1 - p0)) ^ n0).normalize() * d * hermite(2, t) * fullness +
    ((n1 ^ (p1 - p0)) ^ n1).normalize() * d * hermite(3, t) * fullness;
}

static Point3D
transfinite(const std::array<Nielson::Boundary, 3> &boundaries,
            const std::array<Nielson::NormalFence, 3> &fences,
            double fullness, const Point3D &b) {
  Point3D result(0, 0, 0);
  double bsum = 0;
  size_t nonzero = 0;
  for (size_t i = 0; i < 3; ++i) {
    size_t j = (i + 1) % 3, k = (j + 1) % 3;
    auto p =
      hermiteCurve(boundaries[k](0), boundaries[i](b[k] / (1 - b[i])),
                   fences[k](0), fences[i](b[k] / (1 - b[i])),
                   fullness, 1 - b[i]);
    if (b[i] > epsilon)
      nonzero = i;
    double blend = b[j] * b[j] * b[k] * b[k];
    // double blend = b[j] * b[k]; // this can be used when ribbons interpolate all sides
    result += p * blend;
    bsum += blend;
  }
  return bsum > epsilon ? result / bsum : boundaries[(nonzero+1)%3](1);
}


// Boundary constraint generation

static Nielson::Boundary
boundary(const Point3D &p0, const Point3D &p1,
         const Vector3D &n0, const Vector3D &n1,
         double fullness) {
  return [=](double t) { return hermiteCurve(p0, p1, n0, n1, fullness, t); };
}

static Nielson::NormalFence
fence(const Point3D &p0, const Point3D &p1,
      const Vector3D &n0, const Vector3D &n1,
      double fullness) {
  auto b = boundary(p0, p1, n0, n1, fullness);
  double h = 1e-8;
  auto der = [=](double t) {
    if (t + h > 1)
      return (b(t) - b(t - h)) / h;
    return (b(t + h) - b(t)) / h;
  };
  auto d0 = der(0), d1 = der(1);
  return [=](double t) {
    return (der(t) ^ ((n0 ^ d0) * (1 - t) + (n1 ^ d1) * t)).normalize();
  };
}


// Parameterization & tessellation

static PointVector
parameters(size_t resolution) {
  PointVector params, vertices = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
  for (size_t j = 0; j <= resolution; ++j) {
    double u = (double)j / resolution;
    auto p = vertices[0] * u + vertices[2] * (1 - u);
    auto q = vertices[1] * u + vertices[2] * (1 - u);
    for (size_t k = 0; k <= j; ++k) {
      double v = j == 0 ? 1.0 : (double)k / j;
      params.push_back(p * (1 - v) + q * v);
    }
  }
  return params;
}

static TriMesh
triangles(size_t resolution) {
  TriMesh mesh;
  mesh.resizePoints((resolution + 1) * (resolution + 2) / 2);
  size_t prev = 0, current = 1;
  for (size_t i = 0; i < resolution; ++i) {
    for (size_t j = 0; j < i; ++j) {
      mesh.addTriangle(current + j, current + j + 1, prev + j);
      mesh.addTriangle(current + j + 1, prev + j + 1, prev + j);
    }
    mesh.addTriangle(current + i, current + i + 1, prev + i);
    prev = current;
    current += i + 2;
  }
  return mesh;
}


// Main functions

namespace Nielson {

  TriMesh
  evaluate(const std::array<Boundary, 3> &boundaries,
           const std::array<NormalFence, 3> &fences,
           double fullness, size_t resolution) {
    auto result = triangles(resolution);
    auto params = parameters(resolution);
    for (size_t i = 0; i < params.size(); ++i)
      result[i] = transfinite(boundaries, fences, fullness, params[i]);
    return result;
  }

  TriMesh
  evaluate(const Point3D &p1, const Point3D &p2, const Point3D &p3,
           const Vector3D &n1, const Vector3D &n2, const Vector3D &n3,
           double fullness, size_t resolution) {
    return evaluate({
        boundary(p1, p2, n1, n2, fullness),
        boundary(p2, p3, n2, n3, fullness),
        boundary(p3, p1, n3, n1, fullness)
      }, {
        fence(p1, p2, n1, n2, fullness),
        fence(p2, p3, n2, n3, fullness),
        fence(p3, p1, n3, n1, fullness)
      },
      fullness, resolution);
  }

  TriMesh
  evaluate(const TriMesh &mesh, const PointVector &normals,
           double fullness, size_t resolution) {
    TriMesh result;
    for (const auto &tri : mesh.triangles()) {
      auto surface =
        evaluate(mesh[tri[0]], mesh[tri[1]], mesh[tri[2]],
                 normals[tri[0]], normals[tri[1]], normals[tri[2]],
                 fullness, resolution);
      result.append(surface);
    }
    return result;
  }

}
