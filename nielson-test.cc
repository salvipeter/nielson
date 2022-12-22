#include <fstream>
#include <sstream>

#include "nielson.hh"

using namespace Geometry;

void
writeNormals(const std::vector<Point3D> &vertices,
             const std::vector<Vector3D> &normals,
             std::string filename) {
  std::ofstream f(filename);
  f << "# vtk DataFile Version 2.0" << std::endl;
  f << "Vertices with principal curvature values & directions" << std::endl;
  f << "ASCII" << std::endl;
  f << "DATASET POLYDATA" << std::endl;
  f << "POINTS " << vertices.size() << " float" << std::endl;
  for (const auto &v : vertices)
    f << v << std::endl;
  f << "POINT_DATA " << vertices.size() << std::endl;
  f << "NORMALS normal float" << std::endl;
  for (const auto &n : normals)
    f << n << std::endl;
}

// Weights according to:
//   N. Max, Weights for computing vertex normals from facet normals.
//     Journal of Graphics Tools, Vol. 4(2), 1999.
VectorVector
approximateNormals(const TriMesh &mesh) {
  VectorVector normals(mesh.points().size(), Vector3D(0, 0, 0));
  for (const auto &tri : mesh.triangles()) {
    for (size_t i = 0; i < 3; ++i) {
      size_t i0 = tri[i], i1 = tri[(i+1)%3], i2 = tri[(i+2)%3];
      auto v1 = mesh[i0] - mesh[i2], v2 = mesh[i1] - mesh[i0];
      auto w = v1.normSqr() * v2.normSqr();
      normals[i0] += (v1 ^ v2) / (w == 0.0 ? 1.0 : w);
    }
  }
  for (auto &n : normals)
    if (n.norm() > epsilon)
      n.normalize();
  return normals;
}

// Assumes that normal vectors have the same indices
// as the corresponding vertices.
VectorVector extractNormals(std::string filename) {
  VectorVector normals;
  std::ifstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  std::string line;
  std::istringstream ss;
  Vector3D n;
  while (!f.eof()) {
    std::getline(f, line);
    f >> std::ws;
    if (line.empty() || line[0] != 'v' || line[1] != 'n')
      continue;
    ss.str(line);
    ss.seekg(3); // skip the first three characters
    ss >> n[0] >> n[1] >> n[2];
    normals.push_back(n);
  }
  return normals;
}

int
main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    std::cerr << "Usage: " << argv[0]
              << " <input.obj> <output.stl> [fullness] [resolution]" << std::endl;
    std::cerr << "Defaults to fullness = 1 (should be positive) and resolution = 30." << std::endl;
    return 1;
  }
  double fullness = 1.0;
  if (argc >= 4)
    fullness = std::strtod(argv[3], nullptr);
  size_t resolution = 30;
  if (argc == 5)
    resolution = std::atoi(argv[4]);
  auto mesh = TriMesh::readOBJ(argv[1]);
  auto normals = extractNormals(argv[1]);
  if (normals.size() != mesh.points().size())
    normals = approximateNormals(mesh);
  // writeNormals(mesh.points(), normals, "/tmp/normals.vtk");
  Nielson::evaluate(mesh, normals, fullness, resolution).writeSTL(argv[2]);
}
