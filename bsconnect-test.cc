#include <fstream>

#include "bsconnect.hh"

using namespace Geometry;

BSSurface readBSS(std::string filename) {
  std::ifstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);

  size_t nu, nv, du, dv;
  double u, v, x, y, z;
  DoubleVector knots_u, knots_v;
  PointVector cpts;

  f >> du >> dv >> nu >> nv;

  size_t n_knots_u = nu + du + 1, n_knots_v = nv + dv + 1;
  for (size_t i = 0; i < n_knots_u; ++i) {
    f >> u;
    knots_u.push_back(u);
  }
  for (size_t j = 0; j < n_knots_v; ++j) {
    f >> v;
    knots_v.push_back(v);
  }

  for (size_t i = 0; i < nu; ++i)
    for (size_t j = 0; j < nv; ++j) {
      f >> x >> y >> z;
      cpts.emplace_back(x, y, z);
    }
  
  return { du, dv, knots_u, knots_v, cpts };
}

void writeBSS(const BSSurface &s, std::string filename) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);

  f << s.basisU().degree() << ' ' << s.basisV().degree() << std::endl;
  auto num_cp = s.numControlPoints();
  f << num_cp[0] << ' ' << num_cp[1] << std::endl;

  for (double u : s.basisU().knots())
    f << u << ' ';
  f << std::endl;
  for (double v : s.basisV().knots())
    f << v << ' ';
  f << std::endl;

  for (const auto &p : s.controlPoints())
    f << p << std::endl;
}

template<typename T>
void writeType(std::ostream &os, T x) {
  os.write(reinterpret_cast<const char *>(&x), sizeof(T));
}

void writeVector(std::ostream &os, const Vector3D &v) {
  writeType<float>(os, v[0]);
  writeType<float>(os, v[1]);
  writeType<float>(os, v[2]);
}

void writeSTL(const std::vector<BSSurface> &surfaces, std::string filename,
                     size_t resolution = 50) {
  size_t n = surfaces.size();
  std::ofstream f(filename, std::ios::binary);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  std::string comment("Sampled B-spline surfaces.");
  comment.resize(80, ' ');
  f.write(comment.c_str(), 80);
  writeType<uint32_t>(f, n * resolution * resolution * 2); // # of all faces
  PointVector points;
  auto addTriangle = [&](size_t i0, size_t i1, size_t i2) {
                       const auto &p0 = points[i0], &p1 = points[i1], &p2 = points[i2];
                       auto normal = ((p1 - p0) ^ (p2 - p0)).normalize();
                       writeVector(f, normal);
                       writeVector(f, p0);
                       writeVector(f, p1);
                       writeVector(f, p2);
                       writeType<uint16_t>(f, 0);
                     };
  for (const auto &s : surfaces) {
    points.clear();
    for (size_t i = 0; i <= resolution; ++i) {
      double u = (double)i / resolution;
      u = s.basisU().low() * (1 - u) + s.basisU().high() * u;
      for (size_t j = 0; j <= resolution; ++j) {
        double v = (double)j / resolution;
        v = s.basisV().low() * (1 - v) + s.basisV().high() * v;
        points.push_back(s.eval(u, v));
      }
    }
    for (size_t i = 0; i < resolution; ++i) {
      size_t base = i * (resolution + 1);
      for (size_t j = 0; j < resolution; ++j) {
        addTriangle(base + j, base + j + 1, base + j + resolution + 1);
        addTriangle(base + j + resolution + 1, base + j + 1, base + j + resolution + 2);
      }
    }
  }
}

void writeControlNet(const std::vector<BSSurface> &surfaces, std::string filename) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  std::vector<size_t> start_indices;
  size_t index = 1;
  for (const auto &s : surfaces) {
    start_indices.push_back(index);
    for (const auto &p : s.controlPoints())
      f << "v " << p << std::endl;
    index += s.controlPoints().size();
  }
  for (size_t i = 0; i < surfaces.size(); ++i) {
    const auto &s = surfaces[i];
    size_t base = start_indices[i];
    auto [nu, nv] = s.numControlPoints();
    for (size_t j = 0; j < nu; ++j)
      for (size_t k = 1; k < nv; ++k)
        f << "l " << base + j * nv + k - 1 << ' ' << base + j * nv + k << std::endl;
    for (size_t j = 1; j < nu; ++j)
      for (size_t k = 0; k < nv; ++k)
        f << "l " << base + (j - 1) * nv + k << ' ' << base + j * nv + k << std::endl;
  }
}

int main(int argc, char **argv) {
  if (argc < 3 || argc > 4) {
    std::cerr << "Usage: " << argv[0] << " <master.bss> <slave.bss> [resolution]" << std::endl;
    return 1;
  }

  size_t resolution = 100;
  if (argc == 4)
    resolution = std::atoi(argv[3]);

  auto master = readBSS(argv[1]);
  auto slave = readBSS(argv[2]);

  connectBSplineSurfaces(master, slave, true, true, true, { 0, 0, 0 }, resolution);

  writeBSS(slave, "output.bss");

  // For debugging
  slave.reverseV();
  writeSTL({ master, slave }, "output.stl");
  writeControlNet({ master, slave }, "output.obj");
}
