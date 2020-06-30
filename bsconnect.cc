#include "bsconnect.hh"

using namespace Geometry;

static void connectC0(const BSSurface &master, BSSurface &slave) {
  size_t m = master.numControlPoints()[1] - 1;
  for (size_t j = 0; j <= m; ++j)
    slave.controlPoint(0, j) = master.controlPoint(0, j);
}

static void connectG1(const BSSurface &master, BSSurface &slave, size_t resolution) {
  // TODO
}

static void connectG2(const BSSurface &master, BSSurface &slave, size_t resolution) {
  // TODO
}

void connectBSplineSurfaces(const BSSurface &master, BSSurface &slave,
                            bool fix_c0, bool fix_g1, bool fix_g2, size_t resolution) {
  if (fix_c0)
    connectC0(master, slave);
  if (fix_g1)
    connectG1(master, slave, resolution);
  if (fix_g2)
    connectG2(master, slave, resolution);
}
