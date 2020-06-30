#include "bsconnect.hh"

#include "Eigen/LU"

using namespace Geometry;

static void connectC0(const BSSurface &master, BSSurface &slave) {
  size_t m = master.numControlPoints()[1];
  for (size_t j = 0; j < m; ++j)
    slave.controlPoint(0, j) = master.controlPoint(0, j);
}

static void connectG1(const BSSurface &master, BSSurface &slave, size_t resolution) {
  size_t fixed = 1;             // # of fixed control points at both ends
  size_t m = master.numControlPoints()[1];
  size_t p = master.basisV().degree();

  DoubleMatrix der_u;
  slave.basisU().basisFunctionDerivatives(slave.basisU().findSpan(0.0), 0.0, 1, der_u);
  double dbase = der_u[1][1];

  Eigen::MatrixXd A(resolution, m - fixed * 2);
  Eigen::MatrixXd b(resolution, 3);
  using VecMap = Eigen::Map<const Eigen::Vector3d>;

  for (size_t k = 0; k < resolution; ++k) {
    double v = (double)k / (resolution - 1);
    VectorMatrix der;
    master.eval(0, v, 1, der);
    auto n = (der[1][0] ^ der[0][1]).normalize();

    slave.eval(0, v, 1, der);
    b.block<1,3>(k, 0) = VecMap((n * -(der[1][0] * n)).data());

    DoubleVector coeff_v;
    double span = slave.basisV().findSpan(v);
    slave.basisV().basisFunctions(span, v, coeff_v);
    for (size_t j = 0; j <= p; ++j)
      if (span - p + j >= fixed && span - p + j < m - fixed)
        A(k, span - p + j - fixed) = dbase * coeff_v[j];
  }

  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  for (size_t j = 0; j < m - fixed * 2; ++j)
    slave.controlPoint(1, fixed + j) += { x(j, 0), x(j, 1), x(j, 2) };
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
