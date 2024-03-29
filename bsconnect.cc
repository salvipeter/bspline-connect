#include "bsconnect.hh"

#include "Eigen/LU"

using namespace Geometry;

static Vector3D surfaceNormal(const BSSurface &s, double u, double v) {
  VectorMatrix der;
  s.eval(u, v, 1, der);
  return (der[1][0] ^ der[0][1]).normalize();
}

// `p` (supposedly in the `(u,v)` plane) is written as a combination of `u` and `v`.
// Returns `{ a, b }` s.t. `u * a + v * b = p`.
static std::pair<double, double> inSystem(const Vector3D &u, const Vector3D &v, const Vector3D &p) {
  double u2 = u * u, v2 = v * v, pu = p * u, pv = p * v, uv = u * v;
  double denom = u2 * v2 - uv * uv;
  return { (v2 * pu - uv * pv) / denom, (u2 * pv - uv * pu) / denom };
}

static double normalCurvature(const BSSurface &s, double u, double v, const Vector3D &dir) {
  VectorMatrix der;
  s.eval(u, v, 2, der);
  auto n = (der[1][0] ^ der[0][1]).normalize();
  double E = der[1][0] * der[1][0];
  double F = der[1][0] * der[0][1];
  double G = der[0][1] * der[0][1];
  double L = n * der[2][0];
  double M = n * der[1][1];
  double N = n * der[0][2];
  auto [du, dv] = inSystem(der[1][0], der[0][1], dir);
  double du2 = du * du, dudv = du * dv, dv2 = dv * dv;
  return (L * du2 + 2 * M * dudv + N * dv2) / (E * du2 + 2 * F * dudv + G * dv2);
}

static void connectC0(const BSSurface &master, BSSurface &slave, size_t fixed) {
  size_t m = slave.numControlPoints()[1];
  for (size_t j = fixed; j < m - fixed; ++j)
    slave.controlPoint(0, j) = master.controlPoint(0, j);
}

static void connectG1(const BSSurface &master, BSSurface &slave, size_t fixed, size_t resolution) {
  size_t m = slave.numControlPoints()[1];
  size_t p = slave.basisV().degree();

  DoubleMatrix der_u;
  slave.basisU().basisFunctionDerivatives(slave.basisU().findSpan(0.0), 0.0, 1, der_u);
  double dbase = der_u[1][1];

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(resolution, m - fixed * 2);
  Eigen::MatrixXd b(resolution, 3);
  using VecMap = Eigen::Map<const Eigen::Vector3d>;

  for (size_t k = 0; k < resolution; ++k) {
    double v = (double)k / (resolution - 1);
    auto n = surfaceNormal(master, 0, v);

    DoubleVector coeff_v;
    size_t span = slave.basisV().findSpan(v);
    slave.basisV().basisFunctions(span, v, coeff_v);
    for (size_t j = 0; j <= p; ++j)
      if (span - p + j >= fixed && span - p + j < m - fixed)
        A(k, span - p + j - fixed) = dbase * coeff_v[j];

    VectorMatrix der;
    slave.eval(0, v, 1, der);
    b.block<1,3>(k, 0) = VecMap((n * -(der[1][0] * n)).data());
  }

  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  for (size_t j = 0; j < m - fixed * 2; ++j)
    slave.controlPoint(1, fixed + j) += { x(j, 0), x(j, 1), x(j, 2) };
}

static void connectG2(const BSSurface &master, BSSurface &slave, size_t fixed, size_t resolution) {
  size_t m = slave.numControlPoints()[1];
  size_t p = slave.basisV().degree();
  size_t q = slave.basisU().degree();
  const auto &knots = slave.basisU().knots();
  double base = (knots[q+1] - knots[2]) * (knots[q+2] - knots[2]) / (q * (q - 1));

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(resolution, m - fixed * 2);
  Eigen::MatrixXd b(resolution, 3);
  using VecMap = Eigen::Map<const Eigen::Vector3d>;

  for (size_t k = 0; k < resolution; ++k) {
    double v = (double)k / (resolution - 1);
    VectorMatrix der;
    slave.eval(0, v, 1, der);
    auto d1 = der[1][0];

    double k_master = -normalCurvature(master, 0, v, d1);
    double k_slave  =  normalCurvature(slave,  0, v, d1);
    auto n = surfaceNormal(slave, 0, v);

    DoubleVector coeff_v;
    size_t span = slave.basisV().findSpan(v);
    slave.basisV().basisFunctions(span, v, coeff_v);
    for (size_t j = 0; j <= p; ++j)
      if (span - p + j >= fixed && span - p + j < m - fixed)
        A(k, span - p + j - fixed) = coeff_v[j];

    b.block<1,3>(k, 0) = VecMap((n * (k_master - k_slave) * d1.normSqr() * base).data());
  }

  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  for (size_t j = 0; j < m - fixed * 2; ++j)
    slave.controlPoint(2, fixed + j) += { x(j, 0), x(j, 1), x(j, 2) };
}

void connectBSplineSurfaces(const BSSurface &master, BSSurface &slave,
                            bool fix_c0, bool fix_g1, bool fix_g2,
                            const std::array<size_t,3> &fixed, size_t resolution) {
  if (fix_c0)
    connectC0(master, slave, fixed[0]);
  if (fix_g1)
    connectG1(master, slave, fixed[1], resolution);
  if (fix_g2)
    connectG2(master, slave, fixed[2], resolution);
}
