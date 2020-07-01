#pragma once

#include <geometry.hh>

// Connecting B-spline surfaces based on the paper
// P. Salvi, T. Várady: Hierarchical surface fairing with constraints.
// Proceedings of the 14th ACM Symposium on Solid and Physical Modeling, pp. 195-200, 2010.

// Capable of creating exact C0, and approximative G1 / G2 continuity between B-spline patches.

// Here `master` and `slave` are assumed to
// - have the same degree and the same clamped knot vector in the `v` parametric direction
// - parameterized (in both directions) in the range [0,1]
// - share a (possibly C0-discontinuous) common boundary at `u=0`
// - the `v`-orientation is compatible, i.e., C0-continuity is `M(0, v) = S(0, v)`
// `fix_c0`, `fix_g1` and `fix_g2` select the operations.
// `fixed` contains the number of fixed control points at both ends of the boundary,
//   for the C0, G1 and G2 case (in this order).
// Note that `resolution` is used only in G1- and G2-fixing.
void connectBSplineSurfaces(const Geometry::BSSurface &master, Geometry::BSSurface &slave,
                            bool fix_c0, bool fix_g1, bool fix_g2,
                            const std::array<size_t,3> &fixed, size_t resolution);
