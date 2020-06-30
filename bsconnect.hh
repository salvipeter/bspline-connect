#pragma once

#include <geometry.hh>

// Connecting B-spline surfaces based on the paper
// P. Salvi, T. Várady: Hierarchical surface fairing with constraints.
// Proceedings of the 14th ACM Symposium on Solid and Physical Modeling, pp. 195-200, 2010.

// `master' and `slave' are assumed to have the same degree and knot vector
// in the `u' parametric direction. The common boundary is `u=0',
// and the `v'-orientation is compatible, i.e., C0-continuity is `M(0, v) = S(0, v)'.
// `resolution' is used only in G1- and G2-fixing.
void connectBSplineSurfaces(const Geometry::BSSurface &master, Geometry::BSSurface &slave,
                            bool fix_c0, bool fix_g1, bool fix_g2, size_t resolution);
