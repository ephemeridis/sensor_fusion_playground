#ifndef SENSOR_FUSION_PLAYGROUND_TRAJ_H
#define SENSOR_FUSION_PLAYGROUND_TRAJ_H
#include <stdbool.h>
#include <stddef.h>

// traj marries the path and bspline modules, solved through matrix, so it
// owns all three includes and the precision agreement between them.
#include "bspline.h"
#include "matrix.h"
#include "path.h"

// Waypoint data arrives as path_scalar_t and evaluated output is returned in
// the same type, so traj has no scalar type of its own. Scratch buffers
// handed to bspline_* are bspl_scalar_t, and matrix entries are mat_scalar_t;
// traj.c converts at those boundaries.
//
// Mixing the three precisions compiles and runs correctly -- every buffer is
// typed for the module that writes it -- but it silently narrows at each
// boundary, which is rarely what anyone wants. Define
// TRAJ_ALLOW_MIXED_PRECISION if you mean it.
#ifndef TRAJ_ALLOW_MIXED_PRECISION
#if (defined(MAT_DOUBLE_PRECISION)  \
   + defined(BSPL_DOUBLE_PRECISION) \
   + defined(PATH_DOUBLE_PRECISION)) % 3 != 0
#error "precision flags disagree: set all of MAT_/BSPL_/PATH_DOUBLE_PRECISION or none (or define TRAJ_ALLOW_MIXED_PRECISION)"
#endif
#endif

// Degree, and number of extra control points bought by the end conditions.
// K_EXTRA = 4 == (velocity + acceleration) at each of the two ends.
#define TRAJ_P       5
#define TRAJ_K_EXTRA 4

// control points (== matrix dimension) and knot count for N data points
#define TRAJ_NCTRL(N)  ((N) + TRAJ_K_EXTRA)
#define TRAJ_NKNOTS(N) ((N) + TRAJ_P + TRAJ_K_EXTRA + 1)

// Boundary conditions for one scalar axis (x, y or theta).
typedef struct {
  path_scalar_t v0, a0;   // velocity / acceleration at t[0]
  path_scalar_t vn, an;   // velocity / acceleration at t[N-1]
} traj_bc_t;

// Clamped knot vector with the interior knots placed at the interior data
// times. Requires spline->m == TRAJ_NKNOTS(N) and N >= 2.
bool traj_build_knots(bspline_t *s, const path_scalar_t *t, int N);

// Fill A (nctrl x nctrl) with the collocation matrix. Depends only on the
// times, so assemble and factor this ONCE per segment.
bool traj_assemble(mat_t *A, bspline_t *s, const path_scalar_t *t, int N);

// Fill column `col` of a (nctrl x naxes) right-hand side with one axis.
// q is the N data values for that axis.
bool traj_rhs_col(mat_t *B, int col, const path_scalar_t *q, int N, const traj_bc_t *bc);

// Solve A * CTRL = B for every axis at once, given a prior
// mat_lup_decomp(A, L, U, Pm). CTRL and B are both nctrl x naxes and must
// each differ from Pm.
bool traj_solve(mat_t *CTRL, mat_t *B, mat_t *L, mat_t *U, mat_t *Pm);

// Evaluate derivatives 0..n_der of axis `col` at time u.
// out[k] receives the k-th time derivative.
bool traj_eval(path_scalar_t *out, int n_der, path_scalar_t u,
               bspline_t *s, mat_t *CTRL, int col);

// Convenience wrapper: everything needed to interpolate one segment.
bool traj_segment_solve(const segment_t *seg, bspline_t *s,
                        mat_t *A, mat_t *L, mat_t *U, mat_t *Pm,
                        mat_t *B, mat_t *CTRL,
                        const traj_bc_t *bcx, const traj_bc_t *bcy);

#endif //SENSOR_FUSION_PLAYGROUND_TRAJ_H
