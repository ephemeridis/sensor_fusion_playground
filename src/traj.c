#include "traj.h"

#include <stddef.h>

#define DERS_ROWS 3   // value, 1st derivative, 2nd derivative

// build the knot vector
// for p = 5 + 4 extra control points for end derivatives
// there are exactly N-2 interior knots and N-2 interior data points
// so the knots go AT the interior data times
//
//   U[0 .. p]        = t[0]           (p+1 copies, clamped to t[0])
//   U[p+1 .. p+N-2]  = t[1 .. N-2]    (interior knots)
//   U[m-p-1 .. m-1]  = t[N-1]         (p+1 copies, clamped to t[N-1])
bool traj_build_knots(bspline_t *s, const path_scalar_t *t, int N) {
  if (s == NULL || s->U == NULL || t == NULL) return false;
  if (N < 2) return false;
  const int p = s->p;
  if (s->m != TRAJ_NKNOTS(N)) return false;

  for (int i = 0; i <= p; i++) s->U[i] = t[0];
  for (int j = 0; j <= N - 3; j++) s->U[p + 1 + j] = t[1 + j];
  for (int i = 0; i <= p; i++) s->U[s->m - 1 - i] = t[N - 1];

  // times must be monotonic increasing or the system is singular
  for (int i = p; i < s->m - p - 1; i++) {
    if (!(s->U[i] < s->U[i + 1])) return false;
  }
  return true;
}

// Collocation matrix
// to set up the system of linear equations to solve
// row 0      : C  (t[0])   = q[0]        reduces to P_0 = q[0], starting pos constraint
// row 1      : C' (t[0])   = v0          first derivative constraint
// row 2      : C''(t[0])   = a0          second derivative constraint
// rows 3..N  : C  (t[k])   = q[k],  k = 1 .. N-2
// row  nc-3  : C''(t[N-1]) = an
// row  nc-2  : C' (t[N-1]) = vn
// row  nc-1  : C  (t[N-1]) = q[N-1]      reduces to P_{nc-1} = q[N-1]
//
// Rows 0 and nc-1 are written from the value row of DersBasisFuns rather than
// hard-coded to 1 since at a clamped end that row is [1,0,0,0,0,0]
bool traj_assemble(mat_t *A, bspline_t *s, const path_scalar_t *t, int N) {
  if (A == NULL || s == NULL || t == NULL) return false;
  if (N < 2) return false;
  const int p  = s->p;
  const int nc = TRAJ_NCTRL(N);
  if (A->rows != nc || A->cols != nc) return false;
  if (!mat_fill(A, (mat_scalar_t)0)) return false;

  // NOTE: bspline_basis_derivs requires stride == p+1 EXACTLY, so this must be
  // a flat buffer: [DERS_ROWS][BSPLINE_MAX_ORDER] has stride 7, not 6, and the
  // call would silently do nothing. The element type must be bspl_scalar_t --
  // bspline.c writes through this pointer with its own scalar type.
  bspl_scalar_t ders[DERS_ROWS * BSPLINE_MAX_ORDER];
  bspl_scalar_t Nb[BSPLINE_MAX_ORDER];
  const size_t sd = (size_t)(p + 1);

  // start conditions, first 3 rows
  int span = bspline_find_span(t[0], s);
  if (span < 0) return false;
  bspline_basis_derivs(span, t[0], 2, s, ders, DERS_ROWS, sd);
  for (int j = 0; j <= p; j++) {
    const int col = span - p + j;
    A->mat[0][col] = (mat_scalar_t)ders[0 * sd + j];
    A->mat[1][col] = (mat_scalar_t)ders[1 * sd + j];
    A->mat[2][col] = (mat_scalar_t)ders[2 * sd + j];
  }

  // interior data points
  for (int k = 1; k <= N - 2; k++) {
    span = bspline_find_span(t[k], s);
    if (span < 0) return false;
    bspline_basis_funcs(span, t[k], s, Nb, sd);
    for (int j = 0; j <= p; j++) {
      A->mat[k + 2][span - p + j] = (mat_scalar_t)Nb[j];
    }
  }

  // end conditions
  span = bspline_find_span(t[N - 1], s);
  if (span < 0) return false;
  bspline_basis_derivs(span, t[N - 1], 2, s, ders, DERS_ROWS, sd);
  for (int j = 0; j <= p; j++) {
    const int col = span - p + j;
    A->mat[nc - 3][col] = (mat_scalar_t)ders[2 * sd + j];
    A->mat[nc - 2][col] = (mat_scalar_t)ders[1 * sd + j];
    A->mat[nc - 1][col] = (mat_scalar_t)ders[0 * sd + j];
  }

  return true;
}

bool traj_rhs_col(mat_t *B, int col, const path_scalar_t *q, int N, const traj_bc_t *bc) {
  if (B == NULL || q == NULL || bc == NULL) return false;
  if (N < 2) return false;
  const int nc = TRAJ_NCTRL(N);
  if (B->rows != nc || col < 0 || col >= B->cols) return false;

  B->mat[0][col]      = (mat_scalar_t)q[0];
  B->mat[1][col]      = (mat_scalar_t)bc->v0;
  B->mat[2][col]      = (mat_scalar_t)bc->a0;
  for (int k = 1; k <= N - 2; k++) B->mat[k + 2][col] = (mat_scalar_t)q[k];
  B->mat[nc - 3][col] = (mat_scalar_t)bc->an;
  B->mat[nc - 2][col] = (mat_scalar_t)bc->vn;
  B->mat[nc - 1][col] = (mat_scalar_t)q[N - 1];
  return true;
}

// mat_lup_decomp gives P A = L U, so A X = B becomes L U X = P B
// mat_forward_backward handles all RHS columns in one pass and can allow
// x == y, so this is one multiply and one triangular solve in-place
bool traj_solve(mat_t *CTRL, mat_t *B, mat_t *L, mat_t *U, mat_t *Pm) {
  if (CTRL == NULL || B == NULL || L == NULL || U == NULL || Pm == NULL) return false;
  if (CTRL->rows != B->rows || CTRL->cols != B->cols) return false;
  // mat_mul_to rejects aliasing of dst with either operand, and zeroes dst
  // itself, so CTRL must not be Pm.
  if (!mat_mul_to(CTRL, Pm, B)) return false;
  return mat_forward_backward(CTRL, CTRL, L, U);
}

// C^(k)(u) = sum_j ders[k][j] * P_{span-p+j}
bool traj_eval(path_scalar_t *out, int n_der, path_scalar_t u,
               bspline_t *s, mat_t *CTRL, int col) {
  if (out == NULL || s == NULL || CTRL == NULL) return false;
  if (n_der < 0 || n_der > s->p) return false;
  if (col < 0 || col >= CTRL->cols) return false;
  const int p = s->p;
  const size_t sd = (size_t)(p + 1);

  bspl_scalar_t ders[BSPLINE_MAX_ORDER * BSPLINE_MAX_ORDER];
  const int span = bspline_find_span(u, s);
  if (span < 0) return false;
  bspline_basis_derivs(span, u, n_der, s, ders, (size_t)(n_der + 1), sd);

  for (int k = 0; k <= n_der; k++) {
    double acc = 0.0;
    for (int j = 0; j <= p; j++) {
      acc += (double)ders[(size_t)k * sd + j] * (double)CTRL->mat[span - p + j][col];
    }
    out[k] = (path_scalar_t)acc;
  }
  return true;
}

// wrapper to solve one segment
// knots -> assemble -> factor -> rhs -> solve
// the caller owns every matrix so they can be reused across
// segments of the same size; A is preserved (mat_lup_decomp copies into U).
bool traj_segment_solve(const segment_t *seg, bspline_t *s,
                        mat_t *A, mat_t *L, mat_t *U, mat_t *Pm,
                        mat_t *B, mat_t *CTRL,
                        const traj_bc_t *bcx, const traj_bc_t *bcy) {
  if (seg == NULL || seg->n_points < 2) return false;
  const int N = (int)seg->n_points;

  if (!traj_build_knots(s, seg->t, N)) return false;
  if (!traj_assemble(A, s, seg->t, N)) return false;
  if (!mat_lup_decomp(A, L, U, Pm)) return false;
  if (!traj_rhs_col(B, 0, seg->x, N, bcx)) return false;
  if (!traj_rhs_col(B, 1, seg->y, N, bcy)) return false;
  return traj_solve(CTRL, B, L, U, Pm);
}
