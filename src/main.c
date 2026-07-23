#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "traj.h"

#define TOL 1e-4

static int failures = 0;

static void check(const char *name, double got, double want, double tol) {
  double err = fabs(got - want);
  if (err > tol) {
    printf("    FAIL %-28s got % .6f want % .6f (err %.2e)\n", name, got, want, err);
    failures++;
  } else {
    printf("    ok   %-28s % .6f\n", name, got);
  }
}

// sum of N = 1, sum of N' = 0, sum of N'' = 0 at any u
static void test_basis_sums(bspline_t *s, path_scalar_t u) {
  const int p = s->p;
  const size_t sd = (size_t)(p + 1);
  bspl_scalar_t ders[3 * BSPLINE_MAX_ORDER];

  int span = bspline_find_span(u, s);
  if (span < 0) { printf("    FAIL find_span(%g)\n", (double)u); failures++; return; }
  bspline_basis_derivs(span, u, 2, s, ders, 3, sd);

  double s0 = 0, s1 = 0, s2 = 0;
  for (int j = 0; j <= p; j++) {
    s0 += ders[0 * sd + j];
    s1 += ders[1 * sd + j];
    s2 += ders[2 * sd + j];
  }
  printf("  basis sums at u = %g (span %d)\n", (double)u, span);
  check("sum N", s0, 1.0, TOL);
  check("sum N'", s1, 0.0, TOL);
  check("sum N''", s2, 0.0, TOL);
}

static void test_segment(const path_scalar_t *t, const path_scalar_t *x,
                         const path_scalar_t *y,
                         int N, const traj_bc_t *bcx, const traj_bc_t *bcy) {
  const int nc = TRAJ_NCTRL(N);
  printf("  N = %d, control points = %d, knots = %d\n", N, nc, TRAJ_NKNOTS(N));

  bspline_t s;
  if (!bspline_init(&s, TRAJ_P, N, TRAJ_K_EXTRA)) {
    printf("    FAIL bspline_init\n"); failures++; return;
  }
  if (!traj_build_knots(&s, t, N)) {
    printf("    FAIL traj_build_knots\n"); failures++;
    bspline_free(&s); return;
  }

  printf("  U = {");
  for (int i = 0; i < s.m; i++) printf("%s%g", i ? ", " : "", (double)s.U[i]);
  printf("}\n");

  mat_t *A  = mat_new(nc, nc, 0);
  mat_t *L  = mat_new(nc, nc, 0);
  mat_t *U  = mat_new(nc, nc, 0);
  mat_t *Pm = mat_new(nc, nc, 0);
  mat_t *B  = mat_new(nc, 2, 0);
  mat_t *C  = mat_new(nc, 2, 0);
  if (!A || !L || !U || !Pm || !B || !C) {
    printf("    FAIL matrix alloc\n"); failures++; goto cleanup;
  }

  if (!traj_assemble(A, &s, t, N)) {
    printf("    FAIL traj_assemble\n"); failures++; goto cleanup;
  }
  if (nc <= 12) mat_print_prec(A, 3);

  if (!mat_lup_decomp(A, L, U, Pm)) {
    printf("    FAIL mat_lup_decomp (singular?)\n"); failures++; goto cleanup;
  }
  if (!traj_rhs_col(B, 0, x, N, bcx) || !traj_rhs_col(B, 1, y, N, bcy)) {
    printf("    FAIL traj_rhs_col\n"); failures++; goto cleanup;
  }
  if (!traj_solve(C, B, L, U, Pm)) {
    printf("    FAIL traj_solve\n"); failures++; goto cleanup;
  }

  test_basis_sums(&s, t[0]);
  test_basis_sums(&s, 0.5 * (t[0] + t[N - 1]));
  test_basis_sums(&s, t[N - 1]);

  printf("  interpolation residuals\n");
  for (int k = 0; k < N; k++) {
    path_scalar_t ox[1], oy[1];
    char lbl[48];
    if (!traj_eval(ox, 0, t[k], &s, C, 0) || !traj_eval(oy, 0, t[k], &s, C, 1)) {
      printf("    FAIL traj_eval at t[%d]\n", k); failures++; continue;
    }
    snprintf(lbl, sizeof lbl, "x(t[%d])", k); check(lbl, ox[0], x[k], TOL);
    snprintf(lbl, sizeof lbl, "y(t[%d])", k); check(lbl, oy[0], y[k], TOL);
  }

  printf("  endpoint derivatives\n");
  {
    path_scalar_t ox[3], oy[3];
    traj_eval(ox, 2, t[0], &s, C, 0);
    traj_eval(oy, 2, t[0], &s, C, 1);
    check("x'(t0)",  ox[1], bcx->v0, TOL);
    check("x''(t0)", ox[2], bcx->a0, TOL);
    check("y'(t0)",  oy[1], bcy->v0, TOL);
    check("y''(t0)", oy[2], bcy->a0, TOL);

    traj_eval(ox, 2, t[N - 1], &s, C, 0);
    traj_eval(oy, 2, t[N - 1], &s, C, 1);
    check("x'(tn)",  ox[1], bcx->vn, TOL);
    check("x''(tn)", ox[2], bcx->an, TOL);
    check("y'(tn)",  oy[1], bcy->vn, TOL);
    check("y''(tn)", oy[2], bcy->an, TOL);
  }

  printf("  sampled trajectory\n");
  for (int i = 0; i <= 10; i++) {
    path_scalar_t u = t[0] + (t[N - 1] - t[0]) * (path_scalar_t)i / 10;
    path_scalar_t ox[3], oy[3];
    traj_eval(ox, 2, u, &s, C, 0);
    traj_eval(oy, 2, u, &s, C, 1);
    printf("    t %6.3f  p (% .4f, % .4f)  v (% .4f, % .4f)  a (% .4f, % .4f)\n",
           (double)u, (double)ox[0], (double)oy[0], (double)ox[1],
           (double)oy[1], (double)ox[2], (double)oy[2]);
  }

cleanup:
  mat_free(A); mat_free(L); mat_free(U);
  mat_free(Pm); mat_free(B); mat_free(C);
  bspline_free(&s);
}

// straight line, rest to rest: positions must interpolate, all four
// endpoint derivatives must come back zero
static void test_line(void) {
  printf("\n== straight line, rest to rest ==\n");
  path_scalar_t t[4] = {0.0, 1.0, 2.0, 3.0};
  path_scalar_t x[4] = {0.0, 1.0, 2.0, 3.0};
  path_scalar_t y[4] = {0.0, 2.0, 4.0, 6.0};
  traj_bc_t z = {0, 0, 0, 0};
  test_segment(t, x, y, 4, &z, &z);
}

// uneven times and a turn, with nonzero end velocity to exercise the
// derivative rows
static void test_turn(void) {
  printf("\n== uneven times, nonzero end conditions ==\n");
  path_scalar_t t[5] = {0.0, 0.7, 1.9, 2.4, 4.0};
  path_scalar_t x[5] = {0.0, 3.0, -1.0, -4.0, -4.0};
  path_scalar_t y[5] = {0.0, 4.0,  4.0,  0.0, -3.0};
  traj_bc_t bx = {1.5, 0.0, -0.5, 0.0};
  traj_bc_t by = {0.0, 0.0,  2.0, 0.0};
  test_segment(t, x, y, 5, &bx, &by);
}

// minimum case: two points, no interior knots, a single quintic Bezier
static void test_two_points(void) {
  printf("\n== two points ==\n");
  path_scalar_t t[2] = {0.0, 2.0};
  path_scalar_t x[2] = {0.0, 5.0};
  path_scalar_t y[2] = {0.0, 1.0};
  traj_bc_t z = {0, 0, 0, 0};
  test_segment(t, x, y, 2, &z, &z);
}

static void test_path_file(const char *fname) {
  printf("\n== path file: %s ==\n", fname);
  path_t path;
  path.n_segments = 0;
  path.segments = NULL;
  path_load_fromfile(&path, fname);
  if (path.segments == NULL || path.n_segments == 0) {
    printf("  no segments loaded, skipping\n");
    return;
  }
  path_print(&path);

  traj_bc_t z = {0, 0, 0, 0};
  for (size_t i = 0; i < path.n_segments; i++) {
    segment_t *sg = &path.segments[i];
    printf("\n  -- segment %zu --\n", i);
    if (sg->n_points < 2) { printf("  too few points, skipping\n"); continue; }
    test_segment(sg->t, sg->x, sg->y, (int)sg->n_points, &z, &z);
  }
  path_free(&path);
}

int main(int argc, char **argv) {
  printf("scalar sizes: mat %zu, bspl %zu, path %zu\n",
         sizeof(mat_scalar_t), sizeof(bspl_scalar_t), sizeof(path_scalar_t));

  test_line();
  test_turn();
  test_two_points();
  if (argc > 1) test_path_file(argv[1]);

  printf("\n%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}