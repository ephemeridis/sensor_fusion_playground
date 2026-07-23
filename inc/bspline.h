#ifndef SENSOR_FUSION_PLAYGROUND_SPLINE_H
#define SENSOR_FUSION_PLAYGROUND_SPLINE_H
#include <stdbool.h>
#include <stddef.h>
#include <math.h>
#include <float.h>

// scalar type selection
#ifdef BSPL_DOUBLE_PRECISION
typedef double bspl_scalar_t;
#define BSPL_FABS    fabs
#define BSPL_EPSILON DBL_EPSILON
#else
typedef float bspl_scalar_t;
#define BSPL_FABS    fabsf
#define BSPL_EPSILON FLT_EPSILON
#endif

// maximum spline degree supported. sizes the fixed scratch buffers used
// by the basis function and derivative algorithms
#ifndef BSPLINE_MAX_DEGREE
#define BSPLINE_MAX_DEGREE 6
#endif
#define BSPLINE_MAX_ORDER (BSPLINE_MAX_DEGREE + 1)

typedef struct {
  int p; // order of the bspline approximation
  int n; // number of points to interpolate between
  int k; // number of end derivatives to consider
  int m; // size of the knot vector
  bspl_scalar_t *U; // knot vector
  bool heap_alloc;  // only bspline_init sets this true
} bspline_t;

// required length of the knot vector for a given (p, n, k)
#define BSPLINE_KNOTS_LEN(P,N,K) ((N) + (P) + (K) + 1)

// declares a static knot buffer for use with bspline_init_static.
#define BSPLINE_STATIC_KNOTS(name, P, N, K) \
  static bspl_scalar_t name##_knots[BSPLINE_KNOTS_LEN((P), (N), (K))]

bool bspline_init(bspline_t *spline, int p, int n, int k);
bool bspline_init_static(bspline_t *spline, int p, int n, int k,
                         bspl_scalar_t *U_buf, size_t U_len);
void bspline_free(bspline_t *spline);

int bspline_find_span(bspl_scalar_t u, bspline_t *spline);

void bspline_basis_funcs(
  int i, bspl_scalar_t u, bspline_t *spline, bspl_scalar_t *N, size_t len_N);
void bspline_basis_funcs_store_ndu(
  int i, bspl_scalar_t u, bspline_t *spline,
  bspl_scalar_t *ndu_out, size_t rows_ndu_out, size_t stride_ndu_out);

void bspline_basis_derivs(
  int i, bspl_scalar_t u, int N_derivs, bspline_t *spline,
  bspl_scalar_t *deriv_out, size_t rows_deriv_out, size_t stride_deriv_out);
void bspline_basis_derivs_use_ndu(
  int i, bspl_scalar_t u, int N_derivs, bspline_t *spline,
  bspl_scalar_t *deriv_out, size_t rows_deriv_out, size_t stride_deriv_out,
  bspl_scalar_t *ndu_in, size_t rows_ndu_in, size_t stride_ndu_in);

#endif //SENSOR_FUSION_PLAYGROUND_SPLINE_H