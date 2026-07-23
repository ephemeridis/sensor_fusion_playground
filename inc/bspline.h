#ifndef SENSOR_FUSION_PLAYGROUND_SPLINE_H
#define SENSOR_FUSION_PLAYGROUND_SPLINE_H
#include <stdbool.h>
#include <stddef.h>

typedef struct {
  int p; // order of the bspline approximation
  int n; // number of points to interpolate between
  int k; // number of end derivatives to consider
  int m; // size of the knot vector
  double *U; // knot vector
} bspline_t;

bool bspline_init(bspline_t *spline, int p, int n, int k);
void bspline_free(bspline_t *spline);

int bspline_find_span(double u, bspline_t *spline);

void bspline_basis_funcs(
  int i, double u, bspline_t *spline, double *N, size_t len_N);
void bspline_basis_funcs_store_ndu(
  int i, double u, bspline_t *spline, double *ndu_out, size_t rows_ndu_out, size_t stride_ndu_out);

void bspline_basis_derivs(
  int i, double u, int N_derivs, bspline_t *spline,
  double *deriv_out, size_t rows_deriv_out, size_t stride_deriv_out);
void bspline_basis_derivs_use_ndu(
  int i, double u, int N_derivs, bspline_t *spline,
  double *deriv_out, size_t rows_deriv_out, size_t stride_deriv_out,
  double *ndu_in, size_t rows_ndu_in, size_t stride_ndu_in);

#endif //SENSOR_FUSION_PLAYGROUND_SPLINE_H
