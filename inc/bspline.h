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
void bspline_basis_funcs(int i, double u, bspline_t *spline, double *N, size_t len_N);
void bspline_basis_derivs(int i, double u, int N_derivs, bspline_t *spline, double *dN_u, size_t len_dN_u, size_t stride_dN_u);

#endif //SENSOR_FUSION_PLAYGROUND_SPLINE_H
