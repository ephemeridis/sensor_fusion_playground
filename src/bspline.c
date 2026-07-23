#include "bspline.h"

#include <stdlib.h>

bool bspline_init(bspline_t *spline, int p, int n, int k) {
  spline->p = p;
  spline->n = n;
  spline->k = k;

  spline->m = n + p + k + 1;
  spline->U = malloc(sizeof(double) * spline->m);

  return true;
}

void bspline_free(bspline_t *spline) {
  if (spline) {
    if (spline->U) free(spline->U);
  }
}

int bspline_find_span(double u, bspline_t *spline) {
  if (spline == NULL) return -1;

  int last = spline->m - spline->p - 2;
  if (u >= spline->U[last+1]) return last;        // end of domain
  if (u <= spline->U[spline->p]) return spline->p;

  // do binary search
  // with U[low] <= u < U[high]
  int low = spline->p;
  int high = last+1;
  while (low+1 < high) {
    int mid = (low + high)/2;
    if (u < spline->U[mid]) high = mid;
    else low = mid;
  }

  return low;
}

// algorithm A2.2 from the NURBS book
void bspline_basis_funcs(int i, double u, bspline_t *spline, double *N, size_t len_N) {
  if (spline == NULL) return;
  if (spline->U == NULL || N == NULL || len_N < spline->p+1) return;

  // store values to reduce repeat computation
  double left[spline->p], right[spline->p];

  // initial zero-th order basis function
  N[0] = 1.0;
  for (int j = 0; j < spline->p; j++) {
    // each loop computes one more degree of the basis function
    left[j] = u - spline->U[i-j];
    right[j] = spline->U[i+j+1] - u;
    double saved = 0.0;
    for (int r = 0; r <= j; r++) {
      double tmp = N[r] / (right[r] + left[j-r]);
      N[r] = saved + right[r]*tmp;
      saved = left[j-r]*tmp;
    }
    N[j+1] = saved;
  }
}

// algorithm A2.3 from the NURBS book
void bspline_basis_derivs(int i, double u, int N_derivs, bspline_t *spline, double *dN_u, size_t len_dN_u, size_t stride_dN_u) {
  if (spline == NULL) return;
  if (spline->U == NULL || dN_u == NULL || len_dN_u < N_derivs+1 || stride_dN_u < spline->p+1) return;

  // utility to more easily access the output matrix
  double *ders[N_derivs+1];
  for (int k = 0; k <= N_derivs; k++) {
    ders[k] = spline->U + k*stride_dN_u;
  }

}