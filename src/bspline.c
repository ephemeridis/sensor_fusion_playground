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

// utility macro to more easily access the output matrix
#define INDEX_BUF(row,col,stride) ((row) * (stride) + (col))

// algorithm A2.2 from the NURBS book
// modified to store the knot differences in the lower triangle and the
// basis functions in the upper triangle, to be used in A2.3 for computing derivatives
void bspline_basis_funcs_store_ndu(
  int i, double u, bspline_t *spline, double *ndu_out, size_t rows_ndu_out, size_t stride_ndu_out) {
  if (spline == NULL) return;
  if (spline->U == NULL || ndu_out == NULL) return;
  if (rows_ndu_out < (size_t)(spline->p+1) || stride_ndu_out != (size_t)(spline->p+1)) return;

  const size_t s = stride_ndu_out;

  // store values to reduce repeat computation
  double left[spline->p], right[spline->p];

  // initial zero-th order basis function
  ndu_out[INDEX_BUF(0, 0, s)] = 1.0;
  for (int j = 0; j < spline->p; j++) {
    // each loop computes one more degree of the basis function
    left[j] = u - spline->U[i-j];
    right[j] = spline->U[i+j+1] - u;
    double saved = 0.0;
    for (int r = 0; r <= j; r++) {
      // compute knot differences in lower triangle
      // row j+1 means "a difference spanning j+1 knots"
      ndu_out[INDEX_BUF(j+1, r, s)] = right[r] + left[j-r];
      double tmp = ndu_out[INDEX_BUF(r, j, s)] / ndu_out[INDEX_BUF(j+1, r, s)];

      // compute basis functions in upper triangle
      // note: reads column j (degree j), writes column j+1
      // unlike A2.2 the lower-degree values must survive
      // so this is NOT done in place
      ndu_out[INDEX_BUF(r, j+1, s)] = saved + right[r]*tmp;
      saved = left[j-r]*tmp;
    }
    ndu_out[INDEX_BUF(j+1, j+1, s)] = saved;
  }
}

// algorithm A2.3 from the NURBS book
// use the _use_ndu form directly if you already have an ndu.
void bspline_basis_derivs(int i, double u, int N_derivs, bspline_t *spline, double *deriv_out, size_t rows_deriv_out, size_t stride_deriv_out) {
  if (spline == NULL) return;
  if (spline->U == NULL || deriv_out == NULL) return;
  if (spline->p < 1 || N_derivs < 0) return;
  if (rows_deriv_out < (size_t)(N_derivs+1) || stride_deriv_out != (size_t)(spline->p+1)) return;

  const size_t stride_ndu = (size_t)(spline->p + 1);
  double ndu[stride_ndu * stride_ndu];

  bspline_basis_funcs_store_ndu(i, u, spline, ndu, stride_ndu, stride_ndu);
  bspline_basis_derivs_use_ndu(i, u, N_derivs, spline,
                               deriv_out, rows_deriv_out, stride_deriv_out,
                               ndu, stride_ndu, stride_ndu);
}

void bspline_basis_derivs_use_ndu(
  int i, double u, int N_derivs, bspline_t *spline,
  double *deriv_out, size_t rows_deriv_out, size_t stride_deriv_out,
  double *ndu_in, size_t rows_ndu_in, size_t stride_ndu_in) {
  (void)i; (void)u;  // everything needed is already baked into ndu_in
  if (spline == NULL) return;
  if (spline->U == NULL || deriv_out == NULL || ndu_in == NULL) return;
  if (rows_deriv_out < (size_t)(N_derivs+1) || stride_deriv_out != (size_t)(spline->p+1)) return;
  if (rows_ndu_in < (size_t)(spline->p+1) || stride_ndu_in != (size_t)(spline->p+1)) return;

  // we assume that ndu_in is already populated with the basis functions and the knot differences
  // from calling bspline_basis_funcs_store_ndu

  const int p = spline->p;
  const size_t sd = stride_deriv_out, sn = stride_ndu_in;

  // load the basis functions
  // the last ndu column holds the degree-p functions, i.e. the 0th derivative
  for (int j = 0; j <= p; j++) {
    deriv_out[INDEX_BUF(0, j, sd)] = ndu_in[INDEX_BUF(j, p, sn)];
  }

  // a degree-p function is a piecewise polynomial of degree p, so anything
  // past the pth derivative is identically zero
  for (int k = p+1; k <= N_derivs; k++) {
    for (int j = 0; j <= p; j++) deriv_out[INDEX_BUF(k, j, sd)] = 0.0;
  }
  const int kmax = N_derivs < p ? N_derivs : p;

  // compute the derivatives
  // from Eq. 2.9: the kth derivative is a weighted sum of k+1 basis funcs of degree p-k
  // a[][] holds those weights, and only two rows are needed, alternating via s1/s2
  double a[2][p+1];
  for (int r = 0; r <= p; r++) {
    int s1 = 0;   // row a_{k-1,j}, read past layer of weights
    int s2 = 1;   // row a_{k,j},   write current layer of weights
    // these two get swapped back and forth
    a[0][0] = 1.0;

    for (int k = 1; k <= kmax; k++) {
      double d = 0.0;            // accumulates the sum for this derivative
      const int rk = r - k;      // row of the degree p-k function we need
      const int pk = p - k;      // reduced degree == its ndu column

      // j = 0 term; skipped when r < k because that basis func is zero
      if (r >= k) {
        a[s2][0] = a[s1][0] / ndu_in[INDEX_BUF(pk+1, rk, sn)];
        d = a[s2][0] * ndu_in[INDEX_BUF(rk, pk, sn)];
      }

      // clamp the middle band so rk+j stays inside the filled triangle
      const int j1 = (rk >= -1) ? 1 : -rk;
      const int j2 = (r-1 <= pk) ? k-1 : p-r;

      // middle terms j = 1 .. k-1
      for (int j = j1; j <= j2; j++) {
        a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu_in[INDEX_BUF(pk+1, rk+j, sn)];
        d += a[s2][j] * ndu_in[INDEX_BUF(rk+j, pk, sn)];
      }

      // j = k term; skipped when r > p-k, the far end of the same window
      if (r <= pk) {
        a[s2][k] = -a[s1][k-1] / ndu_in[INDEX_BUF(pk+1, r, sn)];
        d += a[s2][k] * ndu_in[INDEX_BUF(r, pk, sn)];
      }

      deriv_out[INDEX_BUF(k, r, sd)] = d;

      int t = s1; s1 = s2; s2 = t;   // switch rows
    }
  }

  // multiply through by the correct factors, p!/(p-k)!  (Eq. 2.9)
  // held out until now so it stays integer arithmetic
  int f = p;
  for (int k = 1; k <= kmax; k++) {
    for (int j = 0; j <= p; j++) deriv_out[INDEX_BUF(k, j, sd)] *= f;
    f *= (p - k);
  }
}