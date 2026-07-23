#include "matrix.h"

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SWAP_CHUNK_SIZE 128

bool mat_init_static(mat_t *m, int rows, int cols, mat_scalar_t *data, mat_scalar_t **rowptrs) {
  if (!m || !data || !rowptrs || rows <= 0 || cols <= 0) return false;
  m->heap_alloc = false;
  m->rows = rows; m->cols = cols;
  m->_mat_data = data; m->mat = rowptrs;
  for (int i = 0; i < rows; i++) rowptrs[i] = data + (size_t)i * cols;
  return true;
}

mat_t *mat_new(int rows, int cols, mat_scalar_t fill_value) {
  if (rows <= 0 || cols <= 0) return NULL;
  mat_t *mat = malloc(sizeof(mat_t));
  if (mat == NULL) return NULL;
  mat->heap_alloc = true;

  mat->rows = rows;
  mat->cols = cols;
  mat->mat = (mat_scalar_t **)malloc(sizeof(mat_scalar_t *) * rows);
  mat->_mat_data = (mat_scalar_t *)malloc(sizeof(mat_scalar_t) * rows * cols);
  if (mat->mat == NULL || mat->_mat_data == NULL) {
    mat_free(mat);
    return NULL;
  }

  for (int i = 0; i < rows; i++) {
    mat->mat[i] = mat->_mat_data + (size_t)i * cols;
    for (int j = 0; j < cols; j++) {
      mat->mat[i][j] = fill_value;
    }
  }

  return mat;
}

mat_t *mat_zero(int rows, int cols) {
  return mat_new(rows, cols, 0);
}

mat_t *mat_identity(int rows, int cols) {
  mat_t *mat = mat_new(rows, cols, 0);   // NULL if dims invalid or alloc fails
  if (mat == NULL) return NULL;
  int n = (rows < cols) ? rows : cols;
  for (int i = 0; i < n; i++) mat->mat[i][i] = 1;
  return mat;
}

mat_t *mat_copy(mat_t *mat) {
  if (mat == NULL) return NULL;
  mat_t *out = mat_new(mat->rows, mat->cols, 0);
  if (out == NULL) return NULL;
  if (!mat_copy_to(out, mat)) {
    mat_free(out);
    return NULL;
  }
  return out;
}

bool mat_copy_to(mat_t *dst, mat_t *src) {
  if (dst == NULL || src == NULL) return false;
  if (dst->rows != src->rows || dst->cols != src->cols) return false;
  if (dst->_mat_data != src->_mat_data) {
    memcpy(dst->_mat_data, src->_mat_data,
           sizeof(mat_scalar_t) * (size_t)dst->rows * dst->cols);
  }
  return true;
}

mat_t *mat_add(mat_t *mat1, mat_t *mat2) {
  if (mat1 == NULL || mat2 == NULL) return NULL;
  if (mat1->rows != mat2->rows) return NULL;
  if (mat1->cols != mat2->cols) return NULL;

  mat_t *out = mat_new(mat1->rows, mat1->cols, 0);
  if (out == NULL) return NULL;
  if (!mat_add_to(out, mat1, mat2)) {
    mat_free(out);
    return NULL;
  }
  return out;
}

bool mat_add_to(mat_t *dst, mat_t *mat1, mat_t *mat2) {
  if (dst == NULL || mat1 == NULL || mat2 == NULL) return false;
  if (mat1->rows != mat2->rows || mat1->cols != mat2->cols) return false;
  if (dst->rows != mat1->rows || dst->cols != mat1->cols) return false;

  for (int i = 0; i < mat1->rows * mat1->cols; i++) {
    dst->_mat_data[i] = mat1->_mat_data[i] + mat2->_mat_data[i];
  }
  return true;
}

mat_t *mat_sub(mat_t *mat1, mat_t *mat2) {
  if (mat1 == NULL || mat2 == NULL) return NULL;
  if (mat1->rows != mat2->rows) return NULL;
  if (mat1->cols != mat2->cols) return NULL;

  mat_t *out = mat_new(mat1->rows, mat1->cols, 0);
  if (out == NULL) return NULL;
  if (!mat_sub_to(out, mat1, mat2)) {
    mat_free(out);
    return NULL;
  }
  return out;
}

bool mat_sub_to(mat_t *dst, mat_t *mat1, mat_t *mat2) {
  if (dst == NULL || mat1 == NULL || mat2 == NULL) return false;
  if (mat1->rows != mat2->rows || mat1->cols != mat2->cols) return false;
  if (dst->rows != mat1->rows || dst->cols != mat1->cols) return false;

  for (int i = 0; i < mat1->rows * mat1->cols; i++) {
    dst->_mat_data[i] = mat1->_mat_data[i] - mat2->_mat_data[i];
  }
  return true;
}

mat_t *mat_mul(mat_t *mat1, mat_t *mat2) {
  if (mat1 == NULL || mat2 == NULL) return NULL;
  if (mat1->cols != mat2->rows) return NULL;

  mat_t *out = mat_new(mat1->rows, mat2->cols, 0);
  if (out == NULL) return NULL;

  if (!mat_mul_to(out, mat1, mat2)) {
    mat_free(out);
    return NULL;
  }
  return out;
}


bool mat_mul_to(mat_t *dst,mat_t *mat1, mat_t *mat2) {
  if (dst == NULL || mat1 == NULL || mat2 == NULL) return false;
  // check matmul condition
  if (mat1->cols != mat2->rows) return false;
  // check output matrix size
  if (dst->rows != mat1->rows || dst->cols != mat2->cols) return false;
  // ensure that the destination is not the same as mat1, mat2
  if (dst->_mat_data == mat1->_mat_data || dst->_mat_data == mat2->_mat_data) return false;

  // zero the output
  memset(dst->_mat_data, 0, sizeof(mat_scalar_t) * (size_t)dst->rows * dst->cols);
  for (int i = 0; i < dst->rows; i++) {
    for (int j = 0; j < dst->cols; j++) {
      for (int k = 0; k < mat1->cols; k++) {
        dst->mat[i][j] += mat1->mat[i][k] * mat2->mat[k][j];
      }
    }
  }

  return true;
}


mat_t *mat_mul_scalar(mat_t *mat, mat_scalar_t scalar) {
  if (mat == NULL) return NULL;

  mat_t *out = mat_new(mat->rows, mat->cols, 0);
  if (out == NULL) return NULL;
  if (!mat_mul_scalar_to(out, mat, scalar)) {
    mat_free(out);
    return NULL;
  }
  return out;
}


bool mat_mul_scalar_to(mat_t *dst, mat_t *mat, mat_scalar_t scalar) {
  if (dst == NULL || mat == NULL) return false;
  if (dst->rows != mat->rows || dst->cols != mat->cols) return false;

  for (int i = 0; i < mat->rows * mat->cols; i++) {
    dst->_mat_data[i] = mat->_mat_data[i] * scalar;
  }
  return true;
}

mat_t *mat_transpose(mat_t *mat) {
  if (mat == NULL) return NULL;
  mat_t *out;
  if (mat->rows != mat->cols) {
    int rows_t = mat->cols;
    int cols_t = mat->rows;
    out = mat_new(rows_t, cols_t, 0);
    if (out == NULL) return NULL;
    if (!mat_transpose_to(out, mat)) {
      mat_free(out);
      return NULL;
    }
  } else {
    out = mat_copy(mat);
    if (out == NULL) return NULL;
    if (!mat_transpose_inplace_square(out)) {
      mat_free(out);
      return NULL;
    }
  }
  return out;
}

bool mat_transpose_to(mat_t *dst, mat_t *src) {
  if (dst == NULL || src == NULL) return false;
  if (dst->rows != src->cols || dst->cols != src->rows) return false;
  // use inplace_square instead
  if (dst->_mat_data == src->_mat_data) return false;

  for (int i = 0; i < src->rows; i++) {
    for (int j = 0; j < src->cols; j++) {
      dst->mat[j][i] = src->mat[i][j];
    }
  }

  return true;
}

bool mat_transpose_inplace_square(mat_t *mat) {
  if (mat == NULL) return false;
  if (mat->rows != mat->cols) return false;

  mat_scalar_t tmp;
  for (int i = 1; i < mat->rows; i++) {
    for (int j = 0; j < i; j++) {
      tmp = mat->mat[i][j];
      mat->mat[i][j] = mat->mat[j][i];
      mat->mat[j][i] = tmp;
    }
  }

  return true;
}


bool mat_lu_decomp(mat_t *mat, mat_t *l, mat_t *u) {
  // lu decomp without pivoting
  if (mat == NULL || l == NULL || u == NULL) return false;
  // users to provide their own L and U matrices
  // but check the dimensions are accurate
  if (mat->rows != mat->cols) return false; // only square matrices
  int N = mat->rows;
  if (l->cols != N || l->rows != N) return false;
  if (u->cols != N || u->rows != N) return false;

  // clear the lower triangular matrix
  memset(l->_mat_data, 0, (size_t)N * N * sizeof(mat_scalar_t));
  // set the upper triangular matrix to the original matrix
  if (!mat_copy_to(u, mat)) return false;

  // perform LU decomp row-wise
  for (int n = 0; n < N; n++) {
    // check if diagonal element is non-zero
    if (MAT_FABS(u->mat[n][n]) <= MAT_EPSILON) return false;
    l->mat[n][n] = 1;

    // iterate through the remaining non-zero rows
    for (int i = n + 1; i < N; i++) {
      mat_scalar_t l_in = u->mat[i][n] / u->mat[n][n];
      l->mat[i][n] = l_in;
      for (int j = n; j < N; j++) {
        u->mat[i][j] -= l_in * u->mat[n][j];
      }

      // set to zero to avoid floating point precision imperfection
      u->mat[i][n] = 0;
    }
  }

  return true;
}

static void row_swap_lazy(mat_t *mat, int a, int b) {
  // lazy swap by swapping the data pointers
  mat_scalar_t *row = mat->mat[a];
  mat->mat[a] = mat->mat[b];
  mat->mat[b] = row;
}

static void row_swap_process(mat_t *mat) {
  // go through a matrix and resolve all the lazy pointer operations
  // so that _mat_data is a contiguous representation
  if (mat == NULL) return;

  mat_scalar_t *row = mat->_mat_data;
  for (int i = 0; i < mat->rows; i++) {
    while (row != mat->mat[i]) {
      // swap necessary
      // figure out the row offset
      size_t ind = (mat->mat[i] - mat->_mat_data) / (size_t)mat->cols;
      // swap the row with this one
      // if the matrix is large, swap in chunks
      mat_scalar_t row_chunk[SWAP_CHUNK_SIZE > mat->cols ? mat->cols : SWAP_CHUNK_SIZE];
      int left = mat->cols;
      int offset = 0;
      while (left > 0) {
        int cpy_amt = (SWAP_CHUNK_SIZE > left) ? left : SWAP_CHUNK_SIZE;
        memcpy(row_chunk, mat->mat[i]+offset, cpy_amt * sizeof(mat_scalar_t));
        memcpy(mat->mat[i]+offset, mat->mat[ind]+offset, cpy_amt * sizeof(mat_scalar_t));
        memcpy(mat->mat[ind]+offset, row_chunk, cpy_amt * sizeof(mat_scalar_t));
        left -= cpy_amt;
        offset += cpy_amt;
      }
      // swap the pointers back again
      // since the data in contiguous memory has been swapped
      mat_scalar_t *tmp = mat->mat[i];
      mat->mat[i] = mat->mat[ind];
      mat->mat[ind] = tmp;
    }
    row += mat->cols;
  }
}

bool mat_lup_decomp(mat_t *mat, mat_t *l, mat_t *u, mat_t *p) {
  // lu decomp with pivoting
  if (mat == NULL || l == NULL || u == NULL || p == NULL) return false;
  // users to provide their own L and U matrices
  // but check the dimensions are accurate
  if (mat->rows != mat->cols) return false; // only square matrices
  int N = mat->rows;
  if (l->cols != N || l->rows != N) return false;
  if (u->cols != N || u->rows != N) return false;
  if (p->cols != N || p->rows != N) return false;

  // clear the lower triangular matrix
  memset(l->_mat_data, 0, (size_t)N * N * sizeof(mat_scalar_t));
  // set the upper triangular matrix to the original matrix
  if (!mat_copy_to(u, mat)) return false;
  // set the pivot matrix to the identity matrix
  memset(p->_mat_data, 0, (size_t)N*N*sizeof(mat_scalar_t));
  for (int n = 0; n < N; n++) p->mat[n][n] = 1;

  // perform LU decomp row-wise
  for (int n = 0; n < N; n++) {
    // find the largest absolute to pivot to
    int row_ind = n;
    mat_scalar_t max_a = 0.0;
    for (int i = n; i < N; i++) {
      mat_scalar_t a = MAT_FABS(u->mat[i][n]);
      if (a > max_a) {
        max_a = a;
        row_ind = i;
      }
    }

    // check if we managed to find a non-zero row first element
    // if not then the matrix is singular
    if (MAT_FABS(max_a) <= MAT_EPSILON) return false;

    // perform row swap using the pointers
    row_swap_lazy(p, row_ind, n);
    row_swap_lazy(l, row_ind, n);
    row_swap_lazy(u, row_ind, n);

    // init l matrix to
    l->mat[n][n] = 1;

    // iterate through the remaining non-zero rows
    for (int i = n + 1; i < N; i++) {
      mat_scalar_t l_in = u->mat[i][n] / u->mat[n][n];
      l->mat[i][n] = l_in;
      for (int j = n; j < N; j++) {
        u->mat[i][j] -= l_in * u->mat[n][j];
      }

      // set to zero to avoid floating point precision imperfection
      u->mat[i][n] = 0;
    }
  }

  // commit the row swaps to preserve the contiguous memory structure
  row_swap_process(p);
  row_swap_process(l);
  row_swap_process(u);

  return true;
}

bool mat_forward_backward(mat_t *x, mat_t *y, mat_t *l, mat_t *u) {
  // y should equal to Pb or b in PAx = Pb = LUx or Ax = b = LUx
  if (x == NULL || y == NULL || l == NULL || u == NULL) return false;
  // set N to the number of linear equations
  int N = x->rows;

  // check that X and Y are column vectors with N rows
  if (x->rows != y->rows || x->cols != 1 || y->cols != 1) return false;

  // check that L and U are NxN square matrices
  if (l->rows != N || l->cols != N) return false;
  if (u->rows != N || u->cols != N) return false;

  // use x as the working for Ux in forward substitution
  // L_i,i * (UX)_i = Y_i - sum from j=0 to i-1 of L_i,j * (UX)_j
  for (int n = 0; n < N; n++) {
    x->mat[n][0] = y->mat[n][0];
    for (int j = 0; j < n; j++) {
      x->mat[n][0] -= l->mat[n][j] * x->mat[j][0];
    }
    x->mat[n][0] /= l->mat[n][n];
  }

  // now x contains UX
  // we begin backward substitution
  for (int n = N-1; n >= 0; n--) {
    for (int j = N-1; j > n; j--) {
      x->mat[n][0] -= u->mat[n][j] * x->mat[j][0];
    }
    x->mat[n][0] /= u->mat[n][n];
  }

  return true;
}

void mat_free(mat_t *mat) {
  if (mat && mat->heap_alloc) {
    if (mat->mat) free(mat->mat);
    if (mat->_mat_data) free(mat->_mat_data);
    free(mat);
  }
}

void mat_print(mat_t *mat) {
  mat_print_prec(mat, 3);
}

void mat_print_prec(mat_t *mat, int precision) {
  assert(mat != NULL);

  int max_widths[mat->cols];
  for (int j = 0; j < mat->cols; j++) max_widths[j] = 0;

  for (int i = 0; i < mat->rows; i++) {
    for (int j = 0; j < mat->cols; j++) {
      int width = snprintf(NULL, 0, "%.*lf", precision, mat->mat[i][j]);
      if (width > max_widths[j]) max_widths[j] = width;
    }
  }

  printf("MAT (R: %d, C: %d)\n", mat->rows, mat->cols);
  for (int i = 0; i < mat->rows; i++) {
    printf("  [ ");
    for (int j = 0; j < mat->cols; j++) {
      printf("%*.*lf%s", max_widths[j], precision, mat->mat[i][j], (j < mat->cols - 1) ? ", " : "");
    }
    printf(" ]\n");
  }
}