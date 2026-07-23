#ifndef SENSOR_FUSION_PLAYGROUND_MATRIX_H
#define SENSOR_FUSION_PLAYGROUND_MATRIX_H
#include <stdbool.h>
#include <math.h>
#include <float.h>

#ifdef MAT_DOUBLE_PRECISION
  typedef double mat_scalar_t;
#define MAT_FABS    fabs
#define MAT_EPSILON DBL_EPSILON
#else
typedef float mat_scalar_t;
#define MAT_FABS    fabsf
#define MAT_EPSILON FLT_EPSILON
#endif

// swap scratch buffer size for
// chunked memory swapping when re-ordering rows in LUP pivoting
#ifndef MAT_SWAP_CHUNK_SIZE
#define MAT_SWAP_CHUNK_SIZE 32
#endif

typedef struct {
  int rows;
  int cols;
  bool heap_alloc;
  mat_scalar_t **mat;
  mat_scalar_t *_mat_data;
} mat_t;

// debug functions for printing the matrix
#ifndef MAT_NO_PRINT

// whether to use uniform sizing for all the cols or to
// dynamically size based on the number to print
#ifndef MAT_PRINT_COLS_UNIFORM_WIDTH
#ifndef MAT_MAX_PRINT_COLS
#define MAT_MAX_PRINT_COLS 32
#endif
#endif

void mat_print(mat_t *mat);
void mat_print_prec(mat_t *mat, int precision);

#endif

#define MAT_STATIC_STORAGE(name, R, C)          \
static mat_scalar_t name##_data[(R) * (C)];   \
static mat_scalar_t *name##_rows[(R)];        \
static mat_t name

#define MAT_STATIC_INIT(name, R, C) \
mat_init_static(&name, (R), (C), name##_data, name##_rows)

bool mat_init_static(mat_t *m, int rows, int cols, mat_scalar_t *data, mat_scalar_t **rowptrs);
mat_t *mat_new(int rows, int cols, mat_scalar_t fill_value);
mat_t *mat_zero(int rows, int cols);
mat_t *mat_identity(int N);
mat_t *mat_copy(mat_t *mat);

bool mat_fill(mat_t *mat, mat_scalar_t v);
bool mat_set_identity(mat_t *mat);
bool mat_copy_to(mat_t *dst, mat_t *src);

mat_t *mat_add(mat_t *mat1, mat_t *mat2);
mat_t *mat_sub(mat_t *mat1, mat_t *mat2);
mat_t *mat_mul(mat_t *mat1, mat_t *mat2);
mat_t *mat_mul_scalar(mat_t *mat, mat_scalar_t scalar);

bool mat_add_to(mat_t *dst, mat_t *mat1, mat_t *mat2);
bool mat_sub_to(mat_t *dst, mat_t *mat1, mat_t *mat2);
bool mat_mul_to(mat_t *dst,mat_t *mat1, mat_t *mat2);
bool mat_mul_scalar_to(mat_t *dst, mat_t *src, mat_scalar_t scalar);

mat_t *mat_transpose(mat_t *mat);
bool mat_transpose_to(mat_t *dst, mat_t *src);
bool mat_mul_transpose_to(mat_t *dst, mat_t *mat1, mat_t *mat2);
bool mat_transpose_inplace_square(mat_t *mat);

bool mat_lu_decomp(mat_t *mat, mat_t *l, mat_t *u);
bool mat_lup_decomp(mat_t *mat, mat_t *l, mat_t *u, mat_t *p);

bool mat_forward_backward(mat_t *x, mat_t *y, mat_t *l, mat_t *u);

void mat_free(mat_t *mat);

#endif //SENSOR_FUSION_PLAYGROUND_MATRIX_H
