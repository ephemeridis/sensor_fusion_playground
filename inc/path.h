#ifndef SENSOR_FUSION_PLAYGROUND_PATH_H
#define SENSOR_FUSION_PLAYGROUND_PATH_H
#include <stddef.h>
#include <stdlib.h>

// scalar type selection, independent of the matrix and bspline modules.
// define PATH_DOUBLE_PRECISION as a build flag (not in a source file)
// so every translation unit agrees on the layout of segment_t.
#ifdef PATH_DOUBLE_PRECISION
typedef double path_scalar_t;
#define PATH_STRTOD strtod
#else
typedef float path_scalar_t;
#define PATH_STRTOD strtof
#endif

typedef struct {
  size_t n_points;
  path_scalar_t *t;
  path_scalar_t *x;
  path_scalar_t *y;
} segment_t;

typedef struct {
  size_t n_segments;
  segment_t *segments;
} path_t;

void path_load_fromfile(path_t *path, const char *filename);
void path_free(path_t *path);

void path_print(path_t *path);
void path_print_prec(path_t *path, int precision);

#endif //SENSOR_FUSION_PLAYGROUND_PATH_H