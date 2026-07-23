#ifndef SENSOR_FUSION_PLAYGROUND_PATH_H
#define SENSOR_FUSION_PLAYGROUND_PATH_H
#include <stddef.h>

typedef struct {
  size_t n_points;
  double *t;
  double *x;
  double *y;
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
