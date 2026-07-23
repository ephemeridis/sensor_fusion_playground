#include "path.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PATH_MAX 1024
#define PATH_STOP_TOKEN "STOP"
#define PATH_RESUME_TOKEN "RESUME"

void __attribute__((weak)) path_error_handler(const char *function, const char *msg) {
  fprintf(stderr, "ERR  (%s): %s\n", function, msg);
}

void __attribute__((weak)) path_warning_handler(const char *function, const char *msg) {
  fprintf(stderr, "WARN (%s): %s\n", function, msg);
}

// pull the next comma-separated field out of the current strtok stream.
// returns false on a missing field as well as an unparseable one, so a
// truncated line can no longer reach PATH_STRTOD with a NULL pointer.
static bool next_num(path_scalar_t *out) {
  char *tok = strtok(NULL, ",");
  if (tok == NULL) return false;
  char *end;
  *out = PATH_STRTOD(tok, &end);
  return end != tok;
}

// grow all three point arrays to n entries. each pointer is committed to the
// segment as soon as its own realloc succeeds, so a later failure never
// leaves a stale pointer to freed memory for path_free to free again.
static bool seg_grow(segment_t *seg, size_t n) {
  path_scalar_t *nt = realloc(seg->t, sizeof(path_scalar_t) * n);
  if (nt == NULL) return false;
  seg->t = nt;
  path_scalar_t *nx = realloc(seg->x, sizeof(path_scalar_t) * n);
  if (nx == NULL) return false;
  seg->x = nx;
  path_scalar_t *ny = realloc(seg->y, sizeof(path_scalar_t) * n);
  if (ny == NULL) return false;
  seg->y = ny;
  return true;
}

void path_load_fromfile(path_t *path, const char *filename) {
  if (path == NULL) {
    path_error_handler(__func__, "path is NULL");
    return;
  }
  path->n_segments = 0;
  path->segments = NULL;

  FILE *fptr = fopen(filename, "r");
  if (fptr == NULL) {
    path_error_handler(__func__, "could not open file");
    return;
  }
  size_t alloc_segments = 0;
  size_t alloc_points = 0;
  segment_t *curr_seg = NULL;

  char buffer[PATH_MAX];
  char *curr;
  while (fgets(buffer, PATH_MAX, fptr) != NULL) {
    size_t len = strlen(buffer);

    if (len > 0 && buffer[len-1] != '\n' && !feof(fptr)) {
      path_warning_handler(__func__, "line too long, skipping");
      int c;
      while ((c = fgetc(fptr)) != '\n' && c != EOF) { }
      continue;
    }

    // strip comments
    strtok(buffer, "#");

    path_scalar_t x = 0, y = 0, t = 0;
    char *endptr;
    bool resume = false;
    bool stop = false;
    // check if string contains special lines
    if (len >= strlen(PATH_STOP_TOKEN)
      && memcmp(buffer, PATH_STOP_TOKEN, strlen(PATH_STOP_TOKEN)) == 0) {
      // stop token detected
      curr = buffer + strlen(PATH_STOP_TOKEN);
      stop = true;
    } else if (len >= strlen(PATH_RESUME_TOKEN)
      && memcmp(buffer, PATH_RESUME_TOKEN, strlen(PATH_RESUME_TOKEN)) == 0) {
      // resume token detected
      curr = buffer + strlen(PATH_RESUME_TOKEN);
      resume = true;
    } else {
      // process as a regular line
      curr = buffer;
    }

    curr = strtok(curr, ",");
    if (curr == NULL) continue;
    t = PATH_STRTOD(curr, &endptr);
    if (endptr == curr) continue;
    if (!resume) {
      if (!next_num(&x)) continue;
      if (!next_num(&y)) continue;
    } else {
      // if we are resuming, begin the X, Y coordinates from the point where the last segment ended
      if (path->n_segments == 0) {
        path_error_handler(__func__, "invalid RESUME token");
        path_free(path); fclose(fptr);
        return;
      }
      segment_t *prev = &path->segments[path->n_segments-1];
      if (prev->n_points == 0) {
        path_error_handler(__func__, "RESUME from empty segment");
        path_free(path); fclose(fptr);
        return;
      }
      x = prev->x[prev->n_points-1];
      y = prev->y[prev->n_points-1];
    }

    if (curr_seg == NULL) {
      // start a new path segment
      // check if we need to reallocate the number of segments
      size_t realloc_n_segments = alloc_segments;
      if (alloc_segments == 0) realloc_n_segments = 1;
      else if (path->n_segments+1 >= alloc_segments) realloc_n_segments *= 2;
      if (realloc_n_segments != alloc_segments) {
        segment_t *new_segments = realloc(path->segments, sizeof(segment_t) * realloc_n_segments);
        if (new_segments == NULL) {
          path_error_handler(__func__, "could not reallocate segments");
          path_free(path); fclose(fptr);
          return;
        }

        path->segments = new_segments;
        alloc_segments = realloc_n_segments;
      }
      // set current working segment
      curr_seg = &path->segments[path->n_segments];
      path->n_segments++;

      // set number of allocated points to zero
      alloc_points = 0;
      curr_seg->n_points = 0;
      curr_seg->t = NULL;
      curr_seg->x = NULL;
      curr_seg->y = NULL;
    }

    size_t realloc_n_pts = alloc_points;
    if (alloc_points == 0) realloc_n_pts = 1;
    else if (curr_seg->n_points+1 >= alloc_points) realloc_n_pts *= 2;
    if (realloc_n_pts != alloc_points) {
      if (!seg_grow(curr_seg, realloc_n_pts)) {
        path_error_handler(__func__, "could not reallocate points in segment");
        path_free(path); fclose(fptr);
        return;
      }
      alloc_points = realloc_n_pts;
    }

    // add point to segment
    curr_seg->t[curr_seg->n_points] = t;
    curr_seg->x[curr_seg->n_points] = x;
    curr_seg->y[curr_seg->n_points] = y;
    curr_seg->n_points++;

    if (stop) {
      // end the current segment
      // trim the spare capacity
      if (!seg_grow(curr_seg, curr_seg->n_points)) {
        path_error_handler(__func__, "could not shrink points in segment");
        path_free(path); fclose(fptr);
        return;
      }
      curr_seg = NULL;
    }

  }

  if (path->n_segments > 0) {
    segment_t *new_segments = realloc(path->segments, sizeof(segment_t) * path->n_segments);
    if (new_segments == NULL) {
      path_error_handler(__func__, "could not reallocate segments");
      path_free(path); fclose(fptr);
      return;
    }
    path->segments = new_segments;
  }

  fclose(fptr);
}

void path_free(path_t *path) {
  if (path == NULL) return;
  if (path->segments) {
    for (size_t i = 0; i < path->n_segments; i++) {
      free(path->segments[i].x);
      free(path->segments[i].y);
      free(path->segments[i].t);
    }
    free(path->segments);
  }
  // leave the struct safe to free again
  path->segments = NULL;
  path->n_segments = 0;
}

void path_print(path_t *path) {
  path_print_prec(path, 3);
}

void path_print_prec(path_t *path, const int precision) {
  if (path == NULL) {
    path_error_handler(__func__, "path is NULL");
    return;
  }

  if (path->segments == NULL) {
    path_error_handler(__func__, "path segments is NULL");
    return;
  }

  // do a first pass first
  // to find the widest formatted value in each column
  int w_t = 0, w_x = 0, w_y = 0;
  for (size_t i = 0; i < path->n_segments; i++) {
    segment_t *seg = &path->segments[i];
    if (seg->t == NULL || seg->x == NULL || seg->y == NULL) {
      path_error_handler(__func__, "path segment points is NULL");
      return;
    }
    for (size_t j = 0; j < seg->n_points; j++) {
      int lt = snprintf(NULL, 0, "%.*lf", precision, (double)seg->t[j]);
      int lx = snprintf(NULL, 0, "%.*lf", precision, (double)seg->x[j]);
      int ly = snprintf(NULL, 0, "%.*lf", precision, (double)seg->y[j]);
      if (lt > w_t) w_t = lt;
      if (lx > w_x) w_x = lx;
      if (ly > w_y) w_y = ly;
    }
  }

  // print with the computed widths
  printf("PATH:\n");
  for (size_t i = 0; i < path->n_segments; i++) {
    segment_t *seg = &path->segments[i];
    printf("  SEGMENT %zu\n", i);
    for (size_t j = 0; j < seg->n_points; j++) {
      printf("    t: %*.*lf, x: %*.*lf, y: %*.*lf\n",
             w_t, precision, (double)seg->t[j],
             w_x, precision, (double)seg->x[j],
             w_y, precision, (double)seg->y[j]);
    }
  }
}