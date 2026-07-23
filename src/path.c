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


void path_load_fromfile(path_t *path, const char *filename) {
  FILE *fptr = fopen(filename, "r");
  if (fptr == NULL) {
    path_error_handler(__func__, "could not open file");
    return;
  }

  if (path == NULL) {
    path_error_handler(__func__, "could not allocate path");
    return;
  }
  path->n_segments = 0;
  path->segments = NULL;
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

    double x, y, t;
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
    t = strtod(curr, &endptr);
    if (endptr == curr) continue;
    if (!resume) {
      curr = strtok(NULL, ",");
      x = strtod(curr, &endptr);
      if (endptr == curr) continue;
      curr = strtok(NULL, ",");
      y = strtod(curr, &endptr);
      if (endptr == curr) continue;
    } else {
      // if we are resuming, begin the X, Y coordinates from the point where the last segment ended
      if (path->n_segments <= 0) {
        path_error_handler(__func__, "invalid RESUME token");
        path_free(path); fclose(fptr);
        return;
      }
      size_t last = path->segments[path->n_segments-1].n_points-1;
      x = path->segments[path->n_segments-1].x[last];
      y = path->segments[path->n_segments-1].y[last];
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
      double *new_t = realloc(curr_seg->t, sizeof(double) * realloc_n_pts);
      double *new_x = realloc(curr_seg->x, sizeof(double) * realloc_n_pts);
      double *new_y = realloc(curr_seg->y, sizeof(double) * realloc_n_pts);
      if (new_t == NULL || new_x == NULL || new_y == NULL) {
        path_error_handler(__func__, "could not reallocate points in segment");
        path_free(path); fclose(fptr);
        if (new_t != NULL) free(new_t);
        if (new_x != NULL) free(new_x);
        if (new_y != NULL) free(new_y);
        return;
      }
      curr_seg->t = new_t;
      curr_seg->x = new_x;
      curr_seg->y = new_y;
      alloc_points = realloc_n_pts;
    }

    // add point to segment
    curr_seg->t[curr_seg->n_points] = t;
    curr_seg->x[curr_seg->n_points] = x;
    curr_seg->y[curr_seg->n_points] = y;
    curr_seg->n_points++;

    if (stop) {
      // end the current segment
      double *new_t = realloc(curr_seg->t, sizeof(double) * curr_seg->n_points);
      double *new_x = realloc(curr_seg->x, sizeof(double) * curr_seg->n_points);
      double *new_y = realloc(curr_seg->y, sizeof(double) * curr_seg->n_points);
      if (new_t == NULL || new_x == NULL || new_y == NULL) {
        path_error_handler(__func__, "could not reallocate points in segment");
        path_free(path); fclose(fptr);
        if (new_t != NULL) free(new_t);
        if (new_x != NULL) free(new_x);
        if (new_y != NULL) free(new_y);
        return;
      }
      curr_seg->t = new_t;
      curr_seg->x = new_x;
      curr_seg->y = new_y;
      curr_seg = NULL;
    }

  }

  segment_t *new_segments = realloc(path->segments, sizeof(segment_t) * path->n_segments);
  if (new_segments == NULL) {
    path_error_handler(__func__, "could not reallocate segments");
    path_free(path); fclose(fptr);
    return;
  }

  path->segments = new_segments;
  fclose(fptr);
}

void path_free(path_t *path) {
  if (path) {
    // path is allocated
    // loop through the segments to free all the points
    if (path->segments) {
      for (size_t i = 0; i < path->n_segments; i++) {
        if (path->segments[i].x) free(path->segments[i].x);
        if (path->segments[i].y) free(path->segments[i].y);
        if (path->segments[i].t) free(path->segments[i].t);
      }
      free(path->segments);
    }
  }
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
      int lt = snprintf(NULL, 0, "%.*lf", precision, seg->t[j]);
      int lx = snprintf(NULL, 0, "%.*lf", precision, seg->x[j]);
      int ly = snprintf(NULL, 0, "%.*lf", precision, seg->y[j]);
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
             w_t, precision, seg->t[j],
             w_x, precision, seg->x[j],
             w_y, precision, seg->y[j]);
    }
  }
}