#include <stdio.h>
#include <string.h>

#include "matrix.h"
#include "path.h"

int main(void) {
  path_t new_path;
  path_load_fromfile(&new_path, "./path");
  path_print(&new_path);

  mat_t *mat = mat_identity(5);
  mat_print(mat);
  return 0;
}
