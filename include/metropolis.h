#include "matrixlib.h"

typedef struct Grid {
  MatrixSU3 *grid;
  const size_t gridsize;
} Grid;

Grid allocateGrid(const size_t gridsize);

MatrixSU3 findGridPoint(const Grid grid, const size_t x, const size_t y, const size_t z, const size_t t);
