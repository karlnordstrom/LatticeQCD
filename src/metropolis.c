#include "metropolis.h"

Grid allocateGrid(const size_t gridsize) {
  MatrixSU3 *matrix = (MatrixSU3 *) malloc ( gridsize * gridsize * gridsize * gridsize * sizeof(MatrixSU3 *) );
  size_t i = 0;
  while(i < gridsize * gridsize * gridsize * gridsize) { *(matrix + i) = allocateSU3(); ++i; }
  Grid grid = { matrix, gridsize };
  return grid;
}

MatrixSU3 findGridPoint(const Grid grid, const size_t x, const size_t y, const size_t z, const size_t t) {
  return *(grid.grid + (x % grid.gridsize) +  grid.gridsize * (y % grid.gridsize) + grid.gridsize * grid.gridsize * (z % grid.gridsize) + grid.gridsize * grid.gridsize * grid.gridsize * (t % grid.gridsize));
}
