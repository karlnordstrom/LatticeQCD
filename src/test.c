#include "metropolis.h"
#include "stdio.h"

int main(void) {
  printf("Attempting to allocate 4D grid... \n");
  Grid grid = allocateGrid(8);
  printf("Successfully allocated grid. \n");
  setRandomSU3(findGridPoint(grid, 0, 0, 0, 0), 0.4);
  printMatrix(findGridPoint(grid, 8, 0, 0, 0));
  return 0;
}
