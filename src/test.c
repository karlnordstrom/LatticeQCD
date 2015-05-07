#include "metropolis.h"
#include "stdio.h"
#include "time.h"

int main(void) {
  init_genrand64(time(NULL));
  printf("Attempting to allocate 4D grid... \n");
  Grid grid = allocateGrid(8);
  printf("Successfully allocated grid. \n");
  setRandomSU3(findLink(grid, 0, 0, 0, 0, 0), 0.4);
  printMatrix(findLink(grid, 0, 0, 0, 0, 0));
  MatrixSU3 inverse = invert(findLink(grid, 8, 0, 0, 0, 0));
  printf("%f + %f i \n", creal(determinant(inverse)), cimag(determinant(inverse)));
  printMatrix(inverse);
  multiply(inverse, findLink(grid, 8, 0, 0, 0, 0));
  printMatrix(inverse);
  printf("%f \n", averagePlaquette(grid, 5.5, 0.24, 50, 200));
  return 0;
}
