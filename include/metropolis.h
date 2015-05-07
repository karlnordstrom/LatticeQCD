#include "matrixlib.h"

typedef struct Grid {
  int *grid;
  MatrixSU3 *links;
  const size_t gridsize;
} Grid;

Grid allocateGrid(const size_t gridsize);

int findGridPoint(const Grid grid, const size_t x, const size_t y, const size_t z, const size_t t);
MatrixSU3 findLink(const Grid grid, const size_t x1, const size_t y1, const size_t z1, const size_t t1, const size_t dir); // 0 = x, 1 = y, 2 = z, 3 = t
double plaquetteAction(const Grid grid, const size_t x, const size_t y, const size_t z, const size_t t);
void updateMetropolis(const Grid grid, double beta, double epsilon);
double averagePlaquette(const Grid grid, double beta, double epsilon, size_t Ncor, size_t Ncf);
double averagePlaquetteNow(const Grid grid, double beta);
