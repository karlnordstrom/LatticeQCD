#include "metropolis.h"

Grid allocateGrid(const size_t gridsize) {
  int *points = (int *) malloc (gridsize * gridsize * gridsize * gridsize * sizeof(int));
  MatrixSU3 *links = (MatrixSU3 *) malloc ( 4*gridsize * gridsize * gridsize * gridsize * sizeof(MatrixSU3 *) );
  for(size_t i = 0; i < gridsize * gridsize * gridsize * gridsize; ++i) {
    *(links + 4*i) = allocateSU3();
    *(links + 4*i + 1) = allocateSU3();
    *(links + 4*i + 2) = allocateSU3();
    *(links + 4*i + 3) = allocateSU3();
  }
  Grid grid = { points, links, gridsize };
  return grid;
}

int findGridPoint(const Grid grid, const size_t x, const size_t y, const size_t z, const size_t t) {
  return *(grid.grid + (x % grid.gridsize) +  grid.gridsize * (y % grid.gridsize) + grid.gridsize * grid.gridsize * (z % grid.gridsize) + grid.gridsize * grid.gridsize * grid.gridsize * (t % grid.gridsize));
}

MatrixSU3 findLink(const Grid grid, const size_t x, const size_t y, const size_t z, const size_t t, const size_t dir) {
  return *(grid.links + 4*((x % grid.gridsize) +  grid.gridsize * (y % grid.gridsize) + grid.gridsize * grid.gridsize * (z % grid.gridsize) + grid.gridsize * grid.gridsize * grid.gridsize * (t % grid.gridsize)) + dir);
}

double plaquetteAction(const Grid grid, const size_t x, const size_t y, const size_t z, const size_t t) {
  double action = 0.;
  MatrixSU3 matrix = allocateSU3();
  for(size_t i = 0; i < 4; ++i) {
    for(size_t j = i+1; j < 4; ++j) {
      setIdentity(matrix);
      size_t tmp_x = 0,tmp_y = 0,tmp_z = 0,tmp_t = 0,tmp_x1 = 0,tmp_y1 = 0,tmp_z1 = 0,tmp_t1 = 0;
      if(i == 0) tmp_x = 1;
      if(i == 1) tmp_y = 1;
      if(i == 2) tmp_z = 1;
      if(i == 3) tmp_t = 1;
      if(j == 1) tmp_y1 = 1;
      if(j == 2) tmp_z1 = 1;
      if(j == 3) tmp_t1 = 1;
      multiply(matrix, findLink(grid, x, y, z, t, i));
      multiply(matrix, findLink(grid, x+tmp_x, y+tmp_y, z+tmp_z, t+tmp_t, j));
      MatrixSU3 tmp = conjugate(findLink(grid, x+tmp_x1, y+tmp_y1, z+tmp_z1, t+tmp_t1, i));
      multiply(matrix, tmp);
      free(tmp);
      tmp = conjugate(findLink(grid, x, y, z, t, j));
      multiply(matrix, tmp);
      free(tmp);
      action += (1/3.)*creal(trace(matrix));
    }
  }
  free(matrix);
  return -action;
}

void updateMetropolis(const Grid grid, double beta, double epsilon) {
  MatrixSU3 *rands = (MatrixSU3 *) malloc ( 50 * sizeof(MatrixSU3 *) );
  MatrixSU3 *rands_conj = (MatrixSU3 *) malloc ( 50 * sizeof(MatrixSU3 *) );
  for(size_t r = 0; r < 50; ++r) { *(rands + r) = allocateRandomSU3(epsilon); *(rands_conj + r) = conjugate(*(rands + r)); }
  for(size_t x = 0; x < grid.gridsize; ++x) {
    for(size_t y = 0; y < grid.gridsize; ++y) {
      for(size_t z = 0; z < grid.gridsize; ++z) {
        for(size_t t = 0; t < grid.gridsize; ++t) {
          double old_S = beta*plaquetteAction(grid, x, y, z, t);
          int rand = floor(uniform(0,49.99));
          multiply(findLink(grid, x, y, z, t, 0), *(rands + rand));
          double dS = beta*plaquetteAction(grid, x, y, z, t) - old_S;
          if(dS > 0 && exp(-dS) < uniform(0,1)) multiply(findLink(grid, x, y, z, t, 0), *(rands_conj + rand));

          rand = floor(uniform(0,49.99));
          multiply(findLink(grid, x, y, z, t, 1), *(rands + rand));
          dS = beta*plaquetteAction(grid, x, y, z, t) - old_S;
          if(dS > 0 && exp(-dS) < uniform(0,1)) multiply(findLink(grid, x, y, z, t, 1), *(rands_conj + rand));

          rand = floor(uniform(0,49.99));
          multiply(findLink(grid, x, y, z, t, 2), *(rands + rand));
          dS = beta*plaquetteAction(grid, x, y, z, t) - old_S;
          if(dS > 0 && exp(-dS) < uniform(0,1)) multiply(findLink(grid, x, y, z, t, 2), *(rands_conj + rand));

          rand = floor(uniform(0,49.99));
          multiply(findLink(grid, x, y, z, t, 3), *(rands + rand));
          dS = beta*plaquetteAction(grid, x, y, z, t) - old_S;
          if(dS > 0 && exp(-dS) < uniform(0,1)) multiply(findLink(grid, x, y, z, t, 3), *(rands_conj + rand));
        }
      }
    }
  }
  free(rands); free(rands_conj);
}

double averagePlaquetteNow(const Grid grid, double beta){
  double ave = 0.;
  for(size_t x = 0; x < grid.gridsize; ++x) {
    for(size_t y = 0; y < grid.gridsize; ++y) {
      for(size_t z = 0; z < grid.gridsize; ++z) {
        for(size_t t = 0; t < grid.gridsize; ++t) {
          ave += beta*plaquetteAction(grid, x, y, z, t);
        }
      }
    }
  }
  return -ave/(grid.gridsize*grid.gridsize*grid.gridsize*grid.gridsize);
}

double averagePlaquette(const Grid grid, double beta, double epsilon, size_t Ncor, size_t Ncf) {
  for(size_t i = 0; i < 20*Ncor; ++i) updateMetropolis(grid, beta, epsilon);
  double averagePlaq = 0.;
  for(size_t it = 0; it < Ncf; ++it) {
    for(size_t i = 0; i < Ncor; ++i) updateMetropolis(grid, beta, epsilon);
//    printf("%f \n", averagePlaquetteNow(grid, beta));
    averagePlaq += averagePlaquetteNow(grid, beta);
  }
  return averagePlaq/Ncf;
}
