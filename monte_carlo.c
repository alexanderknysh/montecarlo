
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <time.h>
#include "random.h"

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  int n = atoi(argv[1]);
  double n_double = atof(argv[1]);

  if (argc != 2) {
    mprintf("Error: number of random points is not chosen.\n");
    exit(1);
  }
  else if (n < 1 || n != n_double) {
    mprintf("Error: number of random points is not positive integer.\n");
    exit(1);
  }
  else if (n % size != 0) {
    mprintf("Error: number of random points is not divisible by number of processes.\n");
    exit(1);
  }
  
  random_seed(rank);
  
  double x0 = -1, x1 = 1, y0 = -1, y1 = 1;
  double dx = (x1 - x0) / size;
  double xb = x0 + dx * rank;
  double xe = xb + dx;
  
  double xmin = random_real(xb, xe);
  double ymin = random_real(y0, y1);
  double fmin = f(xmin, ymin);
  double xcur, ycur;

  for (int i = 1; i < n/size; i++) {
    xcur = random_real(xb,xe);
    ycur = random_real(y0,y1);
    if (f(xcur,ycur) < f(xmin,ymin)) {
      xmin = xcur;
      ymin = ycur;
      fmin = f(xcur,ycur);
    }
  }
  mprintf("Local minimum out of %d points: f(%g, %g) = %g\n", n/size, xmin, ymin, fmin);
 
  struct {
    double minvalue;
    int    minrank;
  } local, global; 
 
  local.minvalue = fmin;
  local.minrank = rank;

  MPI_Reduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    mprintf("Global minimum: %g Corresponding process: %d\n", global.minvalue, global.minrank);
  }

  MPI_Finalize();
  return 0;
}
