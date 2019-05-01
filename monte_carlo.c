
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
  
  double x0 = 0, x1 = 1, y0 = 0, y1 = 1;
  double dx = (x1 - x0) / size;

  double xb = x0 + dx * rank;
  double xe = xb + dx;
  
  mprintf("XB = %g\n", xb);
  mprintf("XE = %g\n", xe);
  
  double xmin = random_real(xb, xe);
  double ymin = random_real(y0, y1);
  double fmin = f(xmin, ymin);
  double xcur, ycur;

  mprintf("per-proc xmin, ymin, fmin are %g , %g, %g\n", xmin, ymin, fmin);
  for (int i = 1; i < n; i++) {
    xcur = random_real(xb,xe);
    ycur = random_real(y0,y1);
    if (f(xcur,ycur) < f(xmin,ymin)) {
      xmin = xcur;
      ymin = ycur;
      fmin = f(xcur,ycur);
    }
  }
  mprintf("per-proc xmin, ymin, fmin are %g , %g, %g\n", xmin, ymin, fmin);

  double global_fmin;
  MPI_Reduce(&fmin, &global_fmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    mprintf("The global minimum is %g\n", global_fmin);
  }

  MPI_Finalize();
  return 0;
}
