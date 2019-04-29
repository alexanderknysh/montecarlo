
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <time.h>
#include "random.h"

void random_seed(int rank)
{
  srand(time(NULL) + rank);
}

double random_real(double low, double high)
{
  double d;

  d = (double) rand() / ((double) RAND_MAX + 1);
  return (low + d * (high - low));
}

double
f(double x, double y)
{
  return x*x+y*y;
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  random_seed(rank);

  int n = 1000;
  
  double x0 = 0, x1 = 1, y0 = 0, y1 = 1;
  double dx = (x1 - x0) / size;

  double xb = x0 + dx * rank;
  double xe = x0 + dx * (rank + 1);
  
  mprintf("XB = %g\n", xb);
  mprintf("XE = %g\n", xe);
  
  double xmin = random_real(xb, xe);
  double ymin = random_real(y0, y1);
  double fmin = f(xmin, ymin);
  double xcur, ycur;

  for (int i = 0; i < n; i++) {
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
  MPI_Reduce(&fmin, &global_fmin, 1, MPI_DOUBLE, MPI_MIN,
             0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    mprintf("The global minimum is %g\n", global_fmin);
  }

  MPI_Finalize();
  return 0;
}
