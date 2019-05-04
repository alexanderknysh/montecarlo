
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
 
  MPI_Barrier(MPI_COMM_WORLD);
  double tbeg = MPI_Wtime();

  // Check if input total number of points is correct.
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
  
  // Seed the random number generator.
  random_seed(rank);
  
  // Initialize global and local boundaries.
  double dx = (X1 - X0) / size;
  double xb = X0 + dx * rank;
  double xe = xb + dx;
  
  // This structure keeps following local/global values:
  //
  // fmin - local/global minimum function value;
  // rmin - corresponding MPI process nubmer;
  // xmin - corresponding X coordinate;
  // ymin - corresponding Y coordinate.  
  struct {
    double fmin;
    int    rmin;
    double xmin;
    double ymin;
  } local, global;
  
  // Initialize local values.
  local.rmin = rank;
  local.xmin = random_real(xb, xe);
  local.ymin = random_real(Y0, Y1);
  local.fmin = f(local.xmin, local.ymin);

  // Find and reassign local minimum values. Print result.
  for (int i = 1; i < n/size; i++) {
    double xcur = random_real(xb, xe);
    double ycur = random_real(Y0, Y1);
    if (f(xcur, ycur) < f(local.xmin, local.ymin)) {
      local.xmin = xcur;
      local.ymin = ycur;
      local.fmin = f(xcur, ycur);
    }
  }
  mprintf("Local minimum at X ϵ [%g, %g], Y ϵ [%g, %g] among %d points:"
          "  F(%g, %g) = %g\n",
              xb, xe, Y0, Y1, n/size, local.xmin, local.ymin, local.fmin);

  // Find global minimum from the local ones. Save corresponding MPI process rank.
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
  
  // Send global minimum data to the first (0) MPI process.
  if (rank == global.rmin) {
    MPI_Send(&local.xmin, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    MPI_Send(&local.ymin, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
  }

  if (rank == 0) {
    MPI_Recv(&global.xmin, 1, MPI_DOUBLE, global.rmin, 1, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    MPI_Recv(&global.ymin, 1, MPI_DOUBLE, global.rmin, 2, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  }

  // Wait for all processes and print final result.
  MPI_Barrier(MPI_COMM_WORLD);
  double tend = MPI_Wtime();

  if (rank == 0) {
    printf("----------------------------------------------------------- -- -"
           "\n[%d] Global minimum at X ϵ [%g,%g], Y ϵ [%g,%g] among %d points:"
           "  F(%g, %g) = %g\n"
           "----------------------------------------------------------- -- -\n",
           global.rmin, X0, X1, Y0, Y1, n, global.xmin, global.ymin, global.fmin);
    printf("    Execution time: %g [s]\n", tend - tbeg);
  }
  
  MPI_Finalize();
  return 0;
}
