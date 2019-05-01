
#ifndef MPI_DEBUG_H
#define MPI_DEBUG_H

#include <mpi.h>

#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)
#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)

// Seeds random based on time and rank.
void random_seed(int rank)
{
  srand(time(NULL) + rank);
}

// Create random real from [low, high].
double random_real(double low, double high)
{
  double d;

  d = (double) rand() / ((double) RAND_MAX + 1);
  return (low + d * (high - low));
}

// Given function.
double
f(double x, double y)
{
  return x*x+y*y;
}

#endif
