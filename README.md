# Final project report

Alexander Knysh, [github repository](https://github.com/alexanderknysh/hpc-final-project)

## Introduction

Solving of the so-called optimization problems usually requires finding minimum or maximum (generally speaking - extremum) of a specific function of one, two or more independent variables. Sometimes, it is mathematically challenging or even impossible to find neither local nor global extremums of the function analytically. In these cases, we use numerical techniques such as Newton's method, gradient descent, bi- and golden-section searches, etc. Unfortunately, they might be not only computationally expensive but have a hard time finding global extremum. Under such circumstances methods like [Monte Carlo method](https://en.wikipedia.org/wiki/Monte_Carlo_method) can be very useful and efficient, since it is based on different principles.

In this project the C program for finding global minimum of a given real function

<a href="https://www.codecogs.com/eqnedit.php?latex=f\left&space;(&space;x,y&space;\right&space;),&space;\quad&space;x&space;\in&space;\left&space;[&space;x_0,&space;x_1&space;\right&space;],&space;\quad&space;y&space;\in&space;\left&space;[&space;y_0,&space;y_1&space;\right&space;]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f\left&space;(&space;x,y&space;\right&space;),&space;\quad&space;x&space;\in&space;\left&space;[&space;x_0,&space;x_1&space;\right&space;],&space;\quad&space;y&space;\in&space;\left&space;[&space;y_0,&space;y_1&space;\right&space;]" title="f\left ( x,y \right ), \quad x \in \left [ x_0, x_1 \right ], \quad y \in \left [ y_0, y_1 \right ]" /></a>

with Monte Carlo method is presented. The MPI parallelization is implemented by [domain decomposition](https://github.com/unh-hpc-2019/class/wiki/Class-20#domain-decomposition-1). The input data are the function itself, its domain parameters and total number of points checked, while the output data are the global minimum of the function with corresponding coordinates. Note that the global maximum could be easily found by changing the sign of the function.   

## Numerical approach

The numerical approach is relatively simple and straightforward. First of all, we divide our rectangular domain along X-axis by several subdomains that correspond a particular MPI process (Fig. 1). Every MPI process generate its own random X- and Y-coordinates and calculate function at this point. Then, obtained value of function is compared to the already assigned initial minimum which is value of function at random point too. The procedure repeats finite number of times prescribed by user. Having the "local" minimums from each subdomain with corresponding coordinates, we can easily extract the global minimum through the communications between MPI processes.

![Problem's domain decomposition with MPI.](https://user-images.githubusercontent.com/46943028/57386810-24966000-7183-11e9-9651-8fad0e178cf2.PNG)

*Figure 1. Domain decomposition with MPI.*

## Description of the code

There are three files included in repository: header *random.h* that includes global variables and neccessary functions, main program *monte_carlo.c* and *CMakeLists.txt* for project building. Let us look at them in detail.

[**random.h**](https://github.com/alexanderknysh/hpc-final-project/blob/master/random.h)

In the header file you can find not only familiar `mprintf()` and `MHERE` structures but also the function to be analized and domain parameters. Parameter's type is `double` but `extern const double` could be a good programming practice too. Capitalization in variable names is forced due to reserved `y0()` and `y1()` functions of *math.h*.  
```c
#ifndef MPI_DEBUG_H
#define MPI_DEBUG_H

#include <mpi.h>

#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)
#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)

// Endpoints of the domain intervals.
double X0 = -1,
       X1 =  1,
       Y0 = -1,
       Y1 =  1;

// Given function.
double f(double x, double y)
{
  return x*x+y*y;
}
```
The `random_seed()` function seeds the random number generator with the current time and rank (MPI process number). It should be called once prior to calling any other random number function. Note that this function must depend on rank in some way, otherwise all processes will generate the same random numbers.
```c
// Seeds random based on time and rank.
void random_seed(int rank)
{
  srand(time(NULL) + rank);
}
```
Finally, the `random_real()` returns a random (preudorandom, to be exact) real number between `low` and `high`, inclusive.
```c
// Create random real from [low, high] interval.
double random_real(double low, double high)
{
  double d;

  d = (double) rand() / ((double) RAND_MAX + 1);
  return (low + d * (high - low));
}

#endif
```

[**monte_carlo.c**](https://github.com/alexanderknysh/hpc-final-project/blob/master/monte_carlo.c)

As it has been already mentioned, this code calculates global minimum with corresponding coordinates. Here is a set of libraries/headers included, declarations of MPI's `rank`, `size` and declaration of `tbeg` time variable that will be used to find execution time:
```c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
```
In order to correctly identify execution time, `MPI_Barrier(MPI_COMM_WORLD)` is placed right before time variable initialization. Next piece of code double check if the number of points we input satisfies general requirenments.
```c
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
```
For example, if the input is `mpirun -n 8 16000`, where 16000 is the total number of points, no error will be shown and each subdomain will test 2000 random points. However, if number of points is not chosen, not positive, not integer or not divisible by number of processes, we will get an error.

Now we can seed random number generator and assign our local `local.fmin`, `local.rmin`, `local.xmin`, `local.ymin` variables:
```c
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
```
The reason behind using of `struct` is that it is quite convenient to way to put data in `Allreduce()` which will be showm later. So now we just go through the prescribed number of points, find local minimum and print it with all corresponding parameters:
```c
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
```
Finally, we can go through all `local` variables, find the lowest `local.fmin` and rewrite it in `global.fmin`. We do it with a following `MPI_Allreduce()` command:
```c
  // Find global minimum from the local ones. Save corresponding MPI process rank.
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
```
This will not only save `global.fmin` value, but also remember the corresponding process number `global.rmin` by putting `MPI_MINLOC` instead of just `MPI_MIN`. Thus, the very last thing to do is to send all neccessary information to the first (`rank == 0`) process 
```c
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
```
and print global results and execution time:
```c
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
```

[**CMakeLists.txt**](https://github.com/alexanderknysh/hpc-final-project/blob/master/CMakeLists.txt)

Nothing special here, just add MPI package and *math.h* library, and build the code:
```c
cmake_minimum_required (VERSION 3.9)

project (MPI)

find_package(MPI REQUIRED)

set(CMAKE_C_STANDARD 99)
link_libraries(m)

add_executable(monte_carlo monte_carlo.c)
target_link_libraries(monte_carlo PUBLIC MPI::MPI_C)
```

## How to run

Make sure that MPI is supported by the compiler you use. In order to run code clone project repository, go in project folder and then build, compile and run:
```
$ git clone git@github.com:alexanderknysh/hpc-final-project.git
$ cd hpc-final-project
$ cmake .
$ make
$ mpirun -n 8 ./monte_carlo 16000
```
This will run code with 8 processes and 16000 points for the whole domain, so it is going to be 2000 points per process.

## Results

Let us first run the sample code included in repository on *fishercat.sr.unh.edu*. The chosen function and domain (Fig. 2) are

<a href="https://www.codecogs.com/eqnedit.php?latex=f\left&space;(&space;x,y&space;\right&space;)=x^2&plus;y^2,&space;\quad&space;x&space;\in&space;\left&space;[&space;-1,&space;1&space;\right&space;],&space;\quad&space;y&space;\in&space;\left&space;[&space;-1,&space;1&space;\right&space;]," target="_blank"><img src="https://latex.codecogs.com/gif.latex?f\left&space;(&space;x,y&space;\right&space;)=x^2&plus;y^2,&space;\quad&space;x&space;\in&space;\left&space;[&space;-1,&space;1&space;\right&space;],&space;\quad&space;y&space;\in&space;\left&space;[&space;-1,&space;1&space;\right&space;]," title="f\left ( x,y \right )=x^2+y^2, \quad x \in \left [ -1, 1 \right ], \quad y \in \left [ -1, 1 \right ]," /></a>

and the output for 8,000,000 points is:
```
[aknysh@fishercat build]$ mpirun -n 8 ./monte_carlo 8000000
[1] Local minimum at X ϵ [-0.75, -0.5], Y ϵ [-1, 1] among 1000000 points:  F(-0.500041, -5.84312e-05) = 0.250041
[3] Local minimum at X ϵ [-0.25, 0], Y ϵ [-1, 1] among 1000000 points:  F(-0.000420367, -2.48933e-05) = 1.77328e-07
[5] Local minimum at X ϵ [0.25, 0.5], Y ϵ [-1, 1] among 1000000 points:  F(0.250011, -0.00183207) = 0.0625091
[0] Local minimum at X ϵ [-1, -0.75], Y ϵ [-1, 1] among 1000000 points:  F(-0.750038, -0.00254893) = 0.562563
[2] Local minimum at X ϵ [-0.5, -0.25], Y ϵ [-1, 1] among 1000000 points:  F(-0.25006, 0.00440658) = 0.0625495
[4] Local minimum at X ϵ [0, 0.25], Y ϵ [-1, 1] among 1000000 points:  F(0.000257945, -0.00046741) = 2.85008e-07
[7] Local minimum at X ϵ [0.75, 1], Y ϵ [-1, 1] among 1000000 points:  F(0.750031, -0.00331392) = 0.562557
[6] Local minimum at X ϵ [0.5, 0.75], Y ϵ [-1, 1] among 1000000 points:  F(0.500025, 0.00241739) = 0.250031
----------------------------------------------------------- -- -
[3] Global minimum at X ϵ [-1,1], Y ϵ [-1,1] among 8000000 points:  F(-0.000420367, -2.48933e-05) = 1.77328e-07
----------------------------------------------------------- -- -
    Execution time: 0.103079 [s]
```
It is obvious that minimum function value is 0 at (0,0) point and either `[3]`<sup>rd</sup> or `[4]`<sup>th</sup> process had to report global minimum.

![x^2 + y^2](https://user-images.githubusercontent.com/46943028/57400548-b9a85180-71a1-11e9-8411-b9c29c660e4d.PNG)

*Figure 2. Plot of sample function.*

The reported solution is correct. Now we change function and domain (Fig. 3) to

<a href="https://www.codecogs.com/eqnedit.php?latex=f\left&space;(&space;x,y&space;\right&space;)=0.1x^2y&space;\enskip&space;sin\left&space;(&space;xy&space;\right&space;),&space;\quad&space;x&space;\in&space;\left&space;[&space;1,&space;2&space;\right&space;],&space;\quad&space;y&space;\in&space;\left&space;[&space;-3,&space;-2&space;\right&space;]," target="_blank"><img src="https://latex.codecogs.com/gif.latex?f\left&space;(&space;x,y&space;\right&space;)=0.1x^2y&space;\enskip&space;sin\left&space;(&space;xy&space;\right&space;),&space;\quad&space;x&space;\in&space;\left&space;[&space;1,&space;2&space;\right&space;],&space;\quad&space;y&space;\in&space;\left&space;[&space;-3,&space;-2&space;\right&space;]," title="f\left ( x,y \right )=0.1x^2y \enskip sin\left ( xy \right ), \quad x \in \left [ 1, 2 \right ], \quad y \in \left [ -3, -2 \right ]," /></a>

and run program for 60,000,000 points
```
[aknysh@fishercat build]$ mpirun -n 6 ./monte_carlo 60000000
[0] Local minimum at X ϵ [1, 1.16667], Y ϵ [-3, -2] among 10000000 points:  F(1.16666, -2.99986) = -0.143159
[2] Local minimum at X ϵ [1.33333, 1.5], Y ϵ [-3, -2] among 10000000 points:  F(1.49988, -2.99999) = -0.659676
[1] Local minimum at X ϵ [1.16667, 1.33333], Y ϵ [-3, -2] among 10000000 points:  F(1.33329, -2.99991) = -0.403508
[3] Local minimum at X ϵ [1.5, 1.66667], Y ϵ [-3, -2] among 10000000 points:  F(1.66666, -2.94859) = -0.802408
[4] Local minimum at X ϵ [1.66667, 1.83333], Y ϵ [-3, -2] among 10000000 points:  F(1.83332, -2.67908) = -0.882646
[5] Local minimum at X ϵ [1.83333, 2], Y ϵ [-3, -2] among 10000000 points:  F(2, -2.45671) = -0.962893
----------------------------------------------------------- -- -
[5] Global minimum at X ϵ [1,2], Y ϵ [-3,-2] among 60000000 points:  F(2, -2.45671) = -0.962893
----------------------------------------------------------- -- -
    Execution time: 4.63637 [s]

```

![0.1 x^2 y Sin(xy)](https://user-images.githubusercontent.com/46943028/57400734-291e4100-71a2-11e9-912c-f15e148d3a85.PNG)

*Figure 3. Plot of the changed function.*

The real minimum of the function is around -0.962892 at (2, -2.4566) and expected to be found by `[5]`<sup>th</sup> process. This is exactly what we get.

## Parallelization

For the strong scaling test we choose the changed function and domain again

<a href="https://www.codecogs.com/eqnedit.php?latex=f\left&space;(&space;x,y&space;\right&space;)=0.1x^2y&space;\enskip&space;sin\left&space;(&space;xy&space;\right&space;),&space;\quad&space;x&space;\in&space;\left&space;[&space;1,&space;2&space;\right&space;],&space;\quad&space;y&space;\in&space;\left&space;[&space;-3,&space;-2&space;\right&space;]," target="_blank"><img src="https://latex.codecogs.com/gif.latex?f\left&space;(&space;x,y&space;\right&space;)=0.1x^2y&space;\enskip&space;sin\left&space;(&space;xy&space;\right&space;),&space;\quad&space;x&space;\in&space;\left&space;[&space;1,&space;2&space;\right&space;],&space;\quad&space;y&space;\in&space;\left&space;[&space;-3,&space;-2&space;\right&space;]," title="f\left ( x,y \right )=0.1x^2y \enskip sin\left ( xy \right ), \quad x \in \left [ 1, 2 \right ], \quad y \in \left [ -3, -2 \right ]," /></a>

and run it on 1,2,...,8 processes and 84,000,000 points (Fig. 4).

![Strong scaling](https://user-images.githubusercontent.com/46943028/57406229-3e00d180-71ae-11e9-8f60-e78c03740142.PNG)

*Figure 4. Results of the strong scaling test. Dashed line represents perfect scaling, black points represent actual timings.*

For the weak scaling the number of points is proportional to the number of processes involved. In the perfect situation we will get a straight horizontal line, but in reality there is a slight deviation (Fig. 5).

![Weak scaling](https://user-images.githubusercontent.com/46943028/57408566-f3825380-71b3-11e9-8fe5-654aac1671bd.PNG)

*Figure 5. Results of the weak scaling test.*

