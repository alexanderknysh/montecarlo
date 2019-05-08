# Final project report

Alexander Knysh, [github repository](https://github.com/alexanderknysh/hpc-final-project)

## Introduction

Solving of so-called optimization problems usually requires finding minimum or maximum (generally speaking - extremum) of a specific function of one, two and more independent variables. Sometimes, it is mathematically challenging or even impossible to find neither local nor global extremums of the function analytically. In these cases, we use numerical techniques such as Newton's method, gradient descent, bi- and golden-section searches, etc. Unfortunately, they might be not only computationally expensive but have a hard time finding global extremum. Under such circumstances methods like [Monte Carlo method](https://en.wikipedia.org/wiki/Monte_Carlo_method) can be very useful and efficient, since it is based on different principles.

In this project the C program for finding global minimum of a given real function

<a href="https://www.codecogs.com/eqnedit.php?latex=f\left&space;(&space;x,y&space;\right&space;),&space;\quad&space;x&space;\in&space;\left&space;[&space;x_0,&space;x_1&space;\right&space;],&space;\quad&space;y&space;\in&space;\left&space;[&space;y_0,&space;y_1&space;\right&space;]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f\left&space;(&space;x,y&space;\right&space;),&space;\quad&space;x&space;\in&space;\left&space;[&space;x_0,&space;x_1&space;\right&space;],&space;\quad&space;y&space;\in&space;\left&space;[&space;y_0,&space;y_1&space;\right&space;]" title="f\left ( x,y \right ), \quad x \in \left [ x_0, x_1 \right ], \quad y \in \left [ y_0, y_1 \right ]" /></a>

with Monte Carlo method is presented. The MPI parallelization is implemented by [domain decomposition](https://github.com/unh-hpc-2019/class/wiki/Class-20#domain-decomposition-1). The input data are the function itself, its domain parameters and total number of points checked, while the output data are the global minimum of the function with corresponding coordinates. Note that the global maximum could be easily found by changing the sign of the function.   

## Numerical approach

The numerical approach is relatively simple and straightforward. First of all, we divide our rectangular domain along X-axis by several subdomains that correspond a particular MPI process, see Figure 1. Every MPI process generate its own random X- and Y-coordinates and calculate function at this point. Then, obtained value of function is compared the already assigned initial minimum which is value of function at random point too. The procedure repeats finite number of times prescribed by user. At this stage we have "local" minimums from each subdomain with corresponding coordinates, so the only thing left is to extract the global minimum through MPI processes communication.

![Problem's domain decomposition with MPI.](https://user-images.githubusercontent.com/46943028/57386810-24966000-7183-11e9-9651-8fad0e178cf2.PNG)

Figure 1. Problem's domain decomposition with MPI.

## Description of the code

## How to run

## Results

## Parallelization
