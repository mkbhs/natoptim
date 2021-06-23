# NatOptim
## Nature-Inspired Optimization Library for C++
NatOptim is an open source C++ library for solving optimization problems by using nature-inspired algorithms. It can be used to solve general non-linear unconstrained optimization problems.

### Installation

Install CMake
```
$ sudo apt-get install cmake
```
Clone this git repository in your $HOME directory
```
$ git clone https://github.com/mkbhs/natoptim.git
```
Build NatOptim:
```
$ cd natoptim
$ mkdir build
$ cd build
$ cmake ..
$ make
```
### How to include NatOptim in your project:
In your project directory create a CMakeLists.txt: 
```
cmake_minimum_required(VERSION 3.5)

project(my_project)

find_library(NatOptim_LIB libNatOptim.so PATHS $ENV{HOME}/natoptim/build/)

include_directories($ENV{HOME}/natoptim/include/)

add_executable(my_project my_project.cpp)

target_link_libraries(my_project ${NatOptim_LIB})
```

### Example using NatOptim Library:
```cpp
#include <iostream>
#include <stdio.h>
#include <bits/stdc++.h>
#include "randsearch.h"

using namespace std;

// Define the cost function to be optimized
double costFunc(double *x)
{
  unsigned int n = sizeof(x)/sizeof(x[0]);
  for(unsigned int i=0; i < n; i++)
  {
    double result = 0.0;
    result = (x[0]*x[0]) + (x[1]*x[1]);
    return result;
  }
}


int main()
{
  srand(time(0));

  // Define the algorithm's parameters
  int num_dim = 2;                                      // Search space dimensions
  double constrains[4] = {-10.0, 10.0, -10.0, 10.0};    // Dimension limits
  double* ptr_constrains = &constrains[0];
  int type_sol = 0;                                     // Minimization
  int max_iter = 500;                                   // Maximum number of iterations
  double sigma = 5.0;                                   // Standard deviation
  int num_workers = 30;                                 // Number of solutions per iteration

  // Create a random search optimization problem
  nptm::RandSearch my_random_search(num_dim, ptr_constrains, type_sol, max_iter, sigma, num_workers);

  // Provide the user defined cost function to the optimization problem
  my_random_search.cost_func = &costFunc;

  // Run the optimization process
  my_random_search.optimizeSolution();

  return 0;
}
```
### Cite Us
If you use NatOptim for a publication, please cite it as:
```
@misc{natoptim,
  author = "Sasha Vila and Others",
  title = "NatOptim",
  howpublished = "\url{https://github.com/mkbhs/natoptim}",
}
```


