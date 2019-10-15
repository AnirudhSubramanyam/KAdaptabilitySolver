
k-adaptability-solver
====

k-adaptability-solver is a numerical optimization package, written in C++, for solving K-adaptability counterparts of two-stage mixed-integer robust optimization problems, based on our paper [K-Adaptability in Two-Stage Mixed-Integer Robust Optimization](https://arxiv.org/abs/1706.07097).

If you use this code in your research, please cite:
```
@article{sgw:19,
  title   = {K-adaptability in two-stage mixed-integer robust optimization},
  author  = {Subramanyam, Anirudh and Gounaris, Chrysanthos E and Wiesemann, Wolfram},
  journal = {Mathematical Programming Computation},
  year    = {2019}
}
```
## Requirements
* Cmake 2.8 or higher
* CPLEX 12.6 or higher
* OS X or Linux OS

## Installation on OS X and Linux
* Obtain a copy of the CPLEX software and license [here](https://www.ibm.com/analytics/cplex-optimizer)
* Specify the CPLEX installation directory in CMakeLists.txt
* Enter the following commands on the terminal
	```
	$ mkdir build
	$ cd build/
	$ cmake Unix Makefiles -DCMAKE_BUILD_TYPE=Release ..
	$ make
	```
* You should find an executable `kadaptability` in the build/ directory

## Note
* You may also compile in debug mode by changing `Release` to `Debug` in the above command.
* You may change the `TIME_LIMIT`, `MEMORY_LIMIT` and `NUM_THREADS` parameters in inc/Constants.h. If you do modify the defaults, please re-build for the changes to take effect.

## Usage
The problem (variables, objective function, constraints, uncertainty set etc) is modeled in the abstract class `KAdaptableInfo` defined in problemInfo.hpp/cpp.

To solve a new problem, you must create a new class which will inherit from this abstract class.
You must then override the following functions:

* `makeVars`
* `makeUncSet`
* `makeConsX`
* `makeConsY`
* `create`
* `clone`

As an example, have a look at problemInfo_knp.hpp/cpp. This file implements examples 4.2 and 4.3 from the above paper.

After you have implemented the above class, you can solve it by creating a `KAdaptableSolver` object and calling the function `solve_KAdaptability(K, heuristic, x)`. The solution vector is returned in `x`, and a summary of the run will be written to standard output.

As an example, please follow the `main` function defined in src/test.cpp.

