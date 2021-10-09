# Quadratic-Assignment-Problem-Solver

The folder 

> Quadratic-Assignment-Problem-Solver/QuadraticAssignmentProblem/

contains the implementation of a Fourier space algorithm for solving instances of the Quadratic Assignment Problem. The algorithm is based on the paper here:

http://people.cs.uchicago.edu/~risi/papers/KondorSODA10.pdf .

My implementation uses Maslen's Fast Fourier transform found in this paper: 

https://www.ams.org/journals/mcom/1998-67-223/S0025-5718-98-00964-8/S0025-5718-98-00964-8.pdf .

The time complexity is a factor-of-n faster than using Clausen's FFT.

The folder 

> Quadratic-Assignment-Problem-Solver/QAPSolverTests/

contains tests for the first project. Both require a C++20 compatible compiler.

## How to use 

An instance of ```QAPSolver``` is used to solve many QAP instances. (Why did I do it this way? Because ```QAPSolver``` caches constants that are reused across QAP instances.) 

Step 1: Initialize matrices of a matrix class that satisfies ```MatrixConcept```.
Note: ```SquareMatrix``` and ```Matrix```, included in this library, satisfy the ```MatrixConcept```.

```
QAPSolver::SquareMatrix mat1(3);

mat1(0, 0) = 0;
mat1(0, 1) = 29;
mat1(0, 2) = -13;
mat1(1, 0) = 2;
mat1(1, 1) = 0;
mat1(1, 2) = 8;
mat1(2, 0) = 17;
mat1(2, 1) = -12;
mat1(2, 2) = 0;
```

Step 2: Initialize a QAPSolver instance.

```
QAPSolver::QAPSolver qap_solver;
```

Step 3: Pass an instance of the QAP to ```qap_solver```.
Note: ```perm``` is a ```Permutation``` defined in ```Permutation.cpp```.
```total``` is a ```double```.

```
auto [perm, total] = qap_solver(mat1, mat2);
```

Step 4: Convert ```perm``` to a ```vector``` using member function ```perm_to_vec```.

```
std::vector<int> vec = perm.perm_to_vec();
```

See examples in ```qapsolvertests.cpp```.


## References
<a id="1">[1]</a> 
R. Kondor. *A Fourier space algorithm for solving quadratic assignment problems*. In SODA, 2010.

<a id="2">[2]</a> 
David K. Maslen. *The efficient computation of Fourier transforms on the symmetric group. Mathematics of Computation*, 67(223):1121–1147, 1998

