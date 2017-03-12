ICQM
====

ICQM is a MATLAB script for approximating the solution to the integer convex quadratic minimization problem.
For more information, see our [companion paper](https://stanford.edu/~boyd/papers/int_least_squares.html).

ICQM requires [CVX](http://cvxr.com/cvx/). Our method is based on a semidefinite programming (SDP) relaxation, and uses the MOSEK solver if possible.

The main routine can be invoked with:
```
[lb, ub, xhat] = ICQM(P, q, r, K)
```
where the first three arguments ``P``, ``q``, and ``r`` describe the convex quadratic function to be minimized on the integer lattice. The fourth argument ``K`` is an optional parameter controlling the number of sample points used in the randomized algorithm.

The output parameters ``lb`` and ``ub`` each gives a lower and upper bound on the optimal value of the problem. The other output ``xhat`` is an integer point that attains the objective value of ``ub``.

See the ``ICQM_example.m`` file for an example.
