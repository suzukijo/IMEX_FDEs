# IMEX_FDEs
Implicit-explicit numerical integration scheme for fractional differential equations (FDEs)

This project contains the first and second-order IMEX time-integration approaches for stiff/nonlinear FDEs presented in:

Zhou, Y., Suzuki, J.L., Zhang, C., Zayernouri, M. "Implicit-explicit time integration of nonlinear fractional differential equations". Applied Numerical Mathematics 156, 2020.

The main features of the presented IMEX approaches are:

* First- and second-order implicit-explicit (IMEX) solver for time-integration of stiff/nonlinear FDEs equations with fractional order &alpha; in (0,1), with proven convergence and linear stability
* The methods are based on a linear multi-step fractional Adams-Moulton method (FAMM), followed by the extrapolation of the nonlinear force terms
* The singularities nearby the initial time are addressed through Lubich-like corrections
* A fast inversion scheme is employed to achieve a computational complexity of O(N log N), where N denotes the number of time-steps
* Currently, the implemented methods support the solution of systems of FDEs with single fractional orders nonlinear force terms
* The example structures involve only a few steps on domain definition, correction terms, the main solver call, and the right-hand-side definition _see src/Ex3_Case3.m_ 

Users can access the _examples/_ folder for the following working examples:

* __Ex2.m__ - Solution of a stiff system of FDEs from a fabricated multi power-law type solution
* __Ex3_Case2.m__ - Solution of a nonlinear FDE from  a fabricated polynomial solution
* __Ex3_Case3.m__ - Solution of a nonlinear FDE with a nonlinear/harmonic RHS

The implemented solvers and auxiliary functions can be found in the _src/_ folder, which contains the following files:

* __IMEX_I.m__     - First-order IMEX solver for single nonlinear FDE problems
* __IMEX_II.m__    - Second-order IMEX solver for single nonlinear FDE problems
* __IMEX_I_A.m__   - First-order IMEX solver for problems involving a nonlinear system of FDEs
* __IMEX_II_A.m__  - Second-order IMEX solver for problems involving a nonlinear system of FDEs

Feel free to contact Jorge Suzuki at suzukijo@msu.edu if you have any questions.
