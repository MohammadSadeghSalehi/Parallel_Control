# Parallel_Control
 Using parallel CPUs to solve a control problem.
 
Project of computing the solution of a control problem using parallel CPUs. The control problem is governed by a Schrodinger-type equation. Several algorithms can be used to solve such problems, one of which is "monotonic" algorithm or "Gradient Descent". By dividing the time interval into smaller intervals, and sending these problems to different CPUs, each CPU solves the provided problem using monotonic algorithm and returns the obtained solution to the main CPU. By simply concatenating these solutions, an approximation of the solution will be achieved after several iterations.

This project has been done as a joint work with [Arian Alavizadeh](https://github.com/alavizadeharyan) in January 2021.

For more detail about theoretical aspect please check [Documentation](https://github.com/MohammadSadeghSalehi/Parallel_Control/blob/main/Documentation.pdf).


# Reference
[Monotonic Parareal Control for Quantum Systems, Y. Maday, J. Salomon, G. Turinici, SIAM Journal On Numerical Analysis (2007)](https://epubs.siam.org/doi/abs/10.1137/050647086)
