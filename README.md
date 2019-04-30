# Generalized Pendulum Software
 
- Describe your pendulums and their linkage.

- Get the Euler-Lagrange equations of your system.

- Try it in the simulator.

Test in simulator:
```
   $ python model.py 0 
```
Generate the equations:
```
   $ python model.py 1
```
# Description

This program uses the Lagrange equation of a two dimensional pendulum system to generate the Euler-Lagrange equation.

It also generates the Jacobian matrix of the generalized position. (Useful for external forces)

The simulator uses the basic equation of motion with a 4th order Runge-Kutta to reduce numerical error.

This software is ideal for beginners in physical simulation.

# Prerequisites

Sympy, Numpy and Matplotlib

