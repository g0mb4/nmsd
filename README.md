# nmsd

Numerical solutions of the mass-spring-damper model.

## Quickstart

```bash
cc nmsd.c -o nmsd -lm
./nmsd && gnuplot -p nmsd.plt
```

> [!WARNING]
> This document lacks many important details. Its purpose is to lay down the basic ideas.
> To fully understand the presented material, the suggested papers/books must be studied.

## Implementation

The choice of the "cryptic," old style is deliberate. It serves three purposes:

+ C is an old but easy-to-learn language that can be compiled even on embedded systems. It is an efficient language, meaning it is favorable for numerical methods and it does not hide the implementation details.
+ Declaring the variables at the top of the function makes the algorithms less noisy, in my opinion.
+ The cryptic variable names prepare the reader: many numerical methods have an old C/FORTRAN implementation that can be very useful to study, but they use the same "cryptic" style because of the limitations of their time.

> [!WARNING]
> The implementation of the methods is **not** general. They are capable of solving **only a second order differential equation**,
> but they can be easily modified to extend their capabilities.
  
## Explanation

The mass-spring-damper model describes a 1-degree-of-freedom mechanical system consisting of a mass, a spring and a damper.
The equation of the motion of the unforced system:
```math
m\ddot{x} + k\dot{x} + cx = 0
```
where $m$ is the **mass** $k$ is the **damping coefficient** and $c$ is the **stiffness of the spring**. To solve this equation, we have to find the $x(t)$ **displacement from the equilibrium position**.

The governing equation is a **2nd order linear homogeneous differential equation with constant coefficients**.
There are many ways to solve this equation numerically, only a few selected will be presented in this document and their implementation can be found in the accompanying *nmsd.c* file.

### Analytical Solution (`an()`)

We can solve this equation analytically, this will be the base-line for the comparison of different numerical methods.

Let's introduce $\beta = \frac{k}{m}$ and the **squared natural frequency** $\omega^2 = \frac{c}{m}$ and assume that the solution $x(t)$ is in a form of $e^{\lambda t}$.
We can rewrite the equation:
```math
\lambda^2 e^{\lambda t} + \beta \lambda e^{\lambda t} + \omega^2 e^{\lambda t} = 0
```
after regrouping we get:
```math
\left( \lambda^2 + \beta \lambda + \omega^2 \right) e^{\lambda t} = 0
```
since $e^{\lambda t}$ cannot be $0$, we can divide by it to arrive to the **characteristic equation**:
```math
\lambda^2 + \beta \lambda + \omega^2 = 0
```
which is a **quadratic equation**. To solve it we can use the **quadratic formula**:
```math
\lambda_{1,2} = \frac {-\beta \pm \sqrt{\beta^2 - 4 \omega^2}}{2}
```
Since the original equation is **linear**, the general solution:
```math
x(t) = c_1 e^{\lambda_1 t} + c_2 e^{\lambda_2 t}
```
The $c_1$ and $c_2$ constants can be determined from the **initial conditions** which are the $x_0$ **initial position** and the $v_0$ **initial velocity**:
```math
x(t=0) = x_0 = c_1 + c_2
```
```math
\dot{x}(t=0) = v_0 =  c_1 \lambda_1 + c_2 \lambda_2
```

### Explicit Eluer method (`ee1()`)

### Classical Runge-Kutta method (`rk4()`)

### Velocity Verlet method (`vv2()`)

### Dormand–Prince method (`dp54()`)

