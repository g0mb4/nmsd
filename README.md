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
There are many ways to solve this equation numerically, only a few will be presented in this document and their implementation can be found in the accompanying *nmsd.c* file.

### Analytical Solution (`an()`)

We can solve this equation analytically, this will be the base-line for the comparison of different numerical methods.

Let's introduce the **damping constant** $\beta = \frac{k}{2m}$ and the **natural frequency** $\omega = \sqrt{\frac{c}{m}}$ and assume that the solution $x(t)$ is in a form of $e^{\lambda t}$.
We can rewrite the equation:
```math
\lambda^2 e^{\lambda t} + 2 \beta \lambda e^{\lambda t} + \omega^2 e^{\lambda t} = 0
```
after regrouping we get:
```math
\left( \lambda^2 + 2 \beta \lambda + \omega^2 \right) e^{\lambda t} = 0
```
since $e^{\lambda t}$ cannot be $0$, we can divide by it to arrive at the **characteristic equation**:
```math
\lambda^2 + 2 \beta \lambda + \omega^2 = 0
```
which is a **quadratic equation**. To solve it we can use the **quadratic formula**:
```math
\lambda_{1,2} = -\beta \pm \sqrt{\beta^2 - \omega^2}
```
Since the original equation is **linear**, the general solution is:
```math
x(t) = c_1 e^{\lambda_1 t} + c_2 e^{\lambda_2 t}
```
The $c_1$ and $c_2$ constants can be determined from the **initial conditions** which are the **initial position** $x_0$  and the **initial velocity** $v_0$:
```math
x(t=0) = x_0 = c_1 + c_2
```
```math
\dot{x}(t=0) = v_0 =  c_1 \lambda_1 + c_2 \lambda_2
```

> [!NOTE]
> Suggested material:
> 
> Boyce W.E., DiPrima R.C., Meade D.B.: Elementary Differential Equations and Boundary Value Problems,
> ISBN 978-1-119-37792-4
>  + Chapter 3: Second-Order Linear Differential Equations
>
> Thomson W.T.: Theory of Vibration with Applications, ISBN 0-7487-4380-4    
>  + Chapter 2.6: Viscously Damped Free Vibration

### Explicit Euler method (`ee1()`)

The Euler method (also called the forward Euler method) is a numerical procedure for solving ordinary differential equations with a given initial value.

Let's recall the defintion of the **derivative**:
```math
\frac{\mathrm{d}x}{\mathrm{d}t} = \lim_{\Delta t \to 0} \frac{x(t+\Delta t) - x(t)} {\Delta t}
```
assuming $\Delta t$ is sufficiently **small**:
```math
\frac{\mathrm{d}x}{\mathrm{d}t} \approx \frac{x(t+\Delta t) - x(t)}{\Delta t}
```
after some reareagement, we get:
```math
x(t+\Delta t) \approx x(t) + \Delta t \frac{\mathrm{d}x}{\mathrm{d}t}
```
where $\Delta t$ is the **step size** or **time step**.

Since $\Delta t$ is **not** infinitesimally small, there is always an error.

The $T_{\Delta t}$ **the local truncation error** is:

```math
T_{\Delta t} = x(t+\Delta t) - x(t) - \Delta t \frac{\mathrm{d}x}{\mathrm{d}t}
```
Using **Taylor’s theorem** to estimate $T_{\Delta t}$ we get:
```math
T_{\Delta t} = \frac {\Delta t^2}{2} \frac{\mathrm{d}^2x}{\mathrm{d}t^2}
```
We say that the local truncation error of Euler’s method is of order $\Delta t^2$, which we write as $O(\Delta t^2)$.

Analizing the overall effect of truncation error we arrive to that the **global truncation error** of Euler’s method is of order 
$\Delta t$, which we write as $O(\Delta t)$. Meaning the Euler method is a **first order** numerical method.

The Eluer method can solve the following problem:

```math
\frac{\mathrm{d}x}{\mathrm{d}t} = f(t, x), x(t_0) = x_0
```
which is a **first order** ordinary differential equation with a given initial value, but we want to solve a **second order** ordinary differential equation.
We can transform our original equation to a **first order system of differential equations** as follows:

We introduce the **state vector**:
```math
\underline{x} =
\begin{bmatrix}
x_0\\
x_1
\end{bmatrix} =
\begin{bmatrix}
x\\
\dot{x}
\end{bmatrix}
```

meaning the drivative of this vector becames:

```math
\underline{\dot{x}} =
\begin{bmatrix}
\dot{x}\\
\ddot{x}
\end{bmatrix}
```

so the original equation can be written as:

```math
\underline{\dot{x}} = f(t, \underline{x})
```

```math
\begin{bmatrix}
\dot{x}\\
\ddot{x}
\end{bmatrix}
=
\begin{bmatrix}
x_1\\
\frac{1}{m} (-kx_1 - c x_0)
\end{bmatrix}
```
This is the `msd()` function in the source code.

This transformation allow us to use the Euler method to solve our problem by aplling it to **each** equation in the system given the
```math
\begin{bmatrix}
x(t_0)\\
\dot{x}(t_0)
\end{bmatrix}
=
\begin{bmatrix}
x_0\\
v_0
\end{bmatrix}
```
initial conditions.

This method is often used in simpler video games as a base of the physics engine.

### Classical Runge-Kutta method (`rk4()`)

### Velocity Verlet method (`vv2()`)

### Dormand–Prince method (`dp54()`)

