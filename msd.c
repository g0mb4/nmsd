#include <stdio.h>

double m = 1;
double k = 0.2;
double c = 0.1;

double x0 = 0;
double v0 = 1;

// mass-spring-damper
void msd(double t, double *x, double *f)
{
    f[0] = x[1];
    f[1] = (1/m) * (-k*x[1] - c*x[0]);
}

// Explicit Eluer (1st order)
double ee1(double t, double *x, double *f, double dt)
{
    msd(t, x, f);

    x[0] += f[0] * dt;
    x[1] += f[1] * dt;
}

// Velocity Verlet (2nd order)
double vv2(double *x, double *f, double dt)
{
    // calculate acceleration
    msd(0, x, f);

    // update velocity (half-step)
    double vh = f[0] + (1/2.0) * f[1] * dt;
    x[1] = vh;

    // update position
    x[0] += x[1] * dt;

    // calculate new acceleration
    msd(0, x, f);

    // update velocity
    x[1] = vh +(1/2.0) * f[1] * dt;
}

int main(void)
{
    double xee[2], fee[2];
    double xvv[2], fvv[2];

    xee[0] = x0;
    xee[1] = v0;

    xvv[0] = x0;
    xvv[1] = v0;

    double dt = 0.1;
    for (double t = 0; t <= 10; t += dt) {
         printf("%G; %G; %G\n", t, xee[0], xvv[0]);

         ee1(t, xee, fee, dt);
         vv2(xvv, fvv, dt);
    }
}
