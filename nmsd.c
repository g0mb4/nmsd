#include <assert.h>
#include <math.h>
#include <stdio.h>

const double m = 1;
const double k = 1;
const double c = 0.1;

const double x0 = 0;
const double v0 = 1;

// mass-spring-damper
void msd(double t, double *x, double *f)
{
    // mx'' + kx' + cx = 0
    //
    // X' = f(t, X)
    //
    // X1 = x
    // X2 = x'
    //
    // X1' = X2
    // X2' = (1/m) (-k * X2 - c * X1)

    f[0] = x[1];
    f[1] = (1/m) * (-k*x[1] - c*x[0]);
}

// Analitic solution
double an(double t)
{
    // mx'' + kx' + cx = 0
    //
    // Z = k/m, O2 = c/m
    //
    // x'' + Zx' + O2x = 0
    // L^2 exp(Lt) + Z L exp(L t) + O2 exp(L t) = 0
    //
    // characteristic equ:
    // L^2 + Z L + O2 = 0
    //
    // L1,2 = (-Z +- sqrt(Z^2-4 O2)) / 2
    //
    // x(t) = c1 exp(L1 t) + c2 exp(L2 t)
    // x0 = c1 + c2
    // v0 = c1 L1 + c2 L2

    static const double Z = k/m;
    static const double O2 = c/m;

    static const double disc = Z*Z - 4 * O2;
    assert(disc >= 0);

    const double L1 = (-Z + sqrt(disc))/2;
    const double L2 = (-Z - sqrt(disc))/2;

    const double c2 = (v0 - x0 * L1) / (-L1 + L2);
    const double c1 = x0 - c2;

    return c1 * exp(L1 * t) + c2 * exp(L2 * t);
}

// Explicit Eluer (1st order)
double ee1(double t, double *x, double dt)
{
    double f[2];

    msd(t, x, f);

    x[0] += f[0] * dt;
    x[1] += f[1] * dt;
}

// Velocity Verlet (2nd order)
double vv2(double t, double *x, double dt)
{
    double f[2];

    // calculate acceleration
    msd(t, x, f);

    // update velocity (half-step)
    double vh = f[0] + (1/2.0) * f[1] * dt;
    x[1] = vh;

    // update position
    x[0] += x[1] * dt;

    // calculate new acceleration
    msd(t, x, f);

    // update velocity
    x[1] = vh +(1/2.0) * f[1] * dt;
}

int main(void)
{
    FILE *ffixed;
    double xee[2], xvv[2];

    xee[0] = x0;
    xee[1] = v0;

    xvv[0] = x0;
    xvv[1] = v0;

    ffixed = fopen("fixed.csv", "w");
    assert(ffixed);

    double dt = 0.1;
    for (double t = 0; t <= 10; t += dt) {
        double a = an(t);
        double eee = fabs(a - xee[0]);
        double evv = fabs(a - xvv[0]);

        fprintf(ffixed, "%G; %G;", t, a);
        fprintf(ffixed, "%G; %G;", xee[0], eee);
        fprintf(ffixed, "%G; %G\n", xvv[0], evv);

        ee1(t, xee, dt);
        vv2(t, xvv, dt);
    }

    fclose(ffixed);
}
