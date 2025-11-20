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

// Runge-Kutta (4th order)
void rk4(double t, double *x, double dt)
{
    double f[2];

    double k10, k11;
    double k20, k21;
    double k30, k31;
    double k40, k41;

    double xk[2];

    // k1
    msd(t, x, f);
    k10 = dt * f[0];
    k11 = dt * f[1];

    // k2
    xk[0] = x[0] + dt*k10/2;
    xk[1] = x[1] + dt*k11/2;

    msd(t + dt/2, xk, f);
    k20 = dt * f[0];
    k21 = dt * f[1];

    // k3
    xk[0] = x[0] + dt*k20/2;
    xk[1] = x[1] + dt*k21/2;

    msd(t + dt/2, xk, f);
    k30 = dt * f[0];
    k31 = dt * f[1];

    // k4
    xk[0] = x[0] + dt*k30;
    xk[1] = x[1] + dt*k31;

    msd(t + dt, xk, f);
    k40 = dt * f[0];
    k41 = dt * f[1];

    x[0] += k10/6 + k20/3 + k30/3 + k40/6;
    x[1] += k11/6 + k21/3 + k31/3 + k41/6;
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

// Dormand–Prince (5th order with 4th order error estimation)
void dp54(double *t_in, double *x, double rtol, double *dt_init)
{
    static const double a21 = 1/5.0;

    static const double a31 = 3/40.0;
    static const double a32 = 9/40.0;

    static const double a41 = 44/55.0;
    static const double a42 = -56/15.0;
    static const double a43 = 32/9.0;

    static const double a51 = 19372/6561.0;
    static const double a52 = -25360/2187.0;
    static const double a53 = 64448/6561.0;
    static const double a54 = -212/729.0;

    static const double a61 = 9017/3186.0;
    static const double a62 = -355/33.0;
    static const double a63 = 46732/5247.0;
    static const double a64 = 49/176.0;
    static const double a65 = -5103/18656.0;

    static const double a71 = 35/384.0;
    static const double a72 = 0.0;
    static const double a73 = 500/1113.0;
    static const double a74 = 125/192.0;
    static const double a75 = -2187/6784.0;
    static const double a76 = 11/84.0;

    static const double b11 = 35/384.0;
    static const double b12 = 0.0;
    static const double b13 = 500/1113.0;
    static const double b14 = 125/192.0;
    static const double b15 = -2187/6784.0;
    static const double b16 = 11/84.0;
    static const double b17 = 0.0;

    static const double b21 = 5179/57600.0;
    static const double b22 = 0.0;
    static const double b23 = 7571/16695.0;
    static const double b24 = 393/640.0;
    static const double b25 = -92097/339200.0;
    static const double b26 = 187/2100.0;
    static const double b27 = 1/40.0;

    static const double c2 = 1/5.0;
    static const double c3 = 3/10.0;
    static const double c4 = 4/5.0;
    static const double c5 = 8/9.0;
    static const double c6 = 1.0;
    static const double c7 = 1.0;

    double k10, k11;
    double k20, k21;
    double k30, k31;
    double k40, k41;
    double k50, k51;
    double k60, k61;
    double k70, k71;

    double f[2];
    double xk[2];

    double t = *t_in;
    double dt = *dt_init;
    double dt_new = 0;

    // TODO: initial step estimation
    if (dt == 0) {
        dt = 0.1;
    }

    for (int tries = 0; tries < 100; ++tries) {
        // k1
        msd(t, x, f);
        k10 = dt * f[0];
        k11 = dt * f[1];

        // k2
        xk[0] = x[0] + dt * a21 * k10;
        xk[1] = x[1] + dt * a21 * k11;

        msd(t + c2 * dt, xk, f);
        k20 = dt * f[0];
        k21 = dt * f[1];

        // k3
        xk[0] = x[0] + dt * (a31 * k10 + a32 * k20);
        xk[1] = x[1] + dt * (a31 * k11 + a32 * k21);

        msd(t + c3 * dt, xk, f);
        k30 = dt * f[0];
        k31 = dt * f[1];

        // k4
        xk[0] = x[0] + dt * (a41 * k10 + a42 * k20 + a43 * k30);
        xk[1] = x[1] + dt * (a41 * k11 + a42 * k21 + a43 * k31);

        msd(t + c4 * dt, xk, f);
        k40 = dt * f[0];
        k41 = dt * f[1];

        // k5
        xk[0] = x[0] + dt * (a51 * k10 + a52 * k20 + a53 * k30 + a54 * k40);
        xk[1] = x[1] + dt * (a51 * k11 + a52 * k21 + a53 * k31 + a54 * k41);

        msd(t + c5 * dt, xk, f);
        k50 = dt * f[0];
        k51 = dt * f[1];

        // k6
        xk[0] = x[0] + dt * (a61 * k10 + a62 * k20 + a63 * k30 + a64 * k40 + a65 * k50);
        xk[1] = x[1] + dt * (a61 * k11 + a62 * k21 + a63 * k31 + a64 * k41 + a65 * k51);

        msd(t + c6 * dt, xk, f);
        k60 = dt * f[0];
        k61 = dt * f[1];

        // k7
        xk[0] = x[0] + dt * (a71 * k10 + a72 * k20 + a73 * k30 + a74 * k40 + a75 * k50 + a76 * k60);
        xk[1] = x[1] + dt * (a71 * k11 + a72 * k21 + a73 * k31 + a74 * k41 + a75 * k51 + a76 * k61);

        msd(t + c7 * dt, xk, f);
        k70 = dt * f[0];
        k71 = dt * f[1];

        // fifth-order solution
        double x05 = x[0] + b11 * k10 + b12 * k20 + b13 * k30 + b14 * k40 + b15 * k50 + b16 * k60 + b17 * k70;
        double x15 = x[1] + b11 * k11 + b12 * k21 + b13 * k31 + b14 * k41 + b15 * k51 + b16 * k61 + b17 * k71;

        // fourth-order solution (error estimation)
        double x04 = x[0] + b21 * k10 + b22 * k20 + b23 * k30 + b24 * k40 + b25 * k50 + b26 * k60 + b27 * k70;
        double x14 = x[1] + b21 * k11 + b22 * k21 + b23 * k31 + b24 * k41 + b25 * k51 + b26 * k61 + b27 * k71;

        double E0 = fabs(x05 - x04);
        double E1 = fabs(x15 - x14);

        double atol = rtol / 1000.0;

        double tol0 = atol + rtol * fabs(x05);
        double tol1 = atol + rtol * fabs(x15);

        // root-mean-square (RMS)
        double E = sqrt(1/2.0 * ((E0/tol0)*(E0/tol0) + (E1/tol1)*(E1/tol1)));
        double tol = sqrt(1/2.0 * (tol0*tol0 + tol1*tol1));

        double S = 0.8;

        // recalculate step size
        // TODO: calculate next step size
        // dt_new = dt * S * pow(E/tol, 1/5.0);
        dt_new = E <= 1.0 ? dt * 2.0 : dt / 2.0;

        // accept solution
        if (E <= 1.0) {
            x[0] = x05;
            x[1] = x15;

            *dt_init = dt_new;
            *t_in = t + dt;
            return;
        }

        // else run calculation with updated dt
        dt = dt_new;
    }

    assert(0 && "Too many tries");
}

int main(void)
{
    FILE *ffixed, *fdp54;
    double xee[2], xvv[2], xrk[2], xdp[2];

    xee[0] = x0;
    xee[1] = v0;

    xvv[0] = x0;
    xvv[1] = v0;

    xrk[0] = x0;
    xrk[1] = v0;

    xdp[0] = x0;
    xdp[1] = v0;

    ffixed = fopen("fixed.csv", "w");
    assert(ffixed);

    fdp54 = fopen("dp54.csv", "w");
    assert(ffixed);

    double dt = 0.1;
    double tend = 10;
    for (double t = 0; t <= tend; t += dt) {
        double a = an(t);
        double eee = fabs(a - xee[0]);
        double evv = fabs(a - xvv[0]);
        double erk = fabs(a - xrk[0]);

        fprintf(ffixed, "%G; %G;", t, a);
        fprintf(ffixed, "%G; %G;", xee[0], eee);
        fprintf(ffixed, "%G; %G;", xvv[0], evv);
        fprintf(ffixed, "%G; %G\n", xrk[0], erk);

        ee1(t, xee, dt);
        vv2(t, xvv, dt);
        rk4(t, xrk, dt);
    }

    double dtdp54 = 0;
    for (double t = 0; t <= tend;) {
        double a = an(t);
        double edp = fabs(a - xdp[0]);

        fprintf(fdp54, "%G; %G; %G; %G\n",
                        t, a, xdp[0], edp);

        dp54(&t, xdp, 1e-6, &dtdp54);
    }

    fclose(ffixed);
    fclose(fdp54);
}
