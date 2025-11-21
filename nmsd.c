#include <assert.h>
#include <math.h>
#include <stdio.h>

#define MAX(a,b)    (a) > (b) ? (a) : (b)
#define MIN(a,b)    (a) < (b) ? (a) : (b)

const double m = 1;     // mass, kg
const double k = 1;     // damping coefficient, Ns/m
const double c = 0.1;   // stiffness of the spring, N/m

const double x0 = 0;    // initial displacement, m
const double v0 = 1;    // initial velocity, m/s

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
    static const double B = k/(2.0 * m);
    static const double O2 = c/m;

    static const double D = B*B - O2;
    // assume real solutions
    assert(D >= 0);

    const double L1 = -B + sqrt(D);
    const double L2 = -B - sqrt(D);

    const double c2 = (v0 - x0 * L1) / (-L1 + L2);
    const double c1 = x0 - c2;

    return c1 * exp(L1 * t) + c2 * exp(L2 * t);
}

// Explicit Eluer (1st order)
double ee1(double t, double *x, double dt)
{
    double f[2];

    msd(t, x, f);

    x[0] += dt * f[0];
    x[1] += dt * f[1];
}

// Runge-Kutta (4th order)
void rk4(double t, double *x, double dt)
{
    // temporary variables
    double f[2];
    double xk[2];

    // stage derivatives
    double k10, k11;
    double k20, k21;
    double k30, k31;
    double k40, k41;

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

    // solution
    x[0] += k10/6 + k20/3 + k30/3 + k40/6;
    x[1] += k11/6 + k21/3 + k31/3 + k41/6;
}

// Velocity Verlet (2nd order)
double vv2(double t, double *x, double dt)
{
    // temporary variables
    double f[2];
    double vh;

    // calculate acceleration
    msd(t, x, f);

    // update velocity (half-step)
    vh = f[0] + (1/2.0) * f[1] * dt;
    x[1] = vh;

    // update position
    x[0] += x[1] * dt;

    // calculate new acceleration
    msd(t, x, f);

    // update velocity
    x[1] = vh +(1/2.0) * f[1] * dt;
}

// Dormand–Prince (5th order with 4th order error estimation)
void dp54(double *t, double *x, double rtol, double *dtinit)
{
    // Butcher tableau values
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

    static const double b1 = 35/384.0;
    static const double b2 = 0.0;
    static const double b3 = 500/1113.0;
    static const double b4 = 125/192.0;
    static const double b5 = -2187/6784.0;
    static const double b6 = 11/84.0;
    static const double b7 = 0.0;

    static const double b1bar = 5179/57600.0;
    static const double b2bar = 0.0;
    static const double b3bar = 7571/16695.0;
    static const double b4bar = 393/640.0;
    static const double b5bar = -92097/339200.0;
    static const double b6bar = 187/2100.0;
    static const double b7bar = 1/40.0;

    static const double c2 = 1/5.0;
    static const double c3 = 3/10.0;
    static const double c4 = 4/5.0;
    static const double c5 = 8/9.0;
    static const double c6 = 1.0;
    static const double c7 = 1.0;

    // temporary variables
    double f[2];
    double xk[2];
    double x0, x1;
    double x0bar, x1bar;

    // stage derivatives
    double k10, k11;
    double k20, k21;
    double k30, k31;
    double k40, k41;
    double k50, k51;
    double k60, k61;
    double k70, k71;

    // starting step estimation
    double d0, d1, d2;
    double dt0, dt1;
    double xf1[2];
    double f0[2], f1[2];

    // automatic step size control
    double dt = *dtinit;
    double dtopt = 0;
    const double atol = rtol / 1000.0;    // practical value
    double sc[2];
    double err;
    const double frac = 0.9;  // practical value

    // starting step estimation
    if (dt == 0) {
        msd(*t, x, f0);

        d0 = sqrt(x[0]*x[0] + x[1]*x[1]);
        d1 = sqrt(f0[0]*f0[0] + f0[1]*f0[1]);

        if (d0 < 1e-5 || d1 < 1e-5)
            dt0 = 1e-6;
        else
            dt0 = 0.01 * (d0/d1);

        xf1[0] = x[0] + dt0 * f0[0];
        xf1[1] = x[1] + dt0 * f0[1];

        msd(*t + dt0, xf1, f1);

        d2 = sqrt(
            ((f1[0]-f0[0]) * (f1[0]-f0[0]))
          + ((f1[1]-f0[1]) * (f1[1]-f0[1]))
        );

        if (MAX(d1, d2) <= 1e-15)
            dt1 = MAX(1e-6, dt0*1e-3);
        else
            dt1 = pow(0.01 / MAX(d1, d2), 1/6.0);

        dt = MIN(100*dt0, dt1);
    }

    // to be sure
    for (int tries = 0; tries < 100; ++tries) {
        // k1
        msd(*t, x, f);
        k10 = dt * f[0];
        k11 = dt * f[1];

        // k2
        xk[0] = x[0] + dt * a21 * k10;
        xk[1] = x[1] + dt * a21 * k11;

        msd(*t + c2 * dt, xk, f);
        k20 = dt * f[0];
        k21 = dt * f[1];

        // k3
        xk[0] = x[0] + dt * (a31 * k10 + a32 * k20);
        xk[1] = x[1] + dt * (a31 * k11 + a32 * k21);

        msd(*t + c3 * dt, xk, f);
        k30 = dt * f[0];
        k31 = dt * f[1];

        // k4
        xk[0] = x[0] + dt * (a41 * k10 + a42 * k20 + a43 * k30);
        xk[1] = x[1] + dt * (a41 * k11 + a42 * k21 + a43 * k31);

        msd(*t + c4 * dt, xk, f);
        k40 = dt * f[0];
        k41 = dt * f[1];

        // k5
        xk[0] = x[0] + dt * (a51 * k10 + a52 * k20 + a53 * k30 + a54 * k40);
        xk[1] = x[1] + dt * (a51 * k11 + a52 * k21 + a53 * k31 + a54 * k41);

        msd(*t + c5 * dt, xk, f);
        k50 = dt * f[0];
        k51 = dt * f[1];

        // k6
        xk[0] = x[0] + dt * (a61 * k10 + a62 * k20 + a63 * k30
                             + a64 * k40 + a65 * k50);
        xk[1] = x[1] + dt * (a61 * k11 + a62 * k21 + a63 * k31
                             + a64 * k41 + a65 * k51);

        msd(*t + c6 * dt, xk, f);
        k60 = dt * f[0];
        k61 = dt * f[1];

        // k7
        xk[0] = x[0] + dt * (a71 * k10 + a72 * k20 + a73 * k30
                             + a74 * k40 + a75 * k50 + a76 * k60);
        xk[1] = x[1] + dt * (a71 * k11 + a72 * k21 + a73 * k31
                             + a74 * k41 + a75 * k51 + a76 * k61);

        msd(*t + c7 * dt, xk, f);
        k70 = dt * f[0];
        k71 = dt * f[1];

        // fifth-order solution
        x0 = x[0] + b1 * k10 + b2 * k20 + b3 * k30 + b4 * k40
             + b5 * k50 + b6 * k60 + b7 * k70;
        x1 = x[1] + b1 * k11 + b2 * k21 + b3 * k31 + b4 * k41
             + b5 * k51 + b6 * k61 + b7 * k71;

        // fourth-order solution (error estimation)
        x0bar = x[0] + b1bar * k10 + b2bar * k20 + b3bar * k30 + b4bar * k40
                + b5bar * k50 + b6bar * k60 + b7bar * k70;
        x1bar = x[1] + b1bar * k11 + b2bar * k21 + b3bar * k31 + b4bar * k41
                + b5bar * k51 + b6bar * k61 + b7bar * k71;

        // automatic step size control
        sc[0] = atol + MAX(x[0],x0)*rtol;
        sc[1] = atol + MAX(x[1],x1)*rtol;

        err = sqrt((1/2.0) * (
              (((x0 - x0bar)/sc[0]) * ((x0 - x0bar)/sc[0]))
            + (((x1 - x1bar)/sc[1]) * ((x1 - x1bar)/sc[1]))
        ));

        dtopt = dt * frac * pow(1/err, 1/5.0);

        // accept solution
        if (err <= 1.0) {
            x[0] = x0;
            x[1] = x1;

            *dtinit = dtopt;
            *t = *t + dt;
            return;
        }

        // else run calculation with updated dt
        dt = dtopt;
    }

    assert(0 && "Too many tries");
}

int main(void)
{
    FILE *ffixed, *fdp54;
    double xee[2], xvv[2], xrk[2], xdp[2];

    double a;
    double eee, evv, erk, edp;

    double dt = 0.1;
    double tend = 10;
    double dtdp54 = 0;
    double rtol = 6e-7;

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

    for (double t = 0; t <= tend; t += dt) {
        a = an(t);
        eee = fabs(a - xee[0]);
        evv = fabs(a - xvv[0]);
        erk = fabs(a - xrk[0]);

        fprintf(ffixed, "%G; %G;", t, a);
        fprintf(ffixed, "%G; %G;", xee[0], eee);
        fprintf(ffixed, "%G; %G;", xvv[0], evv);
        fprintf(ffixed, "%G; %G\n", xrk[0], erk);

        ee1(t, xee, dt);
        vv2(t, xvv, dt);
        rk4(t, xrk, dt);
    }

    for (double t = 0; t <= tend;) {
        a = an(t);
        edp = fabs(a - xdp[0]);

        fprintf(fdp54, "%G; %G; %G; %G\n",
                        t, a, xdp[0], edp);

        dp54(&t, xdp, rtol, &dtdp54);
    }

    fclose(ffixed);
    fclose(fdp54);
}
