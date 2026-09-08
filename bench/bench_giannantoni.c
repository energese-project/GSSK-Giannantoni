/* bench_giannantoni.c — do Giannantoni's claims survive measurement?
 *
 * Two claims are testable here, and they are NOT the same claim. Keeping them
 * apart is most of the value of this benchmark.
 *
 *   SPEED / RELIABILITY (parts A and B).  For a system in exponential form
 *   f(t) = e^phi(t), the incipient route evaluates a closed form: one call,
 *   independent of any step size, exact to machine precision. The traditional
 *   route integrates df/dt = phi'(t) f forward in N steps, paying O(N) and
 *   accumulating truncation error. This is a real and large gain -- but be
 *   precise about what it is evidence *for*. It is closed form beating
 *   numerical stepping, which is unsurprising on its own. The specifically
 *   Giannantoni claim is that the closed form is *available* for a far wider
 *   class of systems than TDC assumes. This benchmark measures the payoff
 *   when that holds; it cannot establish how often it holds. See
 *   docs/giannantoni_assessment.md 5.3 on the Level 1 / Level 2 split.
 *
 *   ANALYTIC DRIFT (part C).  This one is specific to the calculi and owes
 *   nothing to step size. At order n >= 2 the incipient and traditional
 *   derivative *operators* genuinely disagree: incipient gives (phi')^n,
 *   traditional gives the complete Bell polynomial B_n, and the difference
 *   psi_n does not shrink as dt -> 0 because it is not an integration error.
 *   Part C demonstrates that by finite-differencing the exact solution and
 *   watching it converge to the traditional value while the gap to the
 *   incipient value stays pinned at psi_n.
 *
 * The RK4 and Euler baselines are the standard method on the same ODE, not
 * weakened to flatter the comparison, and there is no JSON parsing or graph
 * machinery inside the timing loop. So this measures mathematics against
 * mathematics rather than one program's infrastructure against another's.
 *
 * Built against the real engine (lib/engine.o), not a reimplementation.
 */

#define _POSIX_C_SOURCE 199309L

#include "engine.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* The test system: phi(t) = c1 t + c2 t^2, so phi' = c1 + 2 c2 t is genuinely
 * time-varying. That is the "linear ODE, variable k(t)" row of the assessment
 * table -- the case where the papers claim a real gain, and the case where
 * psi_n is non-zero. A constant-coefficient phi would make both claims
 * trivially uninteresting: the two calculi coincide and RK4 is already exact
 * to its own order. */
#define C1  0.7
#define C2  0.9
#define T_END 2.0

static gia_phi PHI;

static double phi_eval(double t)  { return C1 * t + C2 * t * t; }
static double phi_prime(double t) { return C1 + 2.0 * C2 * t; }

/* Ground truth. */
static double f_exact(double t) { return exp(phi_eval(t)); }

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* ---------------- the traditional route: step it forward ---------------- */

static double rk4_to(double T, long steps) {
    double h = T / (double)steps, t = 0.0, y = 1.0; /* f(0) = e^0 = 1 */
    long   i;
    for (i = 0; i < steps; i++) {
        double k1 = phi_prime(t)             * y;
        double k2 = phi_prime(t + 0.5 * h)   * (y + 0.5 * h * k1);
        double k3 = phi_prime(t + 0.5 * h)   * (y + 0.5 * h * k2);
        double k4 = phi_prime(t + h)         * (y + h * k3);
        y += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t += h;
    }
    return y;
}

static double euler_to(double T, long steps) {
    double h = T / (double)steps, t = 0.0, y = 1.0;
    long   i;
    for (i = 0; i < steps; i++) {
        y += h * phi_prime(t) * y;
        t += h;
    }
    return y;
}

/* ---------------- the incipient route: evaluate the form ---------------- */

static double idc_to(double T) {
    /* Order 0 is the state itself. This goes through the shipped engine. */
    return gia_idc_derivative(&PHI, 0, T);
}

static double rel_err(double got, double want) {
    return fabs(got - want) / fabs(want);
}

/* Repeat until at least `budget` seconds have elapsed, so a single fast
 * evaluation is not lost in clock granularity. Returns seconds per call. */
static double time_per_call(double (*fn)(double, long), double T, long steps,
                            double budget, double *out_value) {
    double t0 = now_sec(), elapsed;
    long   reps = 0;
    volatile double sink = 0.0;
    do {
        sink += fn(T, steps);
        reps++;
        elapsed = now_sec() - t0;
    } while (elapsed < budget);
    if (out_value) *out_value = fn(T, steps);
    (void)sink;
    return elapsed / (double)reps;
}

static double idc_wrapper(double T, long steps) { (void)steps; return idc_to(T); }

/* ================================================================== */

static void part_a(void) {
    static const long STEPS[] = {10, 100, 1000, 10000, 100000, 1000000};
    const int  n = (int)(sizeof(STEPS) / sizeof(STEPS[0]));
    double     truth = f_exact(T_END);
    double     idc_t, idc_v, idc_e;
    double     best_e, best_t = 0.0;
    int        i, best_i;

    printf("\n");
    printf("A. COST OF REACHING THE STATE AT t = %.1f\n", T_END);
    printf("   task: report f(%.1f) for df/dt = phi'(t) f, f(0) = 1\n", T_END);
    printf("   exact: %.15g\n\n", truth);

    idc_t = time_per_call(idc_wrapper, T_END, 0, 0.25, &idc_v);
    idc_e = rel_err(idc_v, truth);

    printf("   %-10s %12s %14s %14s %12s\n",
           "method", "steps", "time/call", "rel error", "vs IDC");
    printf("   %-10s %12s %14s %14s %12s\n",
           "------", "-----", "---------", "---------", "------");
    printf("   %-10s %12s %11.3f us %14.3e %12s\n",
           "incipient", "1", idc_t * 1e6, idc_e, "1.0x");

    best_e = 1.0; best_i = 0;
    for (i = 0; i < n; i++) {
        double v, tt, e;
        tt = time_per_call(rk4_to, T_END, STEPS[i], 0.25, &v);
        e  = rel_err(v, truth);
        if (e < best_e) { best_e = e; best_i = i; best_t = tt; }
        printf("   %-10s %12ld %11.3f us %14.3e %11.0fx\n",
               "RK4", STEPS[i], tt * 1e6, e, tt / idc_t);
    }
    printf("\n");
    for (i = 0; i < n; i++) {
        double v, tt, e;
        tt = time_per_call(euler_to, T_END, STEPS[i], 0.25, &v);
        e  = rel_err(v, truth);
        printf("   %-10s %12ld %11.3f us %14.3e %11.0fx\n",
               "Euler", STEPS[i], tt * 1e6, e, tt / idc_t);
    }

    printf("\n   Read this as cost-to-accuracy, not raw speed. The incipient\n");
    printf("   evaluation is one closed form: its cost does not depend on a\n");
    printf("   step count and its error is roundoff. RK4 buys accuracy with\n");
    printf("   steps and never reaches roundoff before the step cost does.\n");

    printf("\n   RK4's best accuracy is %.3e at %ld steps (%.1f us), which is\n",
           best_e, STEPS[best_i], best_t * 1e6);
    printf("   %.0fx the incipient cost for %.0f decimal digits less accuracy.\n",
           best_t / idc_t,
           floor(log10(best_e / (idc_e > 0.0 ? idc_e : 1e-17))));
    if (best_i < n - 1) {
        printf("   Past %ld steps the error rises again: truncation error has\n",
               STEPS[best_i]);
        printf("   fallen below accumulated roundoff, so more work buys less\n");
        printf("   accuracy. There is a floor here, and it is above exact.\n");
    }
}

static void part_b(void) {
    static const double HORIZON[] = {1.0, 2.0, 4.0, 8.0, 16.0};
    const int    n  = (int)(sizeof(HORIZON) / sizeof(HORIZON[0]));
    const double dt = 1.0e-3;
    int          i;

    printf("\n");
    printf("B. DOES THE ERROR ACCUMULATE WITH THE HORIZON?\n");
    printf("   fixed step dt = %g; horizon grows\n\n", dt);
    printf("   %8s %10s %16s %16s\n", "T", "steps", "RK4 rel err", "IDC rel err");
    printf("   %8s %10s %16s %16s\n", "-", "-----", "-----------", "-----------");

    for (i = 0; i < n; i++) {
        double T     = HORIZON[i];
        long   steps = (long)(T / dt);
        double truth = f_exact(T);
        double r     = rel_err(rk4_to(T, steps), truth);
        double d     = rel_err(idc_to(T), truth);
        printf("   %8.1f %10ld %16.3e %16.3e\n", T, steps, r, d);
    }
    printf("\n   RK4's error is a running total: it is paid per step and never\n");
    printf("   given back. The incipient value is evaluated at T directly, so\n");
    printf("   there is no history for an error to accumulate through.\n");
}

static void part_c(void) {
    static const double H[] = {1e-1, 1e-2, 1e-3, 1e-4};
    const int    n = (int)(sizeof(H) / sizeof(H[0]));
    const double t = 1.0;
    double idc2, tdc2, psi, ef;
    int    i;

    ef   = exp(phi_eval(t));
    idc2 = gia_idc_amplitude(&PHI, 2, t) * ef;
    tdc2 = gia_tdc_amplitude(&PHI, 2, t) * ef;
    psi  = gia_drift(&PHI, 2, t) * ef;

    printf("\n");
    printf("C. THE ANALYTIC DRIFT IS NOT AN INTEGRATION ERROR\n");
    printf("   second derivative of the exact f at t = %.1f\n\n", t);
    printf("   incipient  (phi')^2 e^phi          = %.12g\n", idc2);
    printf("   traditional B_2      e^phi          = %.12g\n", tdc2);
    printf("   psi_2 = B_2 - (phi')^2, times e^phi = %.12g\n\n", psi);

    printf("   Central differences of the EXACT solution, refined:\n\n");
    printf("   %10s %18s %18s\n", "h", "|fd - traditional|", "|fd - incipient|");
    printf("   %10s %18s %18s\n", "-", "------------------", "----------------");
    for (i = 0; i < n; i++) {
        double h  = H[i];
        double fd = (f_exact(t + h) - 2.0 * f_exact(t) + f_exact(t - h)) / (h * h);
        printf("   %10.0e %18.3e %18.3e\n",
               h, fabs(fd - tdc2), fabs(fd - idc2));
    }

    printf("\n   The finite difference converges on the traditional value and\n");
    printf("   stays a fixed distance from the incipient one. So psi is not\n");
    printf("   something a better solver removes -- the two calculi disagree\n");
    printf("   about what the second derivative of this function IS. Whether\n");
    printf("   that disagreement is a defect of TDC or of IDC is a modelling\n");
    printf("   question, not one this benchmark can settle; what it settles is\n");
    printf("   that the gap is structural and survives dt -> 0.\n");
}

int main(void) {
    memset(&PHI, 0, sizeof(PHI));
    PHI.c[1]   = C1;
    PHI.c[2]   = C2;
    PHI.degree = 2;

    printf("============================================================\n");
    printf(" GIANNANTONI vs TRADITIONAL -- STRESS TEST\n");
    printf("============================================================\n");
    printf(" system   f(t) = e^phi(t),  phi(t) = %.2f t + %.2f t^2\n", C1, C2);
    printf(" ODE      df/dt = phi'(t) f,  phi'(t) = %.2f + %.2f t\n", C1, 2.0 * C2);
    printf(" note     phi is quadratic, so phi'' != 0 and psi_n != 0.\n");
    printf("          A constant-coefficient phi would make both calculi\n");
    printf("          agree and there would be nothing here to measure.\n");

    part_a();
    part_b();
    part_c();

    printf("\n============================================================\n");
    return EXIT_SUCCESS;
}
