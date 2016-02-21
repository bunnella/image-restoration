#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <stdlib.h>
#include <math.h>

int curSweep = -1;

// global parameters, y'all
struct {
    double theta; // weight on data term
    double gamma; // microedge penalty
    double alpha; // robustification steepness
    double tau;   // annealing constant
    double step;  // random walk perturbation size
    double kappa; // minimax quadratic/linear crossover
} gp;

// R matrices for (degraded) source and MCMC candidate + dimensions
SEXP y, x;
int  R, C;

// truzt matrix
double *z = NULL;

// macros for 2-D array access
#define X(r, c) REAL(x)[(r) + (c)*R]
#define Y(r, c) REAL(y)[(r) + (c)*R]
#define Z(r, c) z[(r) + (c)*R]

// runif(1)
#define U() ((double) rand() / RAND_MAX)

// stlib abs() is int-only...
inline double ABS(double x) { return (x < 0) ? -x : +x; }

// nonlinearizes each neighbor's contribution to the inconsistency (distrust) measure
inline double minimax(double d) {
    return (ABS(d) < gp.kappa) ? d*d : gp.kappa*gp.kappa + 2*gp.kappa*(ABS(d)-gp.kappa);
} 

// quantifies our level of trust [0, 1] in the given data pixel
double trust(int r, int c) {
    if (r == 0 || c == 0 || r == R-1 || c == C-1)
        return 1.0; // default edge trust

    double distrust = 0;
    for (int i = -1; i <= 1; ++i)
        for (int j = -1; j <= 1; ++j)
            distrust += minimax(Y(r, c) - Y(r+i, c+j));

    return 1 - distrust/(2*gp.kappa-gp.kappa*gp.kappa)/8;
}

// data term
inline double d(double xs, double ys) {
    return (xs-ys)*(xs-ys);
}

// local characteristic term
inline double f(double xs, double xt) {
    return pow(pow(ABS(xs-xt), -2*gp.alpha) + pow(gp.gamma, -gp.alpha), -1/gp.alpha);
}

// combined energy function
double h(int r, int c, double v) {
    double local = 0;
    int n = 0;

    if (r > 0)   { n++; local += f(v, X(r-1, c  )); }
    if (c > 0)   { n++; local += f(v, X(r  , c-1)); }
    if (r < R-1) { n++; local += f(v, X(r+1, c  )); }
    if (c < C-1) { n++; local += f(v, X(r  , c+1)); }

    return local/n + gp.theta*d(v, Y(r, c));
}

// MCMC proposal/update at inverse temperature beta
void sampleXs(int r, int c, double beta) {
    double original = X(r, c);
    double proposal = original + (2*(U()<.5)-1)*gp.step;

    // reflect at extremes
    if (proposal < 0) proposal =   gp.step;
    if (proposal > 1) proposal = 1-gp.step;

    // account for assymmetry at extremes (only possible to get there from one state)
    double edgeAdj = 1.0;
    if (proposal == 0 || proposal == 1) edgeAdj = 2.0;
    if (original == 0 || original == 1) edgeAdj = 0.5;

    // accept proposal if energy is improved or w/ M-H acceptance prob
    double deltaH = h(r, c, proposal) - h(r, c, original);
    if (deltaH < 0 || U() < edgeAdj*(1-Z(r, c))*exp(-beta*deltaH))
        X(r, c) = proposal;
}

// sets up chain parameters and stuff
SEXP R_setupGibbs(SEXP R_y, SEXP R_x, SEXP R_seed, SEXP R_theta, SEXP R_gamma, SEXP R_alpha, SEXP R_tau, SEXP R_step, SEXP R_kappa) {
    curSweep = 0; // ready to start (over)

    // (degraded) source and MCMC starting point
    if (R_x == R_y) error("x and y refer to identical matrices in memory!");
    y = R_y;
    x = R_x;

    // image dimensions (R rows x C columns)
    R = nrows(R_y);
    C = ncols(R_y);

    // seed random engine
    srand((unsigned) asInteger(R_seed));

    // energy function parameters
    gp.theta = asReal(R_theta);
    gp.gamma = asReal(R_gamma);
    gp.alpha = asReal(R_alpha);

    // MCMC/annealing parameters
    gp.tau   = asReal(R_tau);
    gp.step  = asReal(R_step);
    gp.kappa = asReal(R_kappa);

    // init trust matrix
    if (z) free(z);
    z = malloc(sizeof(double)*R*C);
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < C; ++c) 
            Z(r, c) = trust(r, c);

    return R_NilValue;
}

// runs the chain for N (additional) steps
SEXP R_runGibbs(SEXP R_N) {
    if (curSweep == -1) error("unitialized chain; call setupGibbs first");

    int N = asInteger(R_N);

    for (int n = 0; n < N; ++n) {
        printf("\r%.0f%%", (double) n / N * 100); // progress

        double beta = gp.tau*log(curSweep); // annealing inverse temperature

        // select a pixel at random; hope for the best
        for (int s = 0; s < R*C; ++s) {
            int r, c;
            do {
                r = (int) floor(U()*R);
                c = (int) floor(U()*C);
            } while (r == R || c == C); // just in case U() returns 1.0 on the nose

            sampleXs(r, c, beta);
        }

        curSweep++;
        R_CheckUserInterrupt();
    }

    printf("\r100%%\n");
    return x;
}

// returns the number of sweeps completed (-1 if uninitizalized)
SEXP R_numSteps() {
    return ScalarInteger(curSweep);
}
