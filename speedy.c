#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>

int n = 0; // current sweep

// global parameters, y'all
struct {
	int   nlevels; // number of quantization levels
	double theta;   // weight on data term
	double gamma;   // microedge penalty
	double alpha;   // robustification steepness
	double tau;     // max annealing temperature
	double omicron; // annealing steepness
} gp;

// gray values
double *V;

// R matrices for (degraded) source and MCMC "guess" + dimensions
SEXP y, x;
int  R, C;

// macros for 2-D array access
#define Y(r, c) REAL(y)[(r) + (c)*R]
#define X(r, c) REAL(x)[(r) + (c)*R]

// data term
static double d(double xs, double ys) {
	return (xs-ys)*(xs-ys);
}

// local characteristic term
static double f(double xs, double xt) {
	return pow(pow((xs-xt)*(xs-xt), -gp.alpha) + pow(gp.gamma, -gp.alpha), -1/gp.alpha);
}

// combined energy function
static double H(int r, int c, double v) {
	double energy = gp.theta * d(v, Y(r, c));

	if (r > 0)   energy += f(v, X(r-1, c  ));
	if (c > 0)   energy += f(v, X(r,   c-1));
	if (r < R-1) energy += f(v, X(r-1, c  ));
	if (c < C-1) energy += f(v, X(r,   c+1));

	return energy;
}

// returns the index of the largest item in cum which val exceeds
static int bisect(double *cum, double val) {
	int i = 1;
	while (val > cum[i] && i < gp.nlevels)
		i++;
	return i-1;
}

// Gibbs sampler (with annealing parameter)
static void sampleXs(int r, int c, double beta) {
	double *energies = malloc(gp.nlevels*sizeof(double));
	double sum = 0;

	for (int i = 0; i < gp.nlevels; ++i) {
		energies[i] = exp(-beta*H(r, c, V[i]));
		sum += energies[i];		
	}

	double *cumprobs = malloc(gp.nlevels*sizeof(double));
	cumprobs[0] = 0;
	for (int i = 1; i < gp.nlevels; ++i) {
		cumprobs[i] = cumprobs[i-1] + energies[i-1]/sum;
	}
	free(energies);

	double u = (double) rand() / RAND_MAX;
	double v = V[bisect(cumprobs, u)];
	X(r, c) = v;

	free(cumprobs);
}

// sets up chain parameters and stuff
SEXP R_setupGibbs(SEXP R_y, SEXP R_x, SEXP R_V, SEXP R_theta, SEXP R_gamma, SEXP R_alpha, SEXP R_tau, SEXP R_omicron) {
	// TODO: maybe seed better?
	srand(1);

	n = 1; // ready to start (over)

	// (degraded) source and already MCMC starting point
	if (R_x == R_y) error("x and y refer to identical matrices in memory!");
	y = R_y;
	x = R_x;

	// image dimensions (R rows x C columns)
	R = nrows(R_y);
	C = ncols(R_y);

	// setup gray values
	V = REAL(R_V);
	gp.nlevels = length(R_V);

	// energy function parameters
	gp.theta = asReal(R_theta);
	gp.gamma = asReal(R_gamma);
	gp.alpha = asReal(R_alpha);

	// annealing parameters
	gp.tau = asReal(R_tau);
	gp.omicron = asReal(R_omicron);

	return R_NilValue;
}

// runs the chain for N additional steps
SEXP R_runGibbs(SEXP R_N) {
	if (!n) error("Gibbs sampler unitialized; call setupGibbs first");
	
	int N = n + asInteger(R_N);

	for (; n < N; ++n) {
		double invT = gp.tau * (1 - exp(-gp.omicron*n/N)); // THIS IS BETA!
		for (int r = 0; r < R; ++r)
			for (int c = 0; c < C; ++c)
				sampleXs(r, c, invT);
	}

	return x;
}
