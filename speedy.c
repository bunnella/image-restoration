#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>

int curSweep = 0;

// global parameters, y'all
struct {
	double theta;   // weight on data term
	double gamma;   // microedge penalty
	double alpha;   // robustification steepness
	double tau;     // max annealing temperature
	double kappa;   // minimax quadratic/linear crossover
	int    nlevels; // number of levels in first sampling phase
} gp;

// discrete sampling interval centers
double *V = NULL; 

// trust matrix
double *k = NULL;

// R matrices for (degraded) source and MCMC candidate + dimensions
SEXP y, x;
int  R, C;

// macros for 2-D array access
#define Y(r, c) REAL(y)[(r) + (c)*R]
#define X(r, c) REAL(x)[(r) + (c)*R]
#define K(r, c) k[(r) + (c)*R]

// stlib abs() is int-only...
#define ABS(x) ((x) < 0) ? -(x) : +(x)

// runif(1)
#define U() ((double) rand() / RAND_MAX)

// nonlinearizes each neighbor's contribution to the consistency measure
#define minimax(d) ABS(d); // (ABS(d) > gp.kappa) ? d*d : gp.kappa*gp.kappa + 2*gp.kappa*(abs(d)-gp.kappa);

// computes trust values (0-1) for each pixel based on consistency measure
void initTrustMatrix() {
	if (k) free(k);
	k = malloc(sizeof(double)*R*C);

	// assign trust on edges
	double defaultEdgeTrust = 1.0;
	for (int r = 0; r < R; ++r) {
		K(r, 0)   = defaultEdgeTrust;
		K(r, C-1) = defaultEdgeTrust;
	}
	for (int c = 1; c < C-1; ++c) {
		K(0,   c) = defaultEdgeTrust;
		K(R-1, c) = defaultEdgeTrust;
	}

	// do actual calculations on the internal pixels
	for (int r = 1; r < R-1; ++r) {
		for (int c = 1; c < C-1; ++c) {
			double trust = 0;
			for (int i = -1; i <= 1; ++i) {
				for (int j = -1; j <= 1; ++j) {
					trust += minimax(Y(r, c) - Y(r+i, c+j));
				}
			}
			K(r, c) = 1 - trust/8;
		}
	}
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
	if (r < C-1) { n++; local += f(v, X(r  , c+1)); }

	return local/n + gp.theta*d(v, Y(r, c));
}

// returns the index of the largest item in (sorted) arr which val exceeds
static int bisect(double *arr, double val) {
	int i = 0;
	int j = gp.nlevels;
	int k; 

	while (j - i > 1) {
		k = (i + j) / 2;
		if (val > arr[k])
			i = k;
		else
			j = k;
	}

	return i;
}

// M-H sampler (with annealing parameter beta)
void sampleXs(int r, int c, double beta) {
	double *energies = malloc(gp.nlevels*sizeof(double));
	double sum = 0;

	for (int i = 0; i < gp.nlevels; ++i) {
		energies[i] = exp(-beta*h(r, c, V[i]));
		sum += energies[i];
	}

	double *cumprobs = malloc((gp.nlevels+1)*sizeof(double));
	cumprobs[0] = 0;
	cumprobs[gp.nlevels] = 1;
	for (int i = 1; i < gp.nlevels; ++i)
		cumprobs[i] = cumprobs[i-1] + energies[i-1]/sum;

	double original = X(r, c);
	double proposal = V[bisect(cumprobs, U())] + (2*U()-1)/(gp.nlevels+1);
	double deltaH = h(r, c, proposal) - h(r, c, original);

	// accept proposal if energy is improved or w/ M-H acceptance prob
	if (deltaH < 0 || U() < (1-K(r, c))*exp(-beta*deltaH))
		X(r, c) = proposal;

	free(energies);
	free(cumprobs);
}

// sets up chain parameters and stuff
SEXP R_setupGibbs(SEXP R_y, SEXP R_x, SEXP R_seed, SEXP R_nlevels, SEXP R_theta, SEXP R_gamma, SEXP R_alpha, SEXP R_kappa, SEXP R_tau) {

	curSweep = 1; // ready to start (over)

	// (degraded) source and MCMC starting point
	if (R_x == R_y) error("x and y refer to identical matrices in memory!");
	y = R_y;
	x = R_x;

	// image dimensions (R rows x C columns)
	R = nrows(R_y);
	C = ncols(R_y);

	// seed random stuff!
	srand((unsigned) asInteger(R_seed));

	// energy function parameters
	gp.theta = asReal(R_theta);
	gp.gamma = asReal(R_gamma);
	gp.alpha = asReal(R_alpha);
	gp.kappa = asReal(R_kappa); // NOT USED CURRENTLY (for fancier minimax)

	// annealing parameters
	gp.tau = asReal(R_tau);

	// init gray levels for discrete sampling
	gp.nlevels = asInteger(R_nlevels);
	if (V) free(V);
	V = malloc(gp.nlevels*sizeof(double));
	for (int i = 1; i <= gp.nlevels; ++i)
		V[i] = (double) i / (gp.nlevels+1); 

	initTrustMatrix();

	return R_NilValue;
}

// runs the chain for N (additional) steps
SEXP R_runGibbs(SEXP R_N) {
	if (curSweep == 0) error("unitialized chain; call setupGibbs first");

	int N = asInteger(R_N);

	for (int n = 0; n < N; ++n) {
		printf("\r%.0f%%", (double) n / N * 100);
		double invT = gp.tau*log(curSweep); //gp.tau * (1 - exp(-gp.omicron*curSweep)); // THIS IS BETA!

		// select a pixel at random; hope for the best
		for (int s = 0; s < R*C; ++s) {
			int r, c;
			do {
				r = (int) floor(U()*R);
				c = (int) floor(U()*C);
			} while (r == R || c == C); // just in case U() returns 1.0 on the nose
			sampleXs(r, c, invT);
		}

		curSweep++;
	}

	printf("\r100%%\n");
	return x;
}
