#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>

int curSweep = 0;

// global parameters, y'all
struct {
	int    nlevels; // number of quantization levels
	double theta;   // weight on data term
	double gamma;   // microedge penalty
	double alpha;   // robustification steepness
	double tau;     // max annealing temperature
	double omicron; // annealing steepness
	double kappa;   // minimax quadratic/linear crossover
} gp;

// gray values
double *V;

// R matrices for (degraded) source and MCMC candidate + dimensions
SEXP y, x;
int  R, C;

// trust matrix
double *k = NULL;

// macros for 2-D array access
#define Y(r, c) REAL(y)[(r) + (c)*R]
#define X(r, c) REAL(x)[(r) + (c)*R]
#define K(r, c) k[(r) + (c)*R]

// stlib abs() is int-only...
#define ABS(x) ((x) < 0) ? -(x) : +(x)

// data term
static double d(double xs, double ys) {
	return (xs-ys)*(xs-ys);
}

// local characteristic term
static double f(double xs, double xt) {
	return pow(pow(ABS(xs-xt), -2*gp.alpha) + pow(gp.gamma, -gp.alpha), -1/gp.alpha);
}

// nonlinearizes each neighbor's contribution to the consistency measure
static inline double minimax(double d) {
	return ABS(d); // (ABS(d) > gp.kappa) ? d*d : gp.kappa*gp.kappa + 2*gp.kappa*(abs(d)-gp.kappa);
}

// computes trust values (0-1) for each pixel based on consistency measure
static void initTrustMatrix() {
	if (!k) k = malloc(sizeof(double)*R*C);

	// assign trust on edges
	double defaultEdgeTrust = 0.5;
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

#define NEIGHBOR(r, c) { ++n; for (int i = 0; i < gp.nlevels; ++i) energies[i] += f(V[i], X(r, c)); }

#define U() (double) rand() / RAND_MAX // runif(1)

// Gibbs sampler (with annealing parameter beta)
static void sampleXs(int r, int c, double beta) {
	double *energies = calloc(gp.nlevels, sizeof(double));
	double sum = 0;
	int n = 0;

	if (r > 0)   NEIGHBOR(r-1, c  );
	if (c > 0)   NEIGHBOR(r  , c-1);
	if (r < R-1) NEIGHBOR(r+1, c  );
	if (r < C-1) NEIGHBOR(r  , c+1);

	double kappa = K(r, c);
	for (int i = 0; i < gp.nlevels; ++i) {
		energies[i] = energies[i]/n + kappa*gp.theta*d(V[i], Y(r, c));
		energies[i] = exp(-beta*energies[i]);
		sum += energies[i];
	}

// LOCAL MAX CODE COMMENT OUT IF NECESSARY
	int maxIndex = -1;
	double maxVal = -1;
	for (int i = 0; i < gp.nlevels; ++i) {
		if (energies[i] > maxVal) {
			maxIndex = i;
			maxVal = energies[i];
		}
	}

	X(r, c) = V[maxIndex];
	free(energies);
	return;
// END LOCAL MAX CODE

	double *cumprobs = malloc((1+gp.nlevels)*sizeof(double));
	cumprobs[0] = 0;
	cumprobs[gp.nlevels] = 1;
	for (int i = 1; i < gp.nlevels; ++i) {
		cumprobs[i] = cumprobs[i-1] + energies[i-1]/sum;
	}
	free(energies);

	double v = V[bisect(cumprobs, U())];
	X(r, c) = v;

	free(cumprobs);
}

// sets up chain parameters and stuff
SEXP R_setupGibbs(SEXP R_y, SEXP R_x, SEXP R_seed, SEXP R_V, SEXP R_theta, SEXP R_gamma, SEXP R_alpha, SEXP R_kappa, SEXP R_tau, SEXP R_omicron) {

	curSweep = 1; // ready to start (over)

	// (degraded) source and already MCMC starting point
	if (R_x == R_y) error("x and y refer to identical matrices in memory!");
	y = R_y;
	x = R_x;

	// image dimensions (R rows x C columns)
	R = nrows(R_y);
	C = ncols(R_y);

	// seed random stuff!
	srand((unsigned) asInteger(R_seed));

	// setup gray values
	if (R_V == R_NilValue) error("V isn't initialized and that should be different...");
	V = REAL(R_V);
	gp.nlevels = length(R_V);

	// energy function parameters
	gp.theta = asReal(R_theta);
	gp.gamma = asReal(R_gamma);
	gp.alpha = asReal(R_alpha);
	gp.kappa = asReal(R_kappa); // NOT USED CURRENTLY (for fancier minimax)

	// annealing parameters
	gp.tau = asReal(R_tau);
	gp.omicron = asReal(R_omicron);

	initTrustMatrix();

	return R_NilValue;
}

// runs the chain for N (additional) steps
SEXP R_runGibbs(SEXP R_N) {
	if (!curSweep) error("Gibbs sampler unitialized; call setupGibbs first");

	int N = asInteger(R_N);


	for (int n = 0; n < N; ++n) {
		printf("\r%.0f%%", (double) n / N * 100);
		double invT = gp.tau * (1 - exp(-gp.omicron*curSweep)); // THIS IS BETA!

		// select a pixel at random; hope for the best
		for (int s = 0; s < R*C; ++s) {
			int r, c;
			do {
				r = (int) floor(U()*R);
				c = (int) floor(U()*C);
			} while (r == R || c == C); // just in case U() returns 1.0 on the nose
			sampleXs(r, c, invT);
		}

		/* old code where we sweep row-by-row
		for (int r = 0; r < R; ++r)
			for (int c = 0; c < C; ++c)
				sampleXs(r, c, invT);
		*/

		curSweep++;
	}

	printf("\r100%%");
	return x;
}
