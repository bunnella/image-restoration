#include <R.h>
#include <Rinternals.h>

#include <stdlib.h>
#include <assert.h>

SEXP setGlobals(SEXP a, SEXP b) {
	SEXP result = PROTECT(allocVector(REALSXP, 2));
	REAL(result)[0] = asReal(a) + asReal(b) + 500;
	REAL(result)[1] = -24;
	UNPROTECT(1);

	return result;
}

// R matrices for (degraded) source and MCMC guess
SEXP y, x;

// dimensions
int R, C;

// constants
float r_gamma, r_theta, r_beta;

void set_images(SEXP R_x, SEXP R_y) {
	x = R_x;
	y = R_y;
}

void set_dimensions(SEXP R_R, SEXP R_C) {
	R = REAL(R_R)[0];
	C = REAL(R_C)[0];
	init_microedges();
}

void set_constants(SEXP R_beta, SEXP R_gamma, SEXP R_theta) {
	r_beta = REAL(R_beta)[0];
	r_gamma = REAL(R_gamma)[0];
	r_theta = REAL(R_theta)[0];
}

static float d(float xs, float ys) {
	return (xs-ys)*(xs-ys);
}

static float f(float xs, float xt) {
	return min((xs-xt)*(xs-xt), r_gamma);
}

#define ACCESS(mat, r, c) REAL(mat)[r + c*R]

static float H(int r, int c, float v) {
	float energy = r_theta * d(v, ACCESS(y, r, c));

	if (r > 0)   energy += f(v, ACCESS(x, r-1, c  ));
	if (c > 0)   energy += f(v, ACCESS(x, r,   c-1));
	if (r < R-1) energy += f(v, ACCESS(x, r+1, c  ));
	if (c < C-1) energy += f(v, ACCESS(x, r,   c+1));

	return energy;
}

float sampleXs <- function(int r, int c, float beta) {
	float* energies = float[nlevels];
}



/*

# returns a sample from the local characteristic distribution
sampleXs <- function(s, beta = 1) {
	probs <- sapply(V, function(v) exp(-beta*H(s, v)))
	sample(V, 1, prob = probs)
}

# returns a sampled microedge parity between pixels s & t
sampleME <- function(s, t, beta = 1) {
  probs <- sapply(0:1, function(e) exp(-beta*f(x[s], x[t], e)))
  sample(0:1, 1, prob = probs) # TODO: try without beta?
}

*/
