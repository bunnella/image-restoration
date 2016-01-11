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

// microedge data structures
int **eEW, **eNS;

static void init_microedges() {
	eEW = malloc(sizeof(int*) * R);
	eNS = malloc(sizeof(int*) * R);

	for (int i = 0; i < R; i++) {
		eEW[i] = malloc(sizeof(int) * C);
		eNS[i] = malloc(sizeof(int) * C);
		for (int j = 0; j < C; j++) {
			eEW[i][j] = 0;
			eNS[i][j] = 0;
		}
	}
}

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

int getME(int r1, int c1, int r2, int c2) {
	#ifndef SPEEDY
	assert((r1 == r2 && abs(c1-c2) == 1) || (c1 == c2 && abs(r1-r2) == 1));
	#endif

	if (r1 == r2) {
		if (c1 < c2) return eEW[r1][c1];
		else return eEW[r1][c2];
	} else {
		if (r1 < r2) return eNS[r1][c1];
		else return eNS[r2][c1];
	}
}

void setME(int r1, int c1, int r2, int c2, int e) {
	#ifndef SPEEDY
	assert((r1 == r2 && abs(c1-c2) == 1) || (c1 == c2 && abs(r1-r2) == 1));
	#endif

	if (r1 == r2) {
		if (c1 < c2) eEW[r1][c1] = e;
		else eEW[r1][c2] = e;
	} else {
		if (r1 < r2) eNS[r1][c1] = e;
		else eNS[r2][c1] = e;
	}
}

float d(float xs, float ys) {
	return (xs-ys)*(xs-ys);
}

float f(float xs, float xt, int e) {
	return e ? r_gamma : (xs-xt)*(xs-xt);
}

#define ACCESS(mat, r, c) REAL(mat)[r + c*R]

float H(int r, int c, float v) {
	float energy = r_theta * d(v, ACCESS(y, r, c));

	if (r > 0)   energy += f(v, ACCESS(x, r-1, c  ), getME(r, c, r-1, c  ));
	if (c > 0)   energy += f(v, ACCESS(x, r,   c-1), getME(r, c, r  , c-1));
	if (r < R-1) energy += f(v, ACCESS(x, r+1, c  ), getME(r, c, r+1, c  ));
	if (c < C-1) energy += f(v, ACCESS(x, r,   c+1), getME(r, c, r  , c+1));

	return energy;
}
