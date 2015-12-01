# image
# Gibbs Sample/ Simulated Annealing
# RAP
# Output

library(bmp)
library(compiler)
enableJIT(3)

# SET YOUR WORKING DIRECTORY TO THE REPO!!!

Rprof.out <- "Rprof.out" # change to NULL to turn off profiling

display <- function(img, caption = "") {
  image(img/2+.5, col=gray(V/2+.5), zlim=0:1, frame=F, asp=1, xaxt="n", yaxt="n", main=caption)
}

degrade <- function(original, perturb_percentage=.2) {
  static_vector <- rnorm(R*C,0,perturb_percentage/2)
  # 90% within two standard deviations
  # which is +/- perturb_percentage
  # which is perturb_percentage as percentage of dist between -1 and 1
  random_noise <- matrix(static_vector, nrow=R, ncol=C)
  # print(head(random_noise))
  added_noise <- original+random_noise
  #added_noise
  m <- matrix(mapply(function(x) min(max(x,-1), 1),added_noise),nrow=R,ncol=C)
  m
}

neighbors <- function(i) {
	m <- ((i-1) %% R) + 1
	n <- ceiling(i / R)
	return(neighborsHelp(m, n))
}

neighborsHelp <- function(m, n) {
	neighborList <- c()
	if (m > 1) {
		neighborList <- c(neighborList, convertRCtoI(m-1, n))
	}
	if (n > 1) {
		neighborList <- c(neighborList, convertRCtoI(m, n-1))
	}
	if (m < R) {
		neighborList <- c(neighborList, convertRCtoI(m+1, n))
	}
	if (n < C) {
		neighborList <- c(neighborList, convertRCtoI(m, n+1))
	}
	return(neighborList)
}

convertItoRC <- function(i) {
	m <- ((i-1) %% R) + 1
	n <- ceiling(i / R)
	return(c(m,n))
}

convertRCtoI <- function(m, n) {
	return((m + (n-1)*R))
}

# euclidian distance
d <- function(xs, ys) {
	return((xs-ys)^2)
}

f <- function(xs, xt, e) {
	return((xs-xt)^2*(1-e) + gamma*e)
}

# energy function evaluated when x_s = v (only over neighbors of s)
H <- function(s, v) {
  theta*d(v, y[s]) + sum(sapply(neighbors(s), function(t) f(v, x[t], getME(s, t))))
}

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

##############################################################################
## Gibbs Sampler!

N <- 100 # number of sweeps
V <- seq(-1, 1, length.out = 16) # set of discrete gray levels
theta <- 4 # weight on data term
gamma <- .1 # microedge penalty

# read in the test image
picture <- read.bmp("img/small_cat.bmp")
R <- ncol(picture)
C <- nrow(picture)

# transform to image()-ready orientation and degrade
ch.R <- t(picture[C:1, 1:R, 1])
ch.G <- t(picture[C:1, 1:R, 2])
ch.B <- t(picture[C:1, 1:R, 3])
values <- 0.30*ch.R + 0.59*ch.G + 0.11*ch.B
original <- values / 127.5 - 1 # scale from [0..255] -> [-1, 1]
y <- degrade(original, perturb_percentage=.2)

# plot original + degraded, leave room for MAP estimate
par(mfrow = c(1, 3), mar = c(2.6, 1, 2.6, 1))
display(original, "Original")
display(y, "Noisy data")

# data structures/functions for manipulating microedges
eEW <- matrix(0, R, C-1)
eNS <- matrix(0, R-1, C)
getME <- function(s1, s2) {
  r1 <- convertItoRC(s1)[1]; r2 <- convertItoRC(s2)[1];
  c1 <- convertItoRC(s1)[2]; c2 <- convertItoRC(s2)[2];
  if (r1 == r2 && abs(c1 - c2) == 1) eEW[r1, min(c1, c2)] else
    if (c1 == c2 && abs(r1 - r2) == 1) eNS[min(r1, r2), c1] else
      stop("invalid neighbors, stupid")
}
setME <- function(s1, s2, e) {
  r1 <- convertItoRC(s1)[1]; r2 <- convertItoRC(s2)[1];
  c1 <- convertItoRC(s1)[2]; c2 <- convertItoRC(s2)[2];
  if (r1 == r2 && abs(c1 - c2) == 1) eEW[r1, min(c1, c2)] <<- e else
    if (c1 == c2 && abs(r1 - r2) == 1) eNS[min(r1, r2), c1] <<- e else
      stop("invalid neighbors, stupid")
}

# turn on profiling (or possibly not...)
Rprof(filename = Rprof.out)

# Gibbs time ;)
x <- y # init with degraded (given) image
for (n in 1:N) {
  beta = 45 # no annealing for now - was min(exp((n - N/2)/20), 50)
  for (s in 1:(R*C)) {
    x[s] = sampleXs(s, beta)
  }
  for (s in 1:(R*C)) {
    for (t in neighbors(s)) {
      if (t > s) {
        setME(s, t, sampleME(s, t, beta))
      }
    }
  }
  cat(paste0("\r", round(100*n/N), "%\r"))
}

# turn off profiling
Rprof(filename = NULL)

display(x, "MAP estimate")

# profile summary
if (!is.null(Rprof.out)) summaryRprof(filename = Rprof.out)
