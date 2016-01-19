# image
# Gibbs Sample/ Simulated Annealing
# RAP
# Output

library(bmp)
library(compiler)
enableJIT(3)

# SET YOUR WORKING DIRECTORY TO THE REPO!!!

Rprof.out <- "Rprof.out" # change to NULL to turn off profiling

display <- function(img, caption = "", V) {
  image(img/2+.5, col=gray(V/2+.5), zlim=0:1, frame=F, asp=C/R, xaxt="n", yaxt="n", main=caption)
}


# euclidian distance
d <- function(xs, ys) {
  (xs-ys)^2
}

f <- function(xs, xt, alpha, gamma) {
  ((((xs-xt)^2)^-alpha) + gamma^-alpha)^(-1/alpha)
  #min((xs-xt)^2, gamma)
  #atan(gamma*(xs-xt)^2)
}

# energy function evaluated when x_s = v (only over neighbors of s)
H <- function(s, v, y, alpha, gamma) {
  e <- theta*d(v, y[s])
  r <- ((s-1) %% R) + 1
  c <- ceiling(s / R)
  if (r > 1) e <- e + f(v, x[s-1], alpha, gamma)
  if (c > 1) e <- e + f(v, x[s-R], alpha, gamma)
  if (r < R) e <- e + f(v, x[s+1], alpha, gamma)
  if (c < C) e <- e + f(v, x[s+R], alpha, gamma)
  e
}

# returns a sample from the local characteristic distribution
sampleXs <- function(s, beta = 1, V, y, alpha, gamma) {
  probs <- sapply(V, function(v) exp(-beta*H(s, v, y, alpha, gamma)))
  max(probs)/sum(probs)
}

##############################################################################
## Gibbs Sampler!
grayLevelQuantization <- c(16,32,64,128)
inverseTemperatures <- c(1,3.16,10,36.1,100)
theta <- 6 # weight on data term
gamma <- .1 # microedge penalty
alpha <- 1.5 # robustification steepness

# read in the test image
picture <- read.bmp("Anne's stuff!/img/woman_blonde.bmp")
R <- ncol(picture)
C <- nrow(picture)

# transform to image()-ready orientation and degrade
ch.R <- t(picture[C:1, 1:R, 1])
ch.G <- t(picture[C:1, 1:R, 2])
ch.B <- t(picture[C:1, 1:R, 3])
values <- 0.30*ch.R + 0.59*ch.G + 0.11*ch.B
original <- values / 127.5 - 1 # scale from [0..255] -> [-1, 1]
y <- original #degrade(original, perturb_percentage=.2)
x <- y
# plot original + degraded, leave room for MAP estimate
par(mfrow = c(2, 3), mar = c(2.6, 1, 2.6, 1))
for (glq in  grayLevelQuantization) {
	V <- seq(-1, 1, length.out = glq)
	display(original, paste0(glq, " levels"), V)
	#display(y, "Noisy data")
	
	# turn on profiling (or possibly not...)
	Rprof(filename = Rprof.out)
	
	# Gibbs time ;)
	for (invT in inverseTemperatures) {
	  beta = invT # no annealing for now - was min(exp((n - N/2)/20), 50)
	  for (s in 1:(R*C)) {
	    y[s] = sampleXs(s, beta, V, original, alpha, gamma)
	    cat(paste0("\r", round(100*s/R/C), "%\r"))
	  }
	  y <- 2*y-1
	  display(y, paste0("probabilities of max\ninverse temperature: ", invT), V)
	  y <- original 
	}
}
# turn off profiling
Rprof(filename = NULL)
# profile summary
if (!is.null(Rprof.out)) summaryRprof(filename = Rprof.out)
