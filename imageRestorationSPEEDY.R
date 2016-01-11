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
  image(img/2+.5, col=gray(V/2+.5), zlim=0:1, frame=F, asp=C/R, xaxt="n", yaxt="n", main=caption)
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

# euclidian distance
d <- function(xs, ys) {
  (xs-ys)^2
}

f <- function(xs, xt) {
  #1/((xs-xt)^-2/alpha + 1/gamma)
  min((xs-xt)^2, gamma)
  #atan(gamma*(xs-xt)^2)
}

# energy function evaluated when x_s = v (only over neighbors of s)
H <- function(s, v) {
  e <- theta*d(v, y[s])
  r <- ((s-1) %% R) + 1
  c <- ceiling(s / R)
  if (r > 1) e <- e + f(v, x[s-1])
  if (c > 1) e <- e + f(v, x[s-R])
  if (r < R) e <- e + f(v, x[s+1])
  if (c < C) e <- e + f(v, x[s+R])
  e
}

# returns a sample from the local characteristic distribution
sampleXs <- function(s, beta = 1) {
  probs <- sapply(V, function(v) exp(-beta*H(s, v)))
  sample(V, 1, prob = probs)
}

##############################################################################
## Gibbs Sampler!

N <- 100 # number of sweeps
V <- seq(-1, 1, length.out = 32) # set of discrete gray levels
theta <- 4 # weight on data term
gamma <- .4 # microedge penalty
alpha <- 5 # robustification steepness

# read in the test image
picture <- read.bmp("img/papercat69.bmp")
R <- ncol(picture)
C <- nrow(picture)

# transform to image()-ready orientation and degrade
ch.R <- t(picture[C:1, 1:R, 1])
ch.G <- t(picture[C:1, 1:R, 2])
ch.B <- t(picture[C:1, 1:R, 3])
values <- 0.30*ch.R + 0.59*ch.G + 0.11*ch.B
original <- values / 127.5 - 1 # scale from [0..255] -> [-1, 1]
y <- original #degrade(original, perturb_percentage=.2)

# plot original + degraded, leave room for MAP estimate
par(mfrow = c(1, 3), mar = c(2.6, 1, 2.6, 1))
display(original, "Original")
display(y, "Noisy data")

# turn on profiling (or possibly not...)
Rprof(filename = Rprof.out)

# Gibbs time ;)
x <- y # init with degraded (given) image
for (n in 1:N) {
  beta = 40 # no annealing for now - was min(exp((n - N/2)/20), 50)
  for (s in 1:(R*C)) {
    x[s] = sampleXs(s, beta)
  }
  cat(paste0("\r", round(100*n/N), "%\r"))
}

# turn off profiling
Rprof(filename = NULL)

display(x, "MAP estimate")

# profile summary
if (!is.null(Rprof.out)) summaryRprof(filename = Rprof.out)
