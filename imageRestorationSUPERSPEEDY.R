library(png)
library(jpeg)

# SET YOUR WORKING DIRECTORY TO THE REPO!!!
# setwd("K:/math400-10-w16/common/image-restoration")
# setwd("~/Carleton/MATH-COMPS")

##############################################################################
## Function defs

q <- .5

display <- function(img, caption = "") {
  image(t(img[R:1, 1:C]), col=gray(0:255/255), zlim=0:1, frame=F, asp=R/C, xaxt="n", yaxt="n", main=caption)
}

mse <- function(ystar, xstar) {
  mean((ystar-xstar)^2)
}

extractGray <- function(pic) {
  ch.R <- pic[ , , 1]
  ch.G <- pic[ , , 2]
  ch.B <- pic[ , , 3]
  0.30*ch.R + 0.59*ch.G + 0.11*ch.B
}

setupGibbs <- function(
  y, x,
  seed  = 0,
  theta = .13,
  gamma = .005,
  alpha = 1,
  tau   = 100,
  step  = 1/256,
  kappa = 0.5,
  ...) .Call("R_setupGibbs", y, x, seed, theta, gamma, alpha, tau, step, kappa)

runGibbs <- function(N) .Call("R_runGibbs", N)

numSteps <- function() {
  cat(sprintf("Completed %d chain steps.\n", .Call("R_numSteps")))
}

##############################################################################
## Setup

# read in the test image
original <- readPNG("img/milk.png")
R <- nrow(original)
C <- ncol(original)

# degrade away!
ystar <- extractGray(original)
writeJPEG(ystar, target = "tmp.jpg", quality = q)
y <- readJPEG("tmp.jpg")

##############################################################################
## Gibbs sampler!

dyn.load("speedy.dll")

x <- y[] # copy

setupGibbs(y, x, theta = .13, gamma = .005, alpha = 1, kappa = .4, tau = 280)

# plot original + degraded, leave room for MAP estimate
par(mfrow = c(1, 3), mar = c(2.6,  1, 2.6, 1))
display(ystar, "Original image")
display(y, "Noisy data")

x <- runGibbs(1000) # can be called multiple times to proceed further into the chain

display(x, "MAP estimate")
mse(ystar, x) / mse(ystar, y)

dyn.unload("speedy.dll")

##############################################################################
