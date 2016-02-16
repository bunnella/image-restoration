library(png)
library(jpeg)

# SET YOUR WORKING DIRECTORY TO THE REPO!!!
# setwd("K:/math400-10-w16/common/image-restoration")
# setwd("~/Carleton/MATH-COMPS")

##############################################################################
## Function defs

q <- .5

extractGray <- function(pic) {
  ch.R <- pic[ , , 1]
  ch.G <- pic[ , , 2]
  ch.B <- pic[ , , 3]
  0.30*ch.R + 0.59*ch.G + 0.11*ch.B
}

display <- function(img, caption = "") {
  image(t(img[R:1, 1:C]), col=gray(0:255/255), zlim=0:1, frame=F, asp=R/C, xaxt="n", yaxt="n", main=caption)
}

# Mean Square Error
mse <- function(ystar, xstar) {
  mean((ystar-xstar)^2)
}

displayError <- function(ystar, xstar, mse = FALSE) {
  diff <- ystar-xstar
  if (mse) {
    diff <- diff^2
  } else {
    diff <- (diff^2)^.5
  }
  diff <- 1-diff
  display(diff)
}

setupGibbs <- function(
  y, x,
  seed  = 0,
  theta = 4,
  gamma = .1,
  alpha = 1.33,
  kappa = 0.5,
  step  = 1/256,
  tau   = 100,
  ...) .Call("R_setupGibbs", y, x, seed, theta, gamma, alpha, kappa, step, tau)

runGibbs <- function(N) .Call("R_runGibbs", N)

numSteps <- function() {
  cat(sprintf("Completed %d chain steps.\n", .Call("R_numSteps")))
}

##############################################################################
## Setup

# read in the test image
original <- readPNG("img/final/star80.png")
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

setupGibbs(y, x, theta = 0.17, gamma = 0.011, alpha = 0.83, kappa = .41, tau = 282)

# plot original + degraded, leave room for MAP estimate
par(mfrow = c(1, 3), mar = c(2.6,  1, 2.6, 1))
display(ystar, "Original image")
display(y, "Noisy data")

x <- runGibbs(2000) # can be called multiple times to proceed further into the chain

display(x, "MAP estimate")
mse(ystar, x) / mse(ystar, y)

dyn.unload("speedy.dll")

##############################################################################
