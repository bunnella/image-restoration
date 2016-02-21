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
original <- readPNG("img/final/star80.png")
R <- nrow(original)
C <- ncol(original)

# degrade away!
writeJPEG(original, target = "tmp.jpg", quality = q)
y <- readJPEG("tmp.jpg")

##############################################################################
## Gibbs sampler!

dyn.load("speedy.dll")

nsteps <- 2000
fixed <- y

for (ch in 1:3) {
  y.ch <- y[ , , ch]
  x.ch <- y.ch[]

  setupGibbs(y.ch, x.ch, theta = .14, gamma = .01, alpha = 1, kappa = .3, tau = 350)
  fixed[ , , ch] <- runGibbs(nsteps)
}

writePNG(fixed, "fixed.png")

dyn.unload("speedy.dll")

##############################################################################
