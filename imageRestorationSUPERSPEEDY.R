library(png)
library(jpeg)

# SET YOUR WORKING DIRECTORY TO THE REPO!!!
setwd("K:/math400-10-w16/common/image-restoration")

##############################################################################
## Function defs

q <- .8

V <- c()

display <- function(img, caption = "") {
  image(t(img[R:1, 1:C]), col=gray(V), zlim=0:1, frame=F, asp=R/C, xaxt="n", yaxt="n", main=caption)
}

setupGibbs <- function(
  y, x,
  seed    = 0,
  nlevels = 64,
  theta   = 4,
  gamma   = .1,
  alpha   = 1.2,
  kappa   = 0.5,
  tau     = 180,
  omicron = 0.1,
  ...) {
  V <<- seq(0, 1, length.out = nlevels)
  .Call("R_setupGibbs", y, x, seed, V, theta, gamma, alpha, kappa, tau, omicron)
}

runGibbs <- function(N) .Call("R_runGibbs", N)

##############################################################################
## Setup

# read in the test image
picture <- readPNG("img/beiber.png")
R <- nrow(picture)
C <- ncol(picture)

# transform to image()-ready orientation
ch.R <- picture[ , , 1]
ch.G <- picture[ , , 2]
ch.B <- picture[ , , 3]
original <- 0.30*ch.R + 0.59*ch.G + 0.11*ch.B

# degrade away!
writeJPEG(original, target = "tmp.jpg", quality = q)
y <- readJPEG("tmp.jpg")

##############################################################################
## Gibbs sampler!

dyn.load("speedy.dll")

x <- y[] # copy

setupGibbs(y, x, omicron = .075, tau = 1000, theta = 2, gamma = 0.02)

# plot original + degraded, leave room for MAP estimate
par(mfrow = c(1, 3), mar = c(2.6, 1, 2.6, 1))
display(original, "Original image")
display(y, "Noisy data")

x <- runGibbs(10) # can be called multiple times to proceed further into the chain

display(x, "MAP estimate")

dyn.unload("speedy.dll")
