library(png)
library(jpeg)

# SET YOUR WORKING DIRECTORY TO THE REPO!!!
setwd("~/Carleton/MATH-COMPS")

##############################################################################
## Function defs

q <- .2

V <- c()

display <- function(img, caption = "") {
  image(t(img[R:1, 1:C]), col=gray(V), zlim=0:1, frame=F, asp=C/R, xaxt="n", yaxt="n", main=caption)
}

setupGibbs <- function(
  y, x = y[],
  seed    = 0,
  nlevels = 64,
  theta   = 4,
  gamma   = .1,
  alpha   = 1.2,
  tau     = 180,
  omicron = 0.1,
  ...) {
  V <<- seq(0, 1, length.out = nlevels)
  .Call("R_setupGibbs", y, x, seed, V, theta, gamma, alpha, tau, omicron)
}

runGibbs <- function(N) .Call("R_runGibbs", N)

##############################################################################
## Setup

# read in the test image
picture <- readPNG("img/lena_gray_512.png")
R <- ncol(picture)
C <- nrow(picture)

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

setupGibbs(y, omicron = .075, tau = 1000, theta = 1.5, gamma = 0.01)

# plot original + degraded, leave room for MAP estimate
par(mfrow = c(1, 3), mar = c(2.6, 1, 2.6, 1))
display(original, "Original image")
display(y, "Noisy data")

x <- runGibbs(10) # can be called multiple times to proceed further into the chain

display(x, "MAP estimate")

dyn.unload("speedy.dll")
