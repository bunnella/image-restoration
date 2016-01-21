library(bmp)
library(compiler)
invisible(enableJIT(3))

# SET YOUR WORKING DIRECTORY TO THE REPO!!!
setwd("~/Carleton/MATH-COMPS/")

##############################################################################
## Function defs

V <- c()

display <- function(img, caption = "") {
  image(img/2+.5, col=gray(V/2+.5), zlim=0:1, frame=F, asp=C/R, xaxt="n", yaxt="n", main=caption)
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
  V <<- seq(-1, 1, length.out = nlevels)
  .Call("R_setupGibbs", y, x, seed, V, theta, gamma, alpha, tau, omicron)
}

runGibbs <- function(N) .Call("R_runGibbs", N)

##############################################################################
## Setup

# read in the test image
picture <- read.bmp("img/lena_gray_512b.bmp")
R <- ncol(picture)
C <- nrow(picture)

# transform to image()-ready orientation and degrade
ch.R <- t(picture[C:1, 1:R, 1])
ch.G <- t(picture[C:1, 1:R, 2])
ch.B <- t(picture[C:1, 1:R, 3])
values <- 0.30*ch.R + 0.59*ch.G + 0.11*ch.B
original <- values / 127.5 - 1 # scale from [0..255] -> [-1, 1]
y <- original

##############################################################################
## Gibbs sampler!

dyn.load("speedy.dll")

setupGibbs(y, omicron = .09, tau = 200, theta = 3, gamma = .05, alpha = 1)

# plot original + degraded, leave room for MAP estimate
par(mfrow = c(1, 2), mar = c(2.6, 1, 2.6, 1))
display(y, "Noisy data")

x <- runGibbs(50)

display(x, "MAP estimate")

dyn.unload("speedy.dll")
