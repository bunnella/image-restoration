library(png)
library(jpeg)

# SET YOUR WORKING DIRECTORY TO THE REPO!!!
# setwd("K:/math400-10-w16/common/image-restoration")
# setwd("~/Carleton/MATH-COMPS")

##############################################################################
## Function defs

q <- .25

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
mse <- function(originalImg, RestoredImg) {
  n <- length(originalImg)
  returnVal <- 0
  for (i in 1:n) {
    returnVal <- returnVal + (originalImg[i] - RestoredImg[i])^2
  }
  returnVal/n
}

setupMCMC <- function(
  y, x,
  seed    = 0,
  theta   = 4,
  gamma   = .1,
  alpha   = 1.33,
  kappa   = 0.5,
  tau     = 180,
  omicron = 0.1,
  ...) .Call("R_setupMCMC", y, x, seed, theta, gamma, alpha, kappa, tau, omicron)

runMCMC <- function(N) .Call("R_runMCMC", N)

##############################################################################
## Setup

# read in the test image
original <- readPNG("img/startest.png")
R <- nrow(original)
C <- ncol(original)

# degrade away!
writeJPEG(original, target = "tmp.jpg", quality = q)
y <- extractGray(readJPEG("tmp.jpg"))
ystar <- extractGray(original)

##############################################################################
## Gibbs sampler!

dyn.load("speedy.dll")

x <- y[] # copy

setupMCMC(y, x, seed = 100, theta = 0.2, gamma = 0.03333333, alpha = 1.5)

# plot original + degraded, leave room for MAP estimate
par(mfrow = c(1, 3), mar = c(2.6,  1, 2.6, 1))
display(ystar, "Original image")
display(y, "Noisy data")

x <- runMCMC(100) # can be called multiple times to proceed further into the chain

display(x, "MAP estimate")

dyn.unload("speedy.dll")

##############################################################################
