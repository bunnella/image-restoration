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

# Mean Square Error
mse <- function(originalImg, RestoredImg) {
  n <- length(originalImg)
  returnVal <- 0
  for (i in 1:n) {
    returnVal <- returnVal + (originalImg[i] - RestoredImg[i])^2
  }
  returnVal/n
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
picture <- readPNG("img/startest.png")
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
## Search for Gibbs sampler parameters!

dyn.load("speedy.dll")

x <- y[] # copy
#theta, gamma, alpha , kappa, tau, omicron 
numTrials <- 100
numGuesses <- 10
parameters <-  c(2, .02, 1.2, 0.5, 1000, 0.075)
prevBestMSE <- 999999
paramSDs <- parameters/4
lp <- length(parameters)
for (i in 1:(lp*numTrials)) {
  index <- i%%lp
  guesses <- rnorm(numGuesses, parameters[index], paramSDs[index])
  origParam <- parameters[index]
  bestGuess <- 0
  bestMSE <- 999999
  for (j in 1:numGuesses) {
    parameters[index] <- guesses[j]
    setupGibbs(y, x, theta = parameters[1], gamma = parameters[2],
      alpha = parameters[3], kappa = parameters[4]
      tau = parameters[5], omicron = parameters[6])
    x <- runGibbs(100)
    tempMSE <- mse(original, x)
    if (tempMSE < bestMSE) {
      bestGuess <- j
    }
  }
  if (bestMSE < prevBestMSE) {
    parameters[index] <- guesses[bestGuess]
  } else {
    paramSDs[index] <- paramSDs[index]/2
  }
}









# plot original + degraded, leave room for MAP estimate
par(mfrow = c(1, 3), mar = c(2.6,  1, 2.6, 1))
display(original, "Origina1 image")
display(y, "Noisy data")

x <- runGibbs(10) # can be called multiple times to proceed further into the chain

display(x, "MAP estimate")

dyn.unload("speedy.dll") #speedy dot dill


##############################################################################
