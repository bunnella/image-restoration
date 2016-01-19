# image
# Gibbs Sample/ Simulated Annealing
# RAP
# Output

library(bmp)
library(compiler)
invisible(enableJIT(3))

set.seed(57)

# SET YOUR WORKING DIRECTORY TO THE REPO!!!
setwd("~/Carleton/MATH-COMPS/")
dyn.load("speedy.dll")

display <- function(img, caption = "") {
  image(img/2+.5, col=gray(V/2+.5), zlim=0:1, frame=F, asp=C/R, xaxt="n", yaxt="n", main=caption)
}

setupGibbs <- function(
  y, x = y[],
  theta   = 4,
  gamma   = .1,
  alpha   = 1.2,
  tau     = 150,
  omicron = 5,
...) { .Call("R_setupGibbs", y, x, V, theta, gamma, alpha, tau, omicron) }

runGibbs <- function(N) .Call("R_runGibbs", N)

##############################################################################
## Gibbs Sampler!

# read in the test image
picture <- read.bmp("Anne's Stuff!/img/lena_gray_512b.bmp")
R <- ncol(picture)
C <- nrow(picture)

# transform to image()-ready orientation and degrade
ch.R <- t(picture[C:1, 1:R, 1])
ch.G <- t(picture[C:1, 1:R, 2])
ch.B <- t(picture[C:1, 1:R, 3])
values <- 0.30*ch.R + 0.59*ch.G + 0.11*ch.B
original <- values / 127.5 - 1 # scale from [0..255] -> [-1, 1]
y <- original

# setup gray levels
nlevels = 32
V <- seq(-1, 1, length.out = nlevels)

# plot original + degraded, leave room for MAP estimate
par(mfrow = c(1, 3), mar = c(2.6, 1, 2.6, 1))
display(original, "Original")
display(y, "Noisy data")

setupGibbs(y)

x <- runGibbs(1)

display(x, "MAP estimate")

dyn.unload("speedy.dll")
