# image
# Gibbs Sample/ Simulated Annealing
# RAP
# Output

display <- function(img, caption = "", ColVals = (V/2+.5)) {
  img <- -img/2+.5
  image(img, col = gray(ColVals), frame = F, asp = 1, xaxt = "n", yaxt = "n", main = caption)
}

degrade <- function(original, noiseProb=.1, V = createV(1)) {
  degraded <- original
  for (i in 1:(R*C)) {
    flip <- sample(c(T,F), 1, prob=c(noiseProb, 1-noiseProb))
    degraded[i] <- ifelse(flip, sample(V,1), degraded[i])

  }
  degraded
}

neighbors <- function(i) {
	m <- ((i-1) %% R) + 1
	n <- ceiling(i / R)
	return(neighborsHelp(m, n))
}

neighborsHelp <- function(m, n) {
	neighborList <- c()
	if (m > 1) {
		neighborList <- c(neighborList, convertRCtoI(m-1, n))
	}
	if (n > 1) {
		neighborList <- c(neighborList, convertRCtoI(m, n-1))
	}
	if (m < R) {
		neighborList <- c(neighborList, convertRCtoI(m+1, n))
	}
	if (n < C) {
		neighborList <- c(neighborList, convertRCtoI(m, n+1))
	}
	return(neighborList)
}

convertItoRC <- function(i) {
	m <- ((i-1) %% R) + 1
	n <- ceiling(i / R)
	return(c(m,n))
}

convertRCtoI <- function(m, n) {
	return((m + (n-1)*R))
}

# euclidian distance
d <- function(xs, ys) {
	return((xs-ys)^2)
}

f <- function(xs,xt) {
	return((xs-xt)^2)
}

createV <- function(num) {
	seq(-1,1,length.out = num)
}

# data structures/functions for manipulating microedges
eNS <- matrix(0 R-1, C)
eEW <- matrix(0, R, C-1)
getME <- function(s1, s2) {
  r1 <- convertItoRC(s1)[0]; r2 <- convertItoRC(s2)[0];
  c1 <- convertItoRC(s1)[1]; c2 <- convertItoRC(s2)[1];
  if (abs(r1-r2) == 1) eEW(min(r1, r2))
  else if (abs(c1-c2) == 1) eNS(min(c1, c2))
  else stop("invalid neighbors, stupid")
}
setME <- function(s1, s2, v = 1) {
  r1 <- convertItoRC(s1)[0]; r2 <- convertItoRC(s2)[0];
  c1 <- convertItoRC(s1)[1]; c2 <- convertItoRC(s2)[1];
  if (abs(r1-r2) == 1) eEW(min(r1, r2)) <- v
  else if (abs(c1-c2) == 1) eNS(min(c1, c2)) <- v
  else stop("invalid neighbors, stupid")
}


# returns a sample from the local characteristic distribution
sampleXs <- function(s, beta = 1) {
	probs <- sapply(V, function(v)
		exp(-beta*(theta*d(v,y[s]) + sum(sapply(neighbors(s), function(t)
			f(v,x[t])
		))))
	)
	sample(V, 1, prob = probs)
}


#R <- C <- 20
noiseProb <- .2
theta <- 8
N <- 10^2
V <- createV(256)

#blah <- c(rep(-1, times=20), rep(c(-1,1,1,1,-1),times=12), rep(c(-1, rep(1,times=8), -1),times=2))
#original <- matrix(rep(c(blah, rev(blah)), each=2), C, R)

library(bmp)
picture <- read.bmp("hi.bmp")
values <- picture[,,1]
values <- rbind(apply(cbind(apply(values, 2, rev)), 1, I))
original <- values
R <- nrow(original)
C <- ncol(original)

y <- degrade(original, noiseProb, V)

x <- y

par(mfrow = c(1, 3), mar = c(2.6, 1, 2.6, 1))
display(original, "Original")
display(y, "Noisy data")

# Gibbs time ;)
for (n in 1:N) {
  for (s in 1:(R*C)) {
    x[s] = sampleXs(s, min(exp((n - N/2)/20), 50))
  }
  for (s in 1:(R*C)) {
    for (t in neighbors(R*C)) {
      if (t > s) {
        # compute energy with microedge on
        # compute energy with microedge off
        # return lower with exp prob
      }
    }
  }
}

display(x, "MAP")



