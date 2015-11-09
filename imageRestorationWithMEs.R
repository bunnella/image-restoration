# image
# Gibbs Sample/ Simulated Annealing
# RAP
# Output

display <- function(img, caption = "", ColVals = (V/2+.5)) {
  img <- img/2+.5
  image(img, col = gray(ColVals), frame = F, asp = 1, xaxt = "n", yaxt = "n", main = caption)
}

degrade <- function(original, perturb_percentage=.2) {
  static_vector <- rnorm(R*C,0,perturb_percentage/2)
  # 90% within two standard deviations
  # which is +/- perturb_percentage
  # which is perturb_percentage as percentage of dist between -1 and 1
  random_noise <- matrix(static_vector, nrow=R, ncol=C)
  # print(head(random_noise))
  added_noise <- original+random_noise
  #added_noise
  m <- matrix(mapply(function(x) min(max(x,-1), 1),added_noise),nrow=R,ncol=C)
  m
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

f <- function(xs, xt, e) {
	return((xs-xt)^2*(1-e) + gamma*e)
}

createV <- function(num) {
	seq(-1,1,length.out = num)
}

# data structures/functions for manipulating microedges
eNS <- matrix(0, R-1, C)
eEW <- matrix(0, R, C-1)
getME <- function(s1, s2) {
  r1 <- convertItoRC(s1)[1]; r2 <- convertItoRC(s2)[1];
  c1 <- convertItoRC(s1)[2]; c2 <- convertItoRC(s2)[2];
  if (abs(r1-r2) == 1) eNS[min(r1, r2), c1]
  else if (abs(c1-c2) == 1) eEW[r1, min(c1, c2)]
}
setME <- function(s1, s2, v = 1) {
  r1 <- convertItoRC(s1)[1]; r2 <- convertItoRC(s2)[1];
  c1 <- convertItoRC(s1)[2]; c2 <- convertItoRC(s2)[2];
  if (abs(r1-r2) == 1) eNS[min(r1, r2), c1] <- v
  else if (abs(c1-c2) == 1) eEW[r1, min(c1, c2)] <- v
}

# energy function
H <- function(s, v) {
  theta*d(v, y[s]) + sum(sapply(neighbors(s), function(t) f(v, x[t], getME(s, t))))
}

# returns a sample from the local characteristic distribution
sampleXs <- function(s, beta = 1) {
	probs <- sapply(V, function(v) exp(-beta*H(s, v)))
	sample(V, 1, prob = probs)
}

# returns a sampled microedge parity between pixels s & t
sampleME <- function(s, t, beta = 1) {
  probs <- sapply(0:1, function(e) exp(-beta*f(s, t, e)))
  sample(0:1, 1, prob = probs)
}

R <- C <- 20
noiseProb <- .2
theta <- 8
gamma <- .1
N <- 20^1
V <- createV(40)

#blah <- c(rep(-1, times=20), rep(c(-1,1,1,1,-1),times=12), rep(c(-1, rep(1,times=8), -1),times=2))
#original <- matrix(rep(c(blah, rev(blah)), each=2), C, R)

library(bmp)
picture <- read.bmp("small_cat.bmp")
values <- picture[ , , 1]
values <- rbind(apply(cbind(apply(values, 2, rev)), 1, I))
original <- values / 127.5 - 1
R <- nrow(original)
C <- ncol(original)

y <- degrade(original, perturb_percentage=.25)
display(y)

x <- y

par(mfrow = c(1, 3), mar = c(2.6, 1, 2.6, 1))
display(original, "Original")
display(y, "Noisy data")

# Gibbs time ;)
for (n in 1:N) {
  cat("0%")
  beta = min(exp((n - N/2)/20), 50)
  for (s in 1:(R*C)) {
    x[s] = sampleXs(s, beta)
  }
  for (s in 1:(R*C)) {
    for (t in neighbors(R*C)) {
      if (t > s) {
        setME(s, t, sampleME(s, t, beta))
      }
    }
  }
  cat(paste0("\r", round(100*n/N), "%\r"))
}

display(x, "MAP estimate")



