# NOTE: I'm an extension of the normal imageRestoration script; run most of it before me!

error <- sqrt(abs(y - original))

minimax <- function(d) {
  T <- .5
  min(d^2, T^2+2*T^2*(abs(d)-T))
}

trust <- matrix(0, nrow = R, ncol = C)
for (r in 2:(R-1))
  for (c in 2:(C-1))
    for (i in c(-1, 1))
      for (j in c(-1, 1))
        trust[r, c] <- trust[r, c] + minimax(y[r, c] - y[r+i, c+j])


par(mfrow = c(1, 2))
display(error, "JPEG error")
display(sqrt(trust), "Contrast measure")
