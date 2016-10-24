
## Linear Congruential Generator
# Generate N independent samples from the uniform distribution on [0, 1]
unif_gen <- function(N) {
  x0 <- 1
  a <- 39373
  c <- 0
  k <- 2^31 - 1
  u <- NULL;
  for(i in 1:N) {
    x0 = (a*x0 + c) %% k
    u = c(u, x0 / k)
  }
  return(u)
}

## Beasley–Springer–Moro algorithm to compute the inverse function of the std normal
# input u between 0 and 1
stdNormal_cdfInv <- function(u) {
  if(max(u) > 1 | min(u) < 0) {
    stop("Ileagal Input!")
  }
  a0 <- 2.50662823884
  a1 <- -18.61500062529
  a2 <- 41.39119773534
  a3 <- -25.44106049637
  b0 <- -8.47351093090
  b1 <- 23.08336743743
  b2 <- -21.06224101826
  b3 <- 3.13082909833
  c0 <- 0.3374754822726147
  c1 <- 0.9761690190917186
  c2 <- 0.1607979714918209
  c3 <- 0.0276438810333863
  c4 <- 0.0038405729373609
  c5 <- 0.0003951896511919 
  c6 <- 0.0000321767881768 
  c7 <- 0.0000002888167364 
  c8 <- 0.0000003960315187
  
  y <- u - 0.5
  x <- y
  r <- y
  idx <- (abs(y) < 0.42)  # idx - |y| < 0.42
  idx2 <- !idx            # idx2 - |y| > 0.42
  idx3 <- idx2 & (y>0)    # idx3 - |y| > 0.42 & y>0
  idx4 <- idx2 & (y<0)    # idx4 - |y| > 0.42 & y<0
  
  r[idx] <- y[idx]^2
  x[idx] <- y[idx] * (((a3*r[idx]+a2)*r[idx]+a1)*r[idx]+a0)/ 
    ((((b3*r[idx]+b2)*r[idx]+b1)*r[idx]+b0)*r[idx]+1)
  r[idx2] <- u[idx2]
  r[idx3] <- 1 - u[idx3]
  r[idx2] <- log(-log(r[idx2]))
  x[idx2] <- c0 +r[idx2]*(c1 +r[idx2]*(c2 +r[idx2]*(c3 +r[idx2]*(c4+
                                  r[idx2]*(c5 +r[idx2]*(c6 +r[idx2]*(c7 +r[idx2]*c8)))))))
  x[idx4] <- -x[idx4]
  return(x)
}

## Inverse Transform Method
# Generate N independent samples from the inverse function of the cdf F(x)
# F_inv(u) is the inverse function of F(x), must support vector calculation
ITM_gen <- function(N, F_inv = stdNormal_cdfInv) {
  F_u <- match.fun(F_inv)
  u <- unif_gen(N)
  return(F_u(u))
}

## Acceptance–Rejection Method to generate normal distribution, use the method from hw6
# the sample size generated will be less than N
ARM_gen_normal <- function(N) {
  
}
