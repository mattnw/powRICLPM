library(lavaan)

dat <- read.table("C:\\Users\\5879167\\surfdrive\\R packages\\powRICLPM\\examples\\RICLPM.dat",
                  col.names = c("x1", "x2", "x3", "x4", "x5", "y1", "y2", "y3", "y4", "y5"))

RICLPM <- '
  # Create between components (random intercepts)
  RIx =~ 1*x1 + 1*x2 + 1*x3 + 1*x4
  RIy =~ 1*y1 + 1*y2 + 1*y3 + 1*y4

  # Create within-person centered variables
  wx1 =~ 1*x1
  wx2 =~ 1*x2
  wx3 =~ 1*x3
  wx4 =~ 1*x4

  wy1 =~ 1*y1
  wy2 =~ 1*y2
  wy3 =~ 1*y3
  wy4 =~ 1*y4

  # Estimate the lagged effects between the within-person centered variables.
  wx2 ~ 0.4*wx1 + 0.3*wy1
  wy2 ~ 0.2*wx1 + 0.3*wy1
  wx3 ~ 0.4*wx2 + 0.3*wy2
  wy3 ~ 0.2*wx2 + 0.3*wy2
  wx4 ~ 0.4*wx3 + 0.3*wy3
  wy4 ~ 0.2*wx3 + 0.3*wy3

  # Estimate the covariance between the within-person centered variables at the first wave.
  wx1 ~~ 0.35*wy1 # Covariance

  # Estimate the covariances between the residuals of the within-person centered variables (the innovations).
  wx2 ~~ 0.117*wy2
  wx3 ~~ 0.117*wy3
  wx4 ~~ 0.117*wy4

  # Estimate the variance and covariance of the random intercepts.
  RIx ~~ 1.5*RIx
  RIy ~~ 1.5*RIy
  RIx ~~ 0.7875*RIy

  # Estimate the (residual) variance of the within-person centered variables.
  wx1 ~~ 1*wx1 # Variances
  wy1 ~~ 1*wy1
  wx2 ~~ 0.666*wx2 # Residual variances
  wy2 ~~ 0.828*wy2
  wx3 ~~ 0.666*wx3
  wy3 ~~ 0.828*wy3
  wx4 ~~ 0.666*wx4
  wy4 ~~ 0.828*wy4

  wx1 + wx2 + wx3 + wx4 + wy1 + wy2 + wy3 + wy4 ~ 0*1
'
RICLPM.fit <- lavaan(RICLPM, data = dat, missing = 'ML', meanstructure = F, int.ov.free = T)
summary(RICLPM.fit, standardized = T)
sigma <- lavInspect(RICLPM.fit, what = "cov.ov")
eigen(sigma)

t(chol(sigma))%*%chol(sigma)

library(MASS)
x.random <- mvrnorm(n = 1000000, rep(0, 8), diag(8))
x<- x.random%*%chol(sigma)
cov(x)
