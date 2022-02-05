require(mvtnorm)

MVNsampler <- function(N, m, S) {
  #' Return N samples from the a multivariate normal (MVN) density function
  #' @param N the number of realizations (samples)
  #' @param m mu, the mean of the multivariate normal (MVN) density function
  #' @param S sigma, the covariance matrix of the MVN. must be positive definite
  #' @examples
  #' mu = runif(2)
  #' cv.pre = matrix(runif(4,min=0.01,max=1.01),nrow=2)
  #' cv = 1/2 * (cv.pre + cv.pre)
  #' samples = MVNsampler(mu,cv,)
  #' @export

L = t(chol(S))
d = length(m)
samples.pre = mvtnorm::rmvnorm(N,mean = rep(0,d), sigma = diag(d));

  samples = t(m + L %*% t(samples.pre))
  return(samples)
}
