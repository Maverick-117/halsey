
MCMC_beta <- function(x0,N,alph=4,bet=2,sigma=0.125) {
  #' Apply the MCMC to the Beta distribution with a multiplicative log-normal proposal density.
  #'
  #' @param x0 Initial guess
  #' @param N how many iterations for which the MCMC should be ran
  #' @export
  mean_trace = rep(0,N);
  MCMC_trace = rep(0,N);
  stdmean_trace = rep(0,N);
  MCMC_trace[1] = x0;
  mean_trace[1] = x0;
  stdmean_trace[1] = 0;
  for (n in 2:N) {

    x1_temp = rlnorm(1,meanlog=log(MCMC_trace[n-1]),sdlog=sigma);
    hastings = dlnorm(MCMC_trace[n-1],meanlog=log(x1_temp),sdlog=sga)*dbeta(x1_temp,alph,bet)/(dlnorm(x1_temp,meanlog=log(MCMC_trace[n-1]),sdlog=sga)*dbeta(MCMC_trace[n-1],alph,bet));
    #hastings = ((x1_temp^4) * (1-x1_temp))/((MCMC_trace[n-1]^4)*(1-MCMC_trace[n-1]));
    acc = min(hastings,1);
    u = runif(1);
    if (u <= acc) {
      x1 = x1_temp;
    }
    else {
      x1 = MCMC_trace[n-1];
    }
    MCMC_trace[n] = x1;
    mean_trace[n] = mean(MCMC_trace[1:n])
    stdmean_trace[n] = sd(mean_trace[1:n]);
  }
  return(list(MCMC_trace,mean_trace,stdmean_trace))
}
