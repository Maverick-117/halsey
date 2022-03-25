logres_GD <- function(x,y,theta,alpha,T_end) {
  #' Logistic regression optimization with gradient-descent
  #' 
  #' @param x the covariate variables
  #' @param y the response variable
  #' @param theta the initial parameter estimate
  #' @param alpha the learning rate
  #' @param T_end the number of iterations for which the optimization is allowed to run
  #' @export
  LL_history = c();
  N = nrow(x);
  visited_thetas = matrix(rep(0,T_end*length(theta)),nrow=length(theta))
  for (t in 1:T_end) {
    p = 1/(1+exp(-(x %*% theta)))
    
    loss = -(((t(y) %*% log(p)) + (t(1-y) %*% log(1-p)) ));
    LL_history = c(LL_history, loss);
    new_theta = theta - ( alpha/N * ( t(x) %*% (p - y) ))
    theta = new_theta
    visited_thetas[,t] = theta;
  }
  conf_int = rep(0,length(theta));
  for (i in 1:length(theta)){
    conf_int[i] = 1.96 * sd(visited_thetas[i,(T_end-20):T_end])/sqrt(21);
  }
  return(list(theta, conf_int, LL_history))
}