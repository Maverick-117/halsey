logres_FS <- function(x,y,theta,alpha,T_end) {
  #' Logistic regression optimization with Fisher-scoring
  #' 
  #' @param x the covariate variables
  #' @param y the response variable
  #' @param theta the initial parameter estimate
  #' @param alpha the learning rate
  #' @param T_end the number of iterations for which the optimization is allowed to run
  #' @export
  LL_history = c();
  conf_int = matrix(rep(0,T_end*length(theta)),nrow=length(theta));
  for (t in 1:T_end) {
    p = 1/(1+exp(-(x %*% theta)))
    
    loss = -(((t(y) %*% log(p)) + (t(1-y) %*% log(1-p)) ));
    LL_history = c(LL_history, loss);
    # the gradient descent part
    # new_theta = theta - ( alpha/N * ( t(x) %*% (p - y) ))
    
    # fisher scoring part
    var_mat = diag(as.vector(p)*as.vector((1-p)));
    f_i = t(x) %*% var_mat %*% x; f_i_inv = chol2inv(chol(f_i));
    new_theta = theta - alpha * f_i_inv %*% ( t(x) %*% (p - y) );
    theta = new_theta
    conf_int[,t] = 1.96 * diag(sqrt(f_i_inv))
  }
  
  return(list(theta, conf_int, LL_history))
}