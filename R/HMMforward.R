HMMforward <- function(y,P,E,v0) {
  #' Forward algorithm for HMM solving
  #' outputs an N_s x N_T matrix, where N_s is the number of hidden states and
  #' N_T is the number of timesteps
  #'
  #' @param y observable actions taken
  #' @param P transition probabilities for hidden states
  #' @param E emission matrix for observable actions, given hidden actions
  #' @param v0 initial guess
  #'
  #' @export
  n = length(y);
  s = length(v0);
  result = rep(1,s*n,nrow=s);
  a1 = c(v0[1]*E[1,y[1]],v0[2]*E[2,y[1]]);
  a = matrix(rep(0,n*2),nrow=2); a[,1] = a1;
  for (t in 2:n) {
    for (i in 1:s) {
      a[i,t] = E[i,y[t]] * sum(a[,t-1]*P[i,])
    }
  }
  return(a)
}
