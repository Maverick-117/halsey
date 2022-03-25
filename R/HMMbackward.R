HMMbackward <- function(y,P,E,v0) {
  #' Backward algorithm for HMM solving
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
  b = matrix(rep(0,n*2),nrow=2); b[,n] = rep(1,s);
  for (t in seq(from=n-1,to=1,by=-1)) {
    for (i in 1:s) {
      b[i,t] = sum(P[,i]*E[,y[t+1]]*b[,t+1])
    }
  }
  return(b)
}
