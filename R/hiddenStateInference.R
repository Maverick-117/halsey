hiddenStateInference <- function(a,b,T_end){
  #' Marginal probability calculator for hidden states, given forward and backward HMM solver outputs
  #' outputs an N_s x N_T matrix, where N_s is the number of hidden states and
  #' N_T is the number of timesteps
  #'
  #' @param a forward HMM solver output
  #' @param b backward HMM solver output
  #' @param T_end total number of timesteps
  #' @export
  N_s = nrow(a);
  mP = rep(0,T_end);
  Pb = rep(0,T_end);
  marg = matrix(rep(0,2*T_end),nrow=2);
  for (t in 1:T_end) {
    denom = sum(a[,t]*b[,t]);
    for (s in 1:N_s) {
      marg[s,t] = (a[s,t]*b[s,t])/denom
    }
    # mP[t] = which.max(pre_mP)
    # Pb[t] = max(pre_mP)
  }
  return(marg)
}
