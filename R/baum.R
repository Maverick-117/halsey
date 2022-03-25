baum <- function(em_row, y, T_end, K_end) {
  #' Implements the Baum-Welch algorithm to guess initial distribution, transition probabilities, and
  #' second row of the emission matrix
  #'
  #' @param em_row the row for the fair state
  #' @param y the history of observable actions
  #' @param T_end number of timesteps in y
  #' @param K_end number of iterations of times the Baum-Welch algorithm is looped
  #' @export
  em_row_rand = rep(1/6,6);
  tp = matrix(rep(0.5,4),nrow=2);
  em = matrix(c(em_row, em_row_rand));
  v0 = rep(1/2,2);
  v = matrix(rep(0,K_end*2),nrow=2);
  P = array(rep(0,K_end*4),c(2,2,K_end));
  E = array(rep(0,K_end*12),c(2,6,K_end));

  for (k in 1:K_end) {
    E[1,,k] = em_row;
  }

  v[,1] = v0; P[,,1] = tp; E[2,,1] = em_row_rand;

  g = array(rep(0,T_end*4),c(2,2,T_end));

  for (k in 1:(K_end-1)){
    a = HMMforward(y,P[,,k], E[,,k], v[,k]);
    b = HMMbackward(y,P[,,k], E[,,k], v[,k]);
    gamma = hiddenStateInference(a,b,T_end);
    # calculating g
    for (t in 2:T_end) {
      for (i in 1:2){
        for (j in 1:2) {
          g[i,j,t] = a[i,t-1] * P[i,j,k] * E[j,y[t],k] * b[j,t];
        }
      }
      g[,,t] = g[,,t]/sum(g[,,t])
    }
    # calculating v_k and P_k
    for (i in 1:2) {
      v[i,k+1] = gamma[i,1];
      for (j in 1:2) {
        P[i,j,k+1] = sum(g[i,j,]);
      }
      P[i,,k+1] = P[i,,k+1]/sum(P[i,,k+1])
    }
    # calculating E_k for the unfair state
    for (l in 1:6) {
      E[2,l,k+1] = sum(gamma[2,]*(y == l));
    }
    E[2,,k+1] = E[2,,k+1]/sum(gamma[2,]);
  }
  return(list(v,P,E))
}
