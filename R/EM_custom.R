# Problem 2

EM_custom <- function(n, p_vec, T_end) {
  #' Expectation-Maximization function for the blood type example
  #'
  #' @param n is the phenotype vector (n_a, n_b, n_ab, n_o)
  #' @param p_vec is the allele vector (p_a, p_b, p_o)
  #' @param T_end indicates how often the for-loop iterates
  #' @export
  n_a = n[1]; n_ab = n[2]; n_b = n[3]; n_o = n[4]; N = sum(n);
  p_a = p_vec[1]; p_b = p_vec[2]; p_o = p_vec[3];
  for (t in 1:T_end) {
    # new_theta = argmax_theta E(ln(Pr(x|theta))|y,theta_t)
    #p_tt = argmax(E(ln(Pr(g|p))|f,p_t))
    #if (abs(new_theta - theta) < epsilon || abs(Pr(f|p_tt)-Pr(f|p_t)) < sigma){
    #  return(result)
    #}
    # Expectation Step
    m_aa = n_a * (p_a^2)/(p_a^2+2*p_a*p_o); m_ao = n_a * (2*p_a*p_o)/(p_a^2+2*p_a*p_o);
    m_bb = n_b * (p_b^2)/(p_b^2+2*p_b*p_o); m_bo = n_b * (2*p_b*p_o)/(p_b^2+2*p_b*p_o);
    m_ab = n_ab; m_o = n_o;

    # Maximization Step
    p_a = (2*m_aa+m_ao+m_ab)/(2*N);
    p_b = (2*m_bb+m_bo+m_ab)/(2*N);
    p_o = (2*m_o+m_ao+m_bo)/(2*N);
  }
  result = c(p_a,p_b,p_o)
  return(result)
}
