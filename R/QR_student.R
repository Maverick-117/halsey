QR_student <- function(x,y){
  #' QR decomposition-based OLS
  #'
  #' @param X covariates
  #' @param y response
  #' @return regression coefficients c
  #' @export
  q_r <- qr(x)
  Q <- qr.Q(q_r)
  R <- qr.R(q_r)
  c<-solve(R,t(Q)%*%y)
  return(c)
}
