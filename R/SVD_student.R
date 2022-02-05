SVD_student <- function(x,y){
  #' SVD-decomposition based OLS
  #'
  #' @param X covariates
  #' @param y response
  #' @return regression coefficients c
  #' @export
  UDV<-svd(x)
  U<-UDV$u
  D<-UDV$d
  D<-1/D
  Ddiag<-diag(D)
  V<-UDV$v
  c <- V %*% Ddiag %*% t(U) %*% y
  return(c)
}
