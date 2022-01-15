mult2mv <- function(A, B, x, swtch) {
  #' Compare two different linear mappings.
  #'
  #' Multiplying matrices and vectors is associative. This function performs multiplies matrix A (mxp), matrix B(pxn), and vector x (nx1) in two different ways.
  #' @param A The left matrix. (m x p)
  #' @param B The right matrix. (p x n)
  #' @param x The vector. (n x 1)
  #' @param swtch Boolean Flag to use one or the other associative ordering of multiplication.
  #' @examples
  #' mult2mv(matrix(1:4,2,2),matrix(5:8,2,2),c(2,1),TRUE)
  #' mult2mv(matrix(1:8,2,4),matrix(8:19,4,3),c(3,2,1),FALSE)
  #' @export
  if(swtch == TRUE) {
    # A * (B * x)
    result <- A %*% (B %*% x)
  }
  else {
    # (A * B) * x
    result <- (A %*% B)%*% x
  }
  return(result)
}
