#' The col_summarize function
#'
#' Function to quickly return credible intervals
#' @param MAT A matrix
#' @param dig Number of digits to round estimates and CrIs to
#' @param level Confidence level
#' @keywords Credible intervals
#' @export
#' @importFrom stats median quantile
#' @examples
#' M <- matrix(rnorm(1000),ncol = 4)
#' col_summarize(M)

col_summarize <- function(MAT,dig = 2,level = 0.95)
{
  k <- ncol(MAT)
  ret_vec <- rep(0,k)
  for(i in 1:k)
  {
    dat <- MAT[,i]
    est <- round(median(dat),dig)
    cri <- round(unname(quantile(dat,probs = c((1-level)/2,level + (1-level)/2))),dig)
    res <- paste(est," (",cri[1],", ",cri[2],")",sep = "")
    ret_vec[i] <- res
  }
  return(ret_vec)
}
