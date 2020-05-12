#' The mean_CRI function
#'
#' Simple function to return the mean (95% CrI) for a vector
#' @param y A numeric vector
#' @param dig The number of digits to round to
#' @keywords MCMC
#' @export
#' @examples
#' mean_CRI(rnorm(1000))

mean_CRI <- function(y, dig = 2)
{
    m = round(mean(y, na.rm = TRUE),dig)
    cri = round(quantile(y, probs = c(0.025,0.95)), dig)
    ret = paste0(m," (",cri[1],", ",cri[2],")")
    return(ret)
}
