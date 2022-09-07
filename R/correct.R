
#' @export
qq_correct <- function(pvals, ...)
{
  stats <- qchisq(pvals, 1, lower.tail = FALSE)
  
  inflation <- median(stats, na.rm = TRUE) / qchisq(0.5, 1)
  
  ### correction
  pvals_corrected <- pvals
  if(inflation > 1) {
    stats_corrected <- stats / inflation
    pvals_corrected <- pchisq(stats_corrected, 1, lower.tail = FALSE)
  }
  
  return(pvals_corrected)
}

