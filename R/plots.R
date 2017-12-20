
#' @export
qq_inflation <- function(pvals, thr_pval, ...)
{
  #if(!missing(thr_pval)) {
  #  pvals <- pvals[pvals > thr_pval]
  #}
  
  median(qchisq(pvals, 1, lower.tail = FALSE)) / qchisq(0.5, 1, lower.tail = FALSE)
}

#' @export
qq_theme <- function() theme_minimal()

#' @export
qq_plot <- function(pvals, title, size = NULL, linetype = 3, ...)
{
  ### filter `pvals`
  pvals <- pvals[!is.na(pvals)]
  
  num_pvals <- length(pvals)
  
  ### prepare `pvals`
  pvals <- sort(pvals)
  
  pvals0 <- ppoints(length(pvals))
  
  ### plot
  dat <- data_frame(val0 = -log10(pvals0), val = -log10(pvals))
  p <- ggplot(dat, aes(val0, val)) 
  
  if(is.null(size)) {
    p <- p + geom_point() 
  } else {
    p <- p + geom_point(size = size) 
  }
  
  # line
  p <- p + geom_abline(linetype = linetype)
  
  # x/y limits
  lims <- with(dat,
    c(min(val0[num_pvals], val[num_pvals]), max(val0[1], val[1])))
  
  p <- p + xlim(lims) + ylim(lims)
  
  # labs
  if(!missing(title)) {
    p <- p + labs(title = title)
  }
  
  p <- p + labs(x = "Expected -log(pval)", y = "Observed -log(pval)")
  
  # subtitle
  lambda <- qq_inflation(pvals, ...)
  
  p <- p + labs(subtitle = paste0("inflation: ", round(lambda, 2)))
  
  
  # theme
  #p <- p + qq_theme()
              
  return(p)
}
