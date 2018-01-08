
#' @export
qq_inflation <- function(pvals, df = 1, thr_pval,  ...)
{
  #if(!missing(thr_pval)) {
  #  pvals <- pvals[pvals > thr_pval]
  #}
  
  median(qchisq(pvals, df, lower.tail = FALSE)) / qchisq(0.5, df, lower.tail = FALSE)
}

#' @export
qq_theme <- function() theme_minimal()

#' @export
qq_plot <- function(pvals, df = 1,
  title, 
  size = NULL, color = NULL, linetype = 3, lims = NULL, ...)
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
    if(is.null(color)) {
      p <- p + geom_point() 
    } else {
      p <- p + geom_point(color = color) 
    }
  } else {
    if(is.null(color)) {
      p <- p + geom_point(size = size) 
    } else {
      p <- p + geom_point(color = color, size = size) 
    }
  }
  
  # line
  p <- p + geom_abline(linetype = linetype)
  
  # x/y limits
  if(is.null(lims)) {
    lims <- with(dat,
      c(min(val0[num_pvals], val[num_pvals]), max(val0[1], val[1])))
  }
  
  p <- p + xlim(lims) + ylim(lims)
  
  # labs
  if(!missing(title)) {
    p <- p + labs(title = title)
  }
  
  p <- p + labs(x = "Expected -log(pval)", y = "Observed -log(pval)")
  
  # subtitle
  lambda <- qq_inflation(pvals, df = df, ...)
  
  p <- p + labs(subtitle = paste0("inflation: ", round(lambda, 2)))
  
  
  # theme
  #p <- p + qq_theme()
              
  return(p)
}
