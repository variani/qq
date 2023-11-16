
#' @export
qq_inflation <- function(pvals, df = 1, thr_pval,  ...)
{
  #if(!missing(thr_pval)) {
  #  pvals <- pvals[pvals > thr_pval]
  #}
  
  pvals <- pvals[!(is.infinite(pvals) | is.na(pvals))]
  #pvals[pvals == 0] <- min(pvals[pvals != 0])
  
  median(qchisq(pvals, df, lower.tail = FALSE)) / qchisq(0.5, df, lower.tail = FALSE)
}

#' @export
qq_theme <- function() theme_minimal()

fun_format_big = function(x) formatC(x, format = "d", big.mark = ",")

# https://stackoverflow.com/a/47932381
# example call with fontsize:
# > p + annotation_compass(lab, 'NW', gp = grid::gpar(fontsize = 8))
annotation_compass = function(label, position = c('N','NE','E','SE','S','SW','W','NW'), padding = grid::unit(c(0.5,0.5),"line"), ...)
{ 
  position = match.arg(position)
  x = switch(position, N = 0.5, NE = 1, E = 1, SE = 1, S = 0.5, SW = 0, W = 0, NW = 0)
  y = switch(position, N = 1, NE = 1, E = 0.5, SE = 0, S = 0, SW = 0, W = 0.5, NW = 1)
  hjust = switch(position, N = 0.5, NE = 1, E = 1, SE = 1, S = 0.5, SW = 0, W = 0, NW = 0)
  vjust = switch(position, N = 1, NE = 1, E = 0.5, SE = 0, S = 0, SW = 0, W = 0.5, NW = 1)
  f1 = switch(position, N = 0, NE = -1, E = -1, SE = -1, S = 0, SW = 1, W = 1, NW = 1)
  f2 = switch(position, N = -1, NE = -1, E = 0, SE = 1, S = 1, SW = 1, W = 0, NW = -1)
  annotation_custom(grid::textGrob(label,
     x=grid::unit(x,"npc") + f1*padding[1] , y=grid::unit(y,"npc") + f2*padding[2],
     hjust=hjust,vjust=vjust, ...))
}

#' @export
qq_plot <- function(pvals, df = 1, pval_min = NULL,
  title, tests = NULL,
  size = NULL, linetype = 3, lims = NULL, 
  log = FALSE, square = FALSE, thin = TRUE,
  pos_legend = c(0.05, 0.95),
  full_stats = FALSE,
  ...)
{
  ### conver to table
  pvals_input = pvals
  if(NCOL(pvals) == 1) {
    dat = tibble(val = pvals)
  } else {
    dat = pvals
  }
  vals_test = colnames(dat)

  ### log scale?
  if(log) {
    dat = mutate_all(dat, function(x) 10^(-x))
  }
  ### cap?
  if(!is.null(pval_min)) {
    dat = mutate_all(dat, function(x) ifelse(x < pval_min, pval_min, x))
  }

  ## inflation values
  n_dat = nrow(dat)
  vals_lambda = sapply(dat, qq_inflation)

  ### filter `pvals`
  dat = apply(dat, 2, function(pvals) {
    pvals = ifelse(is.infinite(pvals) | is.na(pvals), NA, 
      ifelse(pvals < 0, NA, pvals))
    ifelse(pvals %in% 0, min(pvals, na.rm = TRUE), pvals)
  }) %>% bind_cols
  dat = drop_na(dat)
  
  num_pvals <- nrow(dat)
  if(num_pvals == 0) {
    return(NULL)
  }
  num_test = ncol(dat)
  
  ### prepare `pvals`
  dat = mutate_all(dat, sort)

  dat = pivot_longer(dat, colnames(dat))
  names(dat) = c('test', 'val')

  dat = mutate(dat, val0 = ppoints(n()))

  # # A tibble: 10 × 3
  #  test     val   val0
  #  <chr>  <dbl>  <dbl>
  # 1 val   0.0857 0.0610
  # 2 val   0.104  0.159

  ## test names
  if(!is.null(tests)) {
    dat = mutate(dat, test = factor(test, vals_test, tests))
  } else {
    tests = vals_test
  }
  ## plot CI
  n <- num_pvals
  alpha <- 0.05; np <- 1e3; np = min(np, n - 1)
  gcol <- 'grey80'; galpha = 0.75
  cmat <- matrix(nrow = 2*np, ncol = 2)
  for(i in seq(np)) {
    cmat[i, 1] <- cmat[np*2 + 1 - i, 1] <- -log10((i - 0.5)/n)
    cmat[i, 2] <- -log10(qbeta(1 - alpha/2, i, n - i))
    cmat[np*2 + 1 - i, 2] <- -log10(qbeta(alpha/2, i, n - i))
  }
  cmat <- as.data.frame(cmat)
  colnames(cmat) <- c('x', 'y')
  # cmat = cmat[chull(cmat$x, cmat$y), ]
  # return(cmat)
  p <- ggplot() + geom_polygon(data = cmat, aes(x, y), color = NA, fill = gcol, alpha = galpha)

  ## plot points
  dat = mutate(dat, val0 = -log10(val0), val = -log10(val))
  if(thin & nrow(dat) > 1e3) {
    dat = mutate(dat, val0 = round(val0, 3)) %>% filter(!duplicated(val0))
  }
  
  if(is.null(size)) {
    p <- p + geom_point(data = dat, aes(val0, val, group = test, color = test)) 
  } else {
    p <- p + geom_point(data = dat, aes(val0, val, group = test, color = test), size = size) 
  }
  
  # line
  p <- p + geom_abline(linetype = linetype)
  
  # x/y limits
  if(square) {
  # if(is.null(lims)) {
    lims <- with(dat,
      c(min(val0[num_pvals], val[num_pvals]), max(val0[1], val[1])))
    lims = with(cmat,
      c(min(lims[1], x, y), max(lims[2], x, y)))
    p <- p + xlim(lims) + ylim(lims)
  }
  
  # labs
  if(!missing(title)) {
    p <- p + labs(title = title)
  }
  
  p <- p + labs(x = TeX("Expected $-log_{10}(P)$"), y = TeX("Observed $-log_{10}(P)$"))
  
  # subtitle
  if(pos_legend[1] == 'none') {
    lambda = qq_inflation(pvals, df = df, ...)
    p = p + labs(subtitle = TeX(glue("$\\lambda = $", round(lambda, 2))))
  }

  str_lambda = sapply(vals_lambda, function(x) formatC(x, format = 'f', digits = 3))
  if(num_test == 1) {
    labs_col = lapply(sprintf(r'($\lambda = $%s)', str_lambda), TeX)
    cols = 'black'
    p = p + scale_color_manual(values = cols) + theme(legend.position = 'none') 

    if(full_stats) {
      str_stats = glue(
        'lambda = {round(vals_lambda, 2)}', '\n',
        '# Total = ', length(pvals_input) %>% fun_format_big, '\n',
        '# Displayed = ', n_dat %>% fun_format_big, '\n',
        '# NAs = ', sum(is.na(pvals_input)), '\n',
        '# Inf = ', sum(is.infinite(pvals_input)))
      p = p + annotation_compass(str_stats, 'NW')
    } else {
      str_stats = TeX(sprintf(r'($\lambda = $%s)', str_lambda))
      p = p + annotation_compass(str_stats, 'NW')
    }
  } else {
    labs_col = lapply(sprintf(r'(%s: $\lambda = $%s)', tests, str_lambda), TeX)
    cols = viridis(num_test)
    p = p + scale_color_manual(values = cols, labels = labs_col)
    p = p + theme(legend.position = pos_legend) + labs(color = NULL)
  }

  return(p)
}
