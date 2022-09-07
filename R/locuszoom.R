qq_plot_ld <- function(corr, 
  split = FALSE,
  thr_r2 = 0.01, min_size = 10, max_size = 50, max_K = 50,
  verbose = 0)
{
  require(bigsnpr)
  require(cowplot)
  require(ggsci)

  if(split) {
    out_split <- snp_ldsplit(corr, thr_r2 = thr_r2, min_size = min_size, max_size = min_size, max_K = max_K)
    all_ind <- head(out_split$all_last[[6]], -1)
  }

  # transform sparse representation into (i,j,x) triplets
  corrT <- as(corr, "dgTMatrix")
  upper <- (corrT@i <= corrT@j & corrT@x^2 >= thr_r2)
  tab <- tibble(
    i = corrT@i[upper] + 1L,
    j = corrT@j[upper] + 1L,
    r2 = corrT@x[upper]^2) %>% 
  mutate(y = (j - i)/2)

  r2_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  ptab <- tab %>%
    mutate(r2_cat = cut(r2, breaks = r2_breaks)) %>%
    arrange(r2_cat)

  col_locuszoom <- ggsci::pal_locuszoom()(5)
  p <- ggplot(ptab) +
    geom_point(aes(i + y, y, color = r2_cat), 
      shape = 18, size = rel(0.45), alpha = 1) +
    coord_fixed() +
    scale_color_manual(values = rev(col_locuszoom), guide = "none")

  if(split) {
    p <- p + geom_vline(xintercept = all_ind + 0.5, linetype = 3, size = rel(0.1))
  }

  p
}

