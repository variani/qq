#' @export
manhattan_plot <- function(tab, p = "pval", snp = "predictor", ...)
{
  # assign column names as required by qqman
  tab$p <- tab[[p]]
  tab$snp <- tab[[snp]]
  tab <- tab %>% select(p, snp)

  # extract chr:pos from snp name
  info <- tab$snp %>% strsplit(":") %>% do.call(rbind, .) %>% .[, 1:2] %>% as_tibble
  names(info) <- c("chr", "pos")
  info <- mutate(info, chr = as.integer(chr), pos = as.integer(pos))

  tab <- bind_cols(tab, info)

  # plot with qqman
  qqman::manhattan(tab, chr = "chr", bp = "pos", p = "p", ...)
}

#' Manhatten plots
#'
#' @examples
#' tab <- read_delim("~/Data/SCC/GWAS_processed/gwas/hm3.Consortium_SCC_Meta_SNPonly_impu_filter_anno.txt", delim = " ")
#' tab22 <- filter(tab, Chr == 22)
#'
#' @export
qq_manh <- function(tab, cols = c(snp = 'snp', chr = 'chr', pos = 'pos', pval = 'pval'),
  p_thr = 5e-8, lab_chr = NULL)
{
  stopifnot(require(latex2exp))

  ## check old names from cols
  for(col in cols) {
    if(!(col %in% colnames(tab))) {
      warning(paste('column', col, ' is not present'))
      return(invisible())
    }
  }
  tab <- rename(tab, cols)

  ## rename columns
  cols_new = c('snp', 'chr', 'pos', 'pval')
  for(col in cols_new) {
    if(!(col %in% colnames(tab))) {
      warning(paste('column', col, ' is required; use "cols" argument to specify replacement'))
      return(invisible())
    }
  }
  tab <- select(tab, snp, chr, pos, pval)

  ## convert chr/pos to integers
  tab = mutate_at(tab, c('chr', 'pos'), as.integer)
  stopifnot(!any(is.na(tab$chr))); stopifnot(!any(is.na(tab$pos)))

  tab = mutate(tab, chr_even = factor(chr %% 2))

  ## order SNPs
  tab = arrange(tab, chr, pos) %>% mutate(ord = seq(1, n()))

  ### compute breaks on x axis
  xbreaks_side <- group_by(tab, chr) %>% 
    filter(ord %in% ord[c(1, length(ord))]) %>% ungroup %>% 
    select(ord) %>% mutate(chr = "")
    
  xbreaks_mid <- group_by(tab, chr) %>% 
    filter(ord %in% ord[ceiling(length(ord) / 2)]) %>% ungroup %>% 
    select(ord, chr) %>% mutate(chr = as.character(chr))
    
  xbreaks <- bind_rows(xbreaks_side, xbreaks_mid) %>% arrange(ord)
  xbreaks <- filter(xbreaks, chr != "")
  if(!is.null(lab_chr)) {
    xbreaks <- mutate(xbreaks, chr = ifelse(chr %in% lab_chr, chr, ""))
  }

  ### plot
  p <- ggplot(tab, aes(ord, -log10(pval))) + geom_point(aes(color = chr_even))
  p <- p + geom_hline(yintercept = -log10(p_thr), color = 2)
  p <- p + scale_x_continuous(breaks = xbreaks$ord, labels = xbreaks$chr)
  p <- p + labs(x = 'Chromosome', y = TeX("$-log_{10}(P)$"))
  p = p + scale_color_manual(values = c(8, 1), guide = 'none')

  return(p)
}
