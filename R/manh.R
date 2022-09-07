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
qq_manh <- function(tab)
{
  tab <- tab20
  
  tab <- rename(tab, snp = MarkerName, chr = Chr, pos = Pos, pval = P.value)
  
  tab <- select(tab, snp, chr, pos, pval) %>%
    arrange(chr, pos) %>% mutate(ord = seq(1, n()))

  ### compute breaks on x axis
  xbreaks_side <- group_by(tab, chr) %>% 
    filter(ord %in% ord[c(1, length(ord))]) %>% ungroup %>% 
    select(ord) %>% mutate(chr = "")
    
  xbreaks_mid <- group_by(tab, chr) %>% 
    filter(ord %in% ord[ceiling(length(ord) / 2)]) %>% ungroup %>% 
    select(ord, chr) %>% mutate(chr = as.character(chr))
    
  xbreaks <- bind_rows(xbreaks_side, xbreaks_mid) %>% arrange(ord)

  ### plot
  p <- ggplot(tab, aes(ord, -log10(pval))) + geom_point()
  
  p <- p + scale_x_continuous("Chromosome", breaks = xbreaks$ord, labels = xbreaks$chr)

  return(p)
}
