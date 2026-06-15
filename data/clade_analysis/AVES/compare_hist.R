library(ape)
library(ggplot2)
library(cowplot)
library(dplyr)

compute_distances <- function(aln, pairwise=TRUE){
  
  dna_char <- as.character.DNAbin(aln)
  dna_char[dna_char == "-"] <- NA
  
  aa_slice <- trans(aln)
  aa_char <- as.character.AAbin(aa_slice)
  aa_char[aa_char == "X"] <- NA
  
  dist_nt <- dist.gene(dna_char, method = "percentage", pairwise.deletion = pairwise)
  dist_aa <- dist.gene(aa_char, method = "percentage", pairwise.deletion = pairwise)
  
  dist_nt <- dist2list(dist_nt)
  dist_aa <- dist2list(dist_aa)
  
  colnames(dist_nt)[3] <- "nt_dist"
  colnames(dist_aa)[3] <- "aa_dist"
  
  return(list(nt = dist_nt, aa = dist_aa))
}

plot_hist <- function(df, column, title,
                      xlim_range = NULL,
                      ylim_range = NULL,
                      fill_color = "grey"){
  
  p <- ggplot(df, aes_string(x = column)) +
    geom_histogram(bins = 50, fill = fill_color, color = "black") +
    theme_bw() +
    labs(
      title = title,
      x = column,
      y = "Count"
    )
  
  if (!is.null(xlim_range) || !is.null(ylim_range)) {
    p <- p + coord_cartesian(xlim = xlim_range, ylim = ylim_range)
  }
  
  return(p)
}

aln1 <- read.dna(
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/Aves/Aves_full_seq_removed_1B_align_clear.fasta",
  format = "fasta"
)

aln2 <- read.dna(
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/Aves/Aves_full_seq_removed_1B_align.fasta",
  format = "fasta"
)

d1 <- compute_distances(aln1)
d2 <- compute_distances(aln2)
xlim_nt <- c(0, 0.5)
xlim_aa <- c(0, 0.5)

ylim_nt <- c(0, 1000)
ylim_aa <- c(0, 1000)

p1_nt <- plot_hist(d1$nt, "nt_dist", "new (NT)", xlim_nt, ylim_nt, "#8ECAE6")
p2_nt <- plot_hist(d2$nt, "nt_dist", "old (NT)", xlim_nt, ylim_nt, "#8ECAE6")

p1_aa <- plot_hist(d1$aa, "aa_dist", "new (AA)", xlim_aa, ylim_aa, "#FFB703")
p2_aa <- plot_hist(d2$aa, "aa_dist", "old (AA)", xlim_aa, ylim_aa, "#FFB703")

final_plot <- plot_grid(
  p1_nt, p2_nt,
  p1_aa, p2_aa,
  ncol = 2,
  labels = c("A", "B", "C", "D")
)
print(final_plot)
ggsave(
  filename = "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/Aves/ORF1B_2x2_hist.png",
  plot = final_plot,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

