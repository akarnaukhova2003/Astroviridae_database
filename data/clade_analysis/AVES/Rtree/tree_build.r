library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(randomcoloR)
library(phytools)

plot_annotated_tree_clusters <- function(tree_file, meta_file, taxlabel=TRUE){
  
  tree <- read.tree(tree_file)
  tree$edge.length[is.na(tree$edge.length)] <- 0
  tree$tip.label <- sapply(strsplit(tree$tip.label, "/"), `[`, 1)
  tree_rooted <- midpoint.root(tree)
  
  info <- read.csv(meta_file, sep=",", stringsAsFactors = FALSE, check.names = FALSE)
  colnames(info) <- gsub(" ", "", colnames(info))
  rownames(info) <- info$id
  
  info$aa83_1B[is.na(info$aa83_1B)] <- "NA"
  info$aa90_1B[is.na(info$aa90_1B)] <- "NA"
  info$nt86_1B[is.na(info$nt86_1B)] <- "NA"
  info$nt80_1B[is.na(info$nt80_1B)] <- "NA"
  info$Host[is.na(info$Host)] <- "NA"
  
  cluster1 <- data.frame(AA83 = factor(info$aa83_1B), row.names = info$id)
  cluster2 <- data.frame(AA90 = factor(info$aa90_1B), row.names = info$id)
  cluster3 <- data.frame(NT86 = factor(info$nt86_1B), row.names = info$id)
  cluster4 <- data.frame(NT80 = factor(info$nt80_1B), row.names = info$id)
  host     <- data.frame(Host = factor(info$Host), row.names = info$id)
  
  colors_f1 <- distinctColorPalette(length(unique(cluster1$AA83)))
  colors_f2 <- distinctColorPalette(length(unique(cluster2$AA90)))
  colors_f3 <- distinctColorPalette(length(unique(cluster3$NT86)))
  colors_f4 <- distinctColorPalette(length(unique(cluster4$NT80)))
  colors_host <- distinctColorPalette(length(unique(host$Host)))
  
  if (taxlabel){
    p <- ggtree(tree_rooted, size=0.7) + geom_tiplab(size=2) + geom_treescale()
  } else {
    p <- ggtree(tree_rooted, size=0.7) + geom_treescale()
  }
  
  label_theme <- theme(axis.text.x = element_text(angle=45, hjust=0, vjust=0, size=0.5))
  
  p1 <- gheatmap(p, cluster1, width=0.05, offset=0.12,
                 colnames=TRUE, font.size=2, colnames_position="top") +
    scale_fill_manual(values=colors_f1, guide="none") + label_theme + new_scale_fill()
  
  p1 <- gheatmap(p1, cluster2, width=0.05, offset=0.24,
                 colnames=TRUE, font.size=2, colnames_position="top") +
    scale_fill_manual(values=colors_f2, guide="none") + label_theme + new_scale_fill()
  
  p1 <- gheatmap(p1, cluster3, width=0.05, offset=0.36,
                 colnames=TRUE, font.size=2, colnames_position="top") +
    scale_fill_manual(values=colors_f3, guide="none") + label_theme + new_scale_fill()
  
  p1 <- gheatmap(p1, cluster4, width=0.05, offset=0.48,
                 colnames=TRUE, font.size=2, colnames_position="top") +
    scale_fill_manual(values=colors_f4, guide="none") + label_theme + new_scale_fill()
  
  p1 <- gheatmap(p1, host, width=0.05, offset=0.60,
                 colnames=TRUE, font.size=2, colnames_position="top") +
    scale_fill_manual(values=colors_host, name="Host") + label_theme
  
  return(p1)
}

tree_file <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Aves_full_seq_2_aa_align_new.nwk"
meta_file <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Astroviridae_Aves_05042026_clusters.csv"

p <- plot_annotated_tree_clusters(tree_file, meta_file)
print(p)
ggsave("/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/2_tree.png",
       p, width=12, height=8, dpi=300)

