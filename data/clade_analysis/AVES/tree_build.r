library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(randomcoloR)
library(ggplot2)
library(ggtree)

plot_annotated_tree_clusters = function(tree_file, meta_file, taxlabel=TRUE){
  tree = read.tree(tree_file)
  tree$tip.label <- sapply(strsplit(tree$tip.label, "/"), `[`, 1)
  tree_rooted = tree
  info = read.csv(meta_file, sep=",", stringsAsFactors = FALSE, check.names = FALSE)
  colnames(info) <- gsub(" ", "", colnames(info))
  rownames(info) <- info$id
  info$feat1[is.na(info$feat1)] <- "NA"
  info$feat2[is.na(info$feat2)] <- "NA"
  info$feat3[is.na(info$feat3)] <- "NA"
  info$feat4[is.na(info$feat4)] <- "NA"
  info$Host[is.na(info$Host)] <- "NA"
  cluster1 = data.frame(f1 = as.factor(info$feat1))
  cluster2 = data.frame(f2 = as.factor(info$feat2))
  cluster3 = data.frame(f3 = as.factor(info$feat3))
  cluster4 = data.frame(f4 = as.factor(info$feat4))
  rownames(cluster1) <- info$id
  rownames(cluster2) <- info$id
  rownames(cluster3) <- info$id
  rownames(cluster4) <- info$id
  host = data.frame(host = info$Host)
  rownames(host) <- info$id
  colors_f1 = distinctColorPalette(length(unique(cluster1$f1)))
  colors_f2 = distinctColorPalette(length(unique(cluster2$f2)))
  colors_f3 = distinctColorPalette(length(unique(cluster3$f3)))
  colors_f4 = distinctColorPalette(length(unique(cluster4$f4)))
  colors_host = distinctColorPalette(length(unique(host$host)))
  if (taxlabel){
    p = ggtree(tree_rooted, size=0.7) + geom_tiplab(size=2) + geom_treescale()
  } else {
    p = ggtree(tree_rooted, size=0.7) + geom_treescale()
  }
  p1 = gheatmap(p, cluster1, width=0.05, offset=0.1) + scale_fill_manual(values=colors_f1)
  p1 = p1 + new_scale_fill()
  p1 = gheatmap(p1, cluster2, width=0.05, offset=0.2) + scale_fill_manual(values=colors_f2)
  p1 = p1 + new_scale_fill()
  p1 = gheatmap(p1, cluster3, width=0.05, offset=0.3) + scale_fill_manual(values=colors_f3)
  p1 = p1 + new_scale_fill()
  p1 = gheatmap(p1, cluster4, width=0.05, offset=0.4) + scale_fill_manual(values=colors_f4)
  p1 = p1 + new_scale_fill()
  p1 = gheatmap(p1, host, width=0.05, offset=0.5) + scale_fill_manual(values=colors_host)
  return(p1)
}



tree_file <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Aves_full_seq_1A_aa_align_new.nwk"
meta_file <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Astroviridae_Aves_28032026_clusters.csv"

p <- plot_annotated_tree_clusters(tree_file, meta_file)
print(p)
ggsave("tree_annotated.png", p, width=12, height=8, dpi=300)

install.packages("randomcoloR")
library(randomcoloR)
