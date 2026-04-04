library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(randomcoloR)
library(phytools)

plot_annotated_tree_clusters = function(tree_file, meta_file, taxlabel=TRUE){
  tree = read.tree(tree_file)
  tree$edge.length[is.na(tree$edge.length)] <- 0
  tree$tip.label <- sapply(strsplit(tree$tip.label, "/"), `[`, 1)
  tree_rooted = midpoint.root(tree)
  
  info = read.csv(meta_file, sep=",", stringsAsFactors = FALSE, check.names = FALSE)
  colnames(info) <- gsub(" ", "", colnames(info))
  rownames(info) <- info$id
  
  info$aa83_1B[is.na(info$aa83_1B)] <- "NA"
  info$aa90_1B[is.na(info$aa90_1B)] <- "NA"
  info$nt86_1B[is.na(info$nt86_1B)] <- "NA"
  info$nt80_1B[is.na(info$nt80_1B)] <- "NA"
  info$Host[is.na(info$Host)] <- "NA"
  
  cluster1 = data.frame(f1 = as.factor(info$aa83_1B))
  cluster2 = data.frame(f2 = as.factor(info$aa90_1B))
  cluster3 = data.frame(f3 = as.factor(info$nt86_1B))
  cluster4 = data.frame(f4 = as.factor(info$nt80_1B))
  host = data.frame(host = info$Host)
  
  rownames(cluster1) <- info$id
  rownames(cluster2) <- info$id
  rownames(cluster3) <- info$id
  rownames(cluster4) <- info$id
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
  
  p1 = gheatmap(p, cluster1, width=0.05, offset=0.1, colnames=FALSE, show.legend=FALSE) +
    scale_fill_manual(values=colors_f1)
  p1 = p1 + new_scale_fill()
  
  p1 = gheatmap(p1, cluster2, width=0.05, offset=0.2, colnames=FALSE, show.legend=FALSE) +
    scale_fill_manual(values=colors_f2)
  p1 = p1 + new_scale_fill()
  
  p1 = gheatmap(p1, cluster3, width=0.05, offset=0.3, colnames=FALSE, show.legend=FALSE) +
    scale_fill_manual(values=colors_f3)
  p1 = p1 + new_scale_fill()
  
  p1 = gheatmap(p1, cluster4, width=0.05, offset=0.4, colnames=FALSE, show.legend=FALSE) +
    scale_fill_manual(values=colors_f4)
  p1 = p1 + new_scale_fill()
  
  p1 = gheatmap(p1, host, width=0.05, offset=0.5, colnames=FALSE) +
    scale_fill_manual(values=colors_host, name="Host")
  
  return(p1)
}

tree_file <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Aves_full_seq_1A_aa_align_new.nwk"
meta_file <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Astroviridae_Aves_28032026_clusters.csv"

p <- plot_annotated_tree_clusters(tree_file, meta_file)
print(p)
ggsave("/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/tree_annotated.png", p, width=12, height=8, dpi=300)
