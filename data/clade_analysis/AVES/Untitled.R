library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(phytools)
library(randomcoloR)

meta_path <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Astroviridae_Aves_11042026_clusters.csv"

all_meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)

fix_na <- function(x) {
  x[is.na(x)] <- "NA"
  x
}

all_meta$aa83_1B <- fix_na(all_meta$aa83_1B)
all_meta$aa90_1B <- fix_na(all_meta$aa90_1B)
all_meta$nt86_1B <- fix_na(all_meta$nt86_1B)
all_meta$nt80_1B <- fix_na(all_meta$nt80_1B)
all_meta$Host <- fix_na(all_meta$Host)

levels_f1 <- sort(unique(all_meta$aa83_1B))
levels_f2 <- sort(unique(all_meta$aa90_1B))
levels_f3 <- sort(unique(all_meta$nt86_1B))
levels_f4 <- sort(unique(all_meta$nt80_1B))
levels_host <- sort(unique(all_meta$Host))

set.seed(123)

make_map <- function(levels) {
  cols <- distinctColorPalette(length(levels))
  setNames(cols, levels)
}

color_map_f1 <- make_map(levels_f1)
color_map_f2 <- make_map(levels_f2)
color_map_f3 <- make_map(levels_f3)
color_map_f4 <- make_map(levels_f4)
color_map_host <- make_map(levels_host)

plot_annotated_tree_clusters <- function(tree_file, meta_file, taxlabel = TRUE) {
  
  tree <- read.tree(tree_file)
  tree$edge.length[is.na(tree$edge.length)] <- 0
  tree$tip.label <- sapply(strsplit(tree$tip.label, "/"), `[`, 1)
  
  tree_rooted <- midpoint.root(tree)
  tree_rooted$node.label <- as.numeric(tree_rooted$node.label)
  
  info <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(info) <- gsub(" ", "", colnames(info))
  rownames(info) <- info$id
  
  info$aa83_1B <- fix_na(info$aa83_1B)
  info$aa90_1B <- fix_na(info$aa90_1B)
  info$nt86_1B <- fix_na(info$nt86_1B)
  info$nt80_1B <- fix_na(info$nt80_1B)
  info$Host <- fix_na(info$Host)
  
  info$sequence_nm[is.na(info$sequence_nm)] <- "unknown"
  
  cluster1 <- data.frame(AA83 = factor(info$aa83_1B, levels = levels_f1),
                         row.names = info$id)
  
  cluster2 <- data.frame(AA90 = factor(info$aa90_1B, levels = levels_f2),
                         row.names = info$id)
  
  cluster3 <- data.frame(NT86 = factor(info$nt86_1B, levels = levels_f3),
                         row.names = info$id)
  
  cluster4 <- data.frame(NT80 = factor(info$nt80_1B, levels = levels_f4),
                         row.names = info$id)
  
  host <- data.frame(Host = factor(info$Host, levels = levels_host),
                     row.names = info$id)
  
  p <- ggtree(tree_rooted, size = 0.7)
  
  if (taxlabel) {
    tip_labels_map <- setNames(info$sequence_nm, info$id)
    
    p <- p + geom_tiplab(
      aes(label = tip_labels_map[.data$label]),
      offset = 0.02,
      size = 1.2
    )
  } else {
    p <- p + geom_treescale()
  }
  df <- p$data
  
  df$node <- as.numeric(df$node)
  
  boot_df <- df[!df$isTip, c("node", "label")]
  boot_df$bootstrap <- as.numeric(boot_df$label)
  
  edges <- as.data.frame(tree_rooted$edge)
  colnames(edges) <- c("parent", "child")
  
  edges_boot <- merge(
    edges,
    boot_df,
    by.x = "parent",
    by.y = "node",
    all.x = TRUE
  )
  
  high_parents <- edges_boot$parent[edges_boot$bootstrap > 90]
  
  high_children <- edges_boot$child[edges_boot$parent %in% high_parents]
  
  df$high_boot <- df$node %in% high_children
  
  p$data <- df
  
  p <- p + geom_tippoint(
    aes(subset = high_boot),
    shape = 21,
    size = 2,
    fill = "black"
  )
  
  theme_small_x <- theme(axis.text.x = element_blank())
  
  p1 <- gheatmap(p, cluster1, width = 0.05, offset = 0.12, colnames = FALSE) +
    scale_fill_manual(values = color_map_f1, guide = "none", na.value = "grey90") +
    new_scale_fill()
  
  p1 <- gheatmap(p1, cluster2, width = 0.05, offset = 0.24, colnames = FALSE) +
    scale_fill_manual(values = color_map_f2, guide = "none", na.value = "grey90") +
    new_scale_fill()
  
  p1 <- gheatmap(p1, cluster3, width = 0.05, offset = 0.36, colnames = FALSE) +
    scale_fill_manual(values = color_map_f3, guide = "none", na.value = "grey90") +
    new_scale_fill()
  
  p1 <- gheatmap(p1, cluster4, width = 0.05, offset = 0.48, colnames = FALSE) +
    scale_fill_manual(values = color_map_f4, guide = "none", na.value = "grey90") +
    new_scale_fill()
  
  p1 <- gheatmap(p1, host, width = 0.05, offset = 0.60, colnames = FALSE) +
    scale_fill_manual(values = color_map_host, name = "Host", na.value = "grey90") +
    theme_small_x
  
  return(p1)
}

tree_file <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Aves_2_97_translated_tree.nwk"

p <- plot_annotated_tree_clusters(tree_file, meta_path)
print(p)

ggsave(
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/2_97_tree.png",
  p,
  width = 14,
  height = 10,
  dpi = 300
)

