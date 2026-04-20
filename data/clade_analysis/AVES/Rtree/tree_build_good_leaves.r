library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(phytools)
library(randomcoloR)

meta_path <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Astroviridae_Aves_18042026_clusters.csv"

all_meta <- read.csv(meta_path, sep = ";", stringsAsFactors = FALSE, check.names = FALSE)
all_meta
fix_na <- function(x) {
  x[is.na(x)] <- "NA"
  x
}

all_meta$aa83_1B <- fix_na(all_meta$aa83_1B)
all_meta$aa90_1B <- fix_na(all_meta$aa90_1B)
all_meta$nt86_1B <- fix_na(all_meta$nt86_1B)
all_meta$nt80_1B <- fix_na(all_meta$nt80_1B)
all_meta$Host <- fix_na(all_meta$Host_fam_ru)

levels_f1 <- sort(unique(all_meta$aa83_1B))
levels_f2 <- sort(unique(all_meta$aa90_1B))
levels_f3 <- sort(unique(all_meta$nt86_1B))
levels_f4 <- sort(unique(all_meta$nt80_1B))
levels_host <- sort(unique(all_meta$Host_fam_ru))

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
  
  info <- read.csv(meta_file,sep = ";", stringsAsFactors = FALSE, check.names = FALSE)
  colnames(info) <- gsub(" ", "", colnames(info))
  print(colnames(info))
  rownames(info) <- info$id
  
  fix_na <- function(x) {
    x[is.na(x)] <- "NA"
    x
  }
  
  info$aa83_1B <- fix_na(info$aa83_1B)
  info$aa90_1B <- fix_na(info$aa90_1B)
  info$nt86_1B <- fix_na(info$nt86_1B)
  info$nt80_1B <- fix_na(info$nt80_1B)
  info$Host <- fix_na(info$Host_fam_ru)
  info$sequence_nm[is.na(info$sequence_nm)] <- "unknown"
  
  cluster1 <- data.frame(AA83 = factor(info$aa83_1B, levels = levels_f1), row.names = info$id)
  cluster2 <- data.frame(AA90 = factor(info$aa90_1B, levels = levels_f2), row.names = info$id)
  cluster3 <- data.frame(NT86 = factor(info$nt86_1B, levels = levels_f3), row.names = info$id)
  cluster4 <- data.frame(NT80 = factor(info$nt80_1B, levels = levels_f4), row.names = info$id)
  host <- data.frame(Host = factor(info$Host, levels = levels_host), row.names = info$id)
  
  p <- ggtree(tree_rooted, size = 0.7)
  
  df <- p$data
  df$bootstrap <- suppressWarnings(as.numeric(df$label))
  p$data <- df
  
  if (taxlabel) {
    tip_labels_map <- setNames(info$sequence_nm, info$id)
    
    p <- p + 
      geom_tiplab(
        aes(label = tip_labels_map[.data$label]),
        offset = 0.02,
        size = 1.2
      ) +
      geom_point2(
        aes(subset = !isTip & !is.na(bootstrap) & bootstrap >= 0.95),
        shape = 21,
        size = 2,
        fill = "hotpink",
        color = "black"
      ) +
      new_scale_fill()
  } else {
    p <- p + geom_treescale()
  }
  
  theme_small_x <- theme(axis.text.x = element_blank())
  
  p1 <- gheatmap(p, cluster1, width = 0.05, offset = 0.12,
                 colnames = TRUE, font.size = 2, colnames_position = "top") +
    scale_fill_manual(values = color_map_f1, guide = "none", na.value = "grey90") +
    theme_small_x + new_scale_fill()
  
  p1 <- gheatmap(p1, cluster2, width = 0.05, offset = 0.24,
                 colnames = TRUE, font.size = 2, colnames_position = "top") +
    scale_fill_manual(values = color_map_f2, guide = "none", na.value = "grey90") +
    theme_small_x + new_scale_fill()
  
  p1 <- gheatmap(p1, cluster3, width = 0.05, offset = 0.36,
                 colnames = TRUE, font.size = 2, colnames_position = "top") +
    scale_fill_manual(values = color_map_f3, guide = "none", na.value = "grey90") +
    theme_small_x + new_scale_fill()
  
  p1 <- gheatmap(p1, cluster4, width = 0.05, offset = 0.48,
                 colnames = TRUE, font.size = 2, colnames_position = "top") +
    scale_fill_manual(values = color_map_f4, guide = "none", na.value = "grey90") +
    theme_small_x + new_scale_fill()
  
  p1 <- gheatmap(p1, host, width = 0.05, offset = 0.60,
                 colnames = TRUE, font.size = 2, colnames_position = "top") +
    scale_fill_manual(values = color_map_host, name = "Host", na.value = "grey90") +
    theme_small_x
  
  return(p1)
}

tree_file <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Aves_1B_97_nt_tree.nwk"

p <- plot_annotated_tree_clusters(tree_file, meta_path)
print(p)

ggsave(
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/nt_1B_97_tree.png",
  p,
  width = 14,
  height = 10,
  dpi = 300
)

