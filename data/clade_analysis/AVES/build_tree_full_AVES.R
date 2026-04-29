library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(phytools)
library(randomcoloR)

meta_path <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Astroviridae_Aves_26042026_clusters_shorter.csv"

all_meta <- read.csv(meta_path, sep = ";", stringsAsFactors = FALSE, check.names = FALSE)

fix_na <- function(x) {
  x[is.na(x)] <- "NA"
  x
}

all_meta$aa83_1B <- fix_na(all_meta$aa83_1B)
all_meta$aa90_1B <- fix_na(all_meta$aa90_1B)
all_meta$nt86_1B <- fix_na(all_meta$nt86_1B)
all_meta$nt80_1B <- fix_na(all_meta$nt80_1B)
all_meta$Host_shorter <- fix_na(trimws(all_meta$Host_shorter))

levels_f1 <- sort(unique(all_meta$aa83_1B))
levels_f2 <- sort(unique(all_meta$aa90_1B))
levels_f3 <- sort(unique(all_meta$nt86_1B))
levels_f4 <- sort(unique(all_meta$nt80_1B))

set.seed(123)

make_map <- function(levels) {
  cols <- distinctColorPalette(length(levels))
  setNames(cols, levels)
}

color_map_f1 <- make_map(levels_f1)
color_map_f2 <- make_map(levels_f2)
color_map_f3 <- make_map(levels_f3)
color_map_f4 <- make_map(levels_f4)

bird_colors <- c(
  "Anatidae (утиные)" = "#F4D35E",
  "Phasianidae (фазановые)" = "#800080",
  "Columbidae (голубиные)" = "#71A9F7",
  "Corvidae (врановые)" = "#380036",
  "Passeridae (воробьиные)" = "#0CBABA",
  "Scolopacidae (бекасовые)" = "#C7F2A7",
  "Paridae (синицевые)" = "#D9D0DE",
  "Prunellidae (завирушковые)" = "#800000",
  "Petroicidae (австралийские малиновки)" = "#F40076",
  "Recurvirostridae (шилоклювковые)" = "#FFA630",
  "Cacatuidae (какаду)" = "#380036",
  "Psittacidae (попугаевые)" = "#1C6E8C",
  "Gruidae (журавлиные)" = "#FFAFF0",
  "Acanthizidae (шипоклювковые)" = "#88B7B5",
  "Rhipiduridae (веерохвостковые)" = "#4B0082"
)

manual_host_colors <- c(
  "Mammalia (Млекопитающие)" = "#F8333C",
  "Reptilia (Пресмыкающиеся)" = "#0A2472",
  "Amphibia (Амфибии)" = "#06BA63",
  bird_colors
)

host_order <- c(
  "Mammalia (Млекопитающие)",
  "Reptilia (Пресмыкающиеся)",
  "Amphibia (Амфибии)",
  names(bird_colors),
  "Other"
)

plot_annotated_tree_clusters <- function(tree_file, meta_file) {
  
  tree <- read.tree(tree_file)
  tree$edge.length[is.na(tree$edge.length)] <- 0
  
  # 🔥 чистим ID
  tree$tip.label <- trimws(sapply(strsplit(tree$tip.label, "/"), `[`, 1))
  
  tree_rooted <- midpoint.root(tree)
  
  info <- read.csv(meta_file, sep = ";", stringsAsFactors = FALSE, check.names = FALSE)
  colnames(info) <- gsub(" ", "", colnames(info))
  
  info$id <- trimws(info$id)
  
  fix_na <- function(x) { x[is.na(x)] <- "NA"; x }
  
  info$aa83_1B <- fix_na(info$aa83_1B)
  info$aa90_1B <- fix_na(info$aa90_1B)
  info$nt86_1B <- fix_na(info$nt86_1B)
  info$nt80_1B <- fix_na(info$nt80_1B)
  
  info$Host <- trimws(info$Host_shorter)
  info$Host <- gsub("\\s+", " ", info$Host)
  info$Host <- fix_na(info$Host)
  
  info$sequence_nm[is.na(info$sequence_nm)] <- "unknown"
  
  
  highlight_ids <- c("NC_002470","Y15936", "PP623814", "NC_005790","EU143845", "AF206663", "AB033998")
  
  info$label_color <- ifelse(info$id %in% highlight_ids, info$id, "other")
  
  highlight_colors <- c(
    "NC_002470" = "#FF0000",
    "Y15936" = "#FF0000",
    "PP623814" = "#00BFC4",
    "NC_005790" = "#00FF00",
    "EU143845" = "#00FF00",
    "AF206663" = "#00FF00",
    "AB033998" = "#C77CFF",
    "other" = "black"
  )
  
  
  p <- ggtree(tree_rooted, size = 0.7)
  
  df <- p$data
  
  df <- merge(df,
              info[, c("id", "label_color")],
              by.x = "label",
              by.y = "id",
              all.x = TRUE)
  
  df$label_color[is.na(df$label_color)] <- "other"
  
  p$data <- df
  
  df$fontface <- ifelse(df$label_color == "other", "plain", "bold")
  p$data <- df
  
  p <- p +
    geom_tiplab(
      aes(color = label_color, fontface = fontface),
      size = 2,
      align = FALSE
    ) +
    scale_color_manual(
      values = highlight_colors,
      guide = "none"
    )
  
  # bootstrap
  df <- p$data
  df$bootstrap <- suppressWarnings(as.numeric(df$label))
  p$data <- df
  
  p <- p +
    geom_point2(
      aes(subset = !isTip & !is.na(bootstrap) & bootstrap >= 0.95),
      shape = 21,
      size = 2,
      fill = "pink",
      color = "black"
    ) +
    new_scale_fill()
  
  cluster_list <- list(
    AA83 = data.frame(AA83 = factor(info$aa83_1B, levels = levels_f1), row.names = info$id),
    AA90 = data.frame(AA90 = factor(info$aa90_1B, levels = levels_f2), row.names = info$id),
    NT86 = data.frame(NT86 = factor(info$nt86_1B, levels = levels_f3), row.names = info$id),
    NT80 = data.frame(NT80 = factor(info$nt80_1B, levels = levels_f4), row.names = info$id)
  )
  
  color_maps <- list(color_map_f1, color_map_f2, color_map_f3, color_map_f4)
  
  info$Host[!info$Host %in% names(manual_host_colors)] <- "Other"
  manual_host_colors["Other"] <- "grey85"
  
  host <- data.frame(
    Host = factor(info$Host, levels = host_order),
    row.names = info$id
  )
  
  base_offset <- 0.03
  step <- 0.07
  width <- 0.07
  
  p1 <- p
  
  for (i in seq_along(cluster_list)) {
    p1 <- gheatmap(
      p1,
      cluster_list[[i]],
      width = width,
      offset = base_offset + (i - 1) * step,
      colnames = TRUE,
      font.size = 2,
      colnames_position = "top"
    ) +
      scale_fill_manual(
        values = color_maps[[i]],
        guide = "none",
        na.value = "grey90"
      ) +
      theme(axis.text.x = element_blank()) +
      new_scale_fill()
  }
  
  p1 <- gheatmap(
    p1,
    host,
    width = width,
    offset = base_offset + length(cluster_list) * step,
    colnames = TRUE,
    font.size = 2,
    colnames_position = "top"
  ) +
    scale_fill_manual(
      values = manual_host_colors,
      breaks = host_order,
      name = "Host",
      drop = FALSE,
      na.value = "grey85"
    ) +
    theme(axis.text.x = element_blank())
  
  return(p1)
}

#tree_file <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Aves_1A_removed_97_nt_tree.nwk"
tree_file <- "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Aves_full_seq_removed_1A_align.nwk" 
p <- plot_annotated_tree_clusters(tree_file, meta_path)

print(p)

ggsave(
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/nt_1A_tree_2704.png",
  p,
  width = 14,
  height = 10,
  dpi = 300
)