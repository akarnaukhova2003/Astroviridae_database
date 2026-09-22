library(ape)
library(ggtree)
library(ggplot2)
library(randomcoloR)
library(colorspace)
library(dplyr)

setwd(
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/snakemake/results/trees/"
)

reference_tree_file = "P_Y/P_Y_ORF1b.treefile"
metadata_file = "Astroviridae_Aves_Amphibia_14072026.csv"

outgroup_id = "MG599917"

tree_files = c(
  "P_Y/P_Y_ORF1a.treefile",
  "P_Y/P_Y_ORF1b.treefile",
  "P_Y/P_Y_ORF2_1.treefile",
  "P_Y/P_Y_ORF2_2.treefile"
)

clades = c(51, 76, 49, 63)

normalize_id = function(x) {
  x = trimws(x)
  x = gsub("'", "", x)
  x = sub("/.*$", "", x)
  x = gsub("_", "-", x)
  x
}

metadata = read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata$ID = normalize_id(metadata$ID)

tree = read.tree(reference_tree_file)

tree$tip.label = normalize_id(tree$tip.label)

tree = root(
  tree,
  outgroup = outgroup_id,
  resolve.root = TRUE
)

p = ggtree(tree)

clade_taxa = lapply(
  clades,
  function(x) normalize_id(get_taxa_name(p, x))
)

base_colors = distinctColorPalette(length(clades))

color_table = data.frame()

for (i in seq_along(clades)) {
  
  ids = clade_taxa[[i]]
  
  colors = colorRampPalette(
    c(
      lighten(base_colors[i], 0.4),
      darken(base_colors[i], 0.4)
    )
  )(length(ids))
  
  color_table = rbind(
    color_table,
    data.frame(
      ID = ids,
      color = colors
    )
  )
}

metadata = left_join(
  metadata,
  color_table,
  by = "ID"
)

write.csv(
  metadata,
  "metadata_upd.csv",
  row.names = FALSE
)

for (file in tree_files) {
  
  tree = read.tree(file)
  
  tree$tip.label = normalize_id(tree$tip.label)
  
  tree = root(
    tree,
    outgroup = outgroup_id,
    resolve.root = TRUE
  )
  
  p = ggtree(tree)
  
  p$data$ID = normalize_id(p$data$label)
  
  p$data = left_join(
    p$data,
    color_table,
    by = "ID"
  )
  
  p$data$color[is.na(p$data$color)] = "black"
  
  p = p +
    geom_tiplab(
      aes(color = color),
      size = 4
    ) +
    scale_color_identity()
  
  bootstrap_data = p$data[
    !p$data$isTip &
      !is.na(p$data$label) &
      suppressWarnings(
        as.numeric(p$data$label)
      ) > 90,
  ]
  
  p = p +
    geom_nodepoint(
      data = bootstrap_data,
      size = 3
    ) +
    geom_treescale()
  
  name = tools::file_path_sans_ext(
    basename(file)
  )
  
  ggsave(
    paste0(name, "_gradient.png"),
    p,
    width = 20,
    height = 15,
    dpi = 300
  )
  
  ggsave(
    paste0(name, "_gradient.pdf"),
    p,
    width = 20,
    height = 15
  )
  
  ggsave(
    paste0(name, "_gradient.svg"),
    p,
    width = 20,
    height = 15
  )
}

cat("Все деревья обработаны.\n")