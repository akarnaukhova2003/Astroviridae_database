library(ape)
library(ggtree)
library(ggplot2)
library(phytools)
library(randomcoloR)
library(colorspace)
library(dplyr)

setwd(
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/Aves_Amphibia_trees/gradient_tree"
)

reference_tree_file = "B_clade/B_ORF1a.nwk"

metadata_file = "Astroviridae_Aves_Amphibia_14072026.csv"

tree_files = c(
  "B_clade/B_ORF1a.nwk",
  "B_clade/B_ORF1b.nwk",
  "B_clade/B_ORF2_trim.nwk"
)

clades = c(
  95,
  97, 
  98,
  101,
  58,
  109,
  69,
  75,
  63
  
)

reference_tree = read.nexus(
  reference_tree_file
)

reference_tree = midpoint.root(
  reference_tree
)

reference_plot = ggtree(
  reference_tree,
  size = 0.75
)

order_list = list()

for (node in clades) {
  
  taxa = get_taxa_name(
    reference_plot,
    node
  )
  
  taxa = gsub(
    "'",
    "",
    taxa
  )
  
  taxa = sub(
    "/.*$",
    "",
    taxa
  )
  
  taxa = gsub(
    "_",
    "-",
    taxa
  )
  
  order_list[[as.character(node)]] = taxa
}

base_colors = distinctColorPalette(
  length(order_list)
)

color_tables = list()

for (i in seq_along(order_list)) {
  
  taxa = order_list[[i]]
  
  base = base_colors[i]
  
  color_light = lighten(
    base,
    0.4
  )
  
  color_dark = darken(
    base,
    0.4
  )
  
  color_function = colorRampPalette(
    c(
      color_light,
      color_dark
    )
  )
  
  colors_clade = color_function(
    length(taxa)
  )
  
  color_tables[[i]] = data.frame(
    ID = taxa,
    color = colors_clade,
    stringsAsFactors = FALSE
  )
}

color_table = bind_rows(
  color_tables
) %>%
  distinct(
    ID,
    .keep_all = TRUE
  )

info = read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

info$ID = gsub(
  "'",
  "",
  info$ID
)

info$ID = sub(
  "/.*$",
  "",
  info$ID
)

info$ID = gsub(
  "_",
  "-",
  info$ID
)

info_upd = info %>%
  left_join(
    color_table,
    by = "ID"
  )

write.csv(
  info_upd,
  "metadata_upd.csv",
  row.names = FALSE
)

plot_gradient_tree = function(
    tree_file,
    metadata
) {
  
  tree = read.nexus(
    tree_file
  )
  
  tree = midpoint.root(
    tree
  )
  
  full_labels = gsub(
    "'",
    "",
    tree$tip.label
  )
  
  tree_ids = sub(
    "/.*$",
    "",
    full_labels
  )
  
  tree_ids = gsub(
    "_",
    "-",
    tree_ids
  )
  
  p = ggtree(
    tree,
    size = 0.75
  )
  
  p$data$ID = NA_character_
  
  p$data$ID[
    p$data$isTip
  ] = tree_ids
  
  p$data$full_label = NA_character_
  
  p$data$full_label[
    p$data$isTip
  ] = full_labels
  
  color_info = metadata %>%
    select(
      ID,
      color
    ) %>%
    filter(
      !is.na(color)
    ) %>%
    distinct(
      ID,
      .keep_all = TRUE
    )
  
  p$data = left_join(
    p$data,
    color_info,
    by = "ID"
  )
  
  tips_colored = p$data %>%
    filter(
      isTip,
      !is.na(color)
    )
  
  tips_uncolored = p$data %>%
    filter(
      isTip,
      is.na(color)
    )
  
  p =
    p +
    geom_tiplab(
      data = tips_uncolored,
      aes(
        label = full_label
      ),
      color = "black",
      size = 6,
      hjust = 0
    ) +
    geom_tiplab(
      data = tips_colored,
      aes(
        label = full_label,
        color = color
      ),
      size = 6,
      hjust = 0
    ) +
    scale_color_identity() +
    geom_treescale(
      fontsize = 6
    ) +
    theme(
      legend.position = "none",
      plot.margin = margin(
        10,
        300,
        10,
        10
      ),
      plot.title = element_text(
        size = 28
      )
    )
  
  return(p)
}

for (tree_file in tree_files) {
  
  p = plot_gradient_tree(
    tree_file,
    info_upd
  )
  
  tree_name = tools::file_path_sans_ext(
    basename(tree_file)
  )
  
  p = p +
    ggtitle(
      tree_name
    )
  
  ggsave(
    filename = paste0(
      tree_name,
      "_gradient.png"
    ),
    plot = p,
    width = 45,
    height = 30,
    dpi = 300
  )
  
  ggsave(
    filename = paste0(
      tree_name,
      "_gradient.pdf"
    ),
    plot = p,
    width = 45,
    height = 30
  )
  
  ggsave(
    filename = paste0(
      tree_name,
      "_gradient.svg"
    ),
    plot = p,
    width = 45,
    height = 30
  )
  
  cat(
    tree_name,
    ":",
    sum(p$data$isTip),
    "листьев;",
    sum(
      p$data$isTip &
        !is.na(p$data$color)
    ),
    "окрашенных\n"
  )
}