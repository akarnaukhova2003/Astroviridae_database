library(ape)
library(ggtree)
library(ggplot2)
library(phytools)
library(randomcoloR)
library(colorspace)
library(dplyr)

setwd(
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/snakemake/results/trees/"
)

reference_tree_file = "P_Y/P_Y_ORF1b.nwk"
metadata_file = "Astroviridae_Aves_Amphibia_14072026.csv"
outgroup_id = "MG599917"

tree_files = c(
  "P_Y/P_Y_ORF1a.nwk",
  "P_Y/P_Y_ORF1b.nwk", 
  "P_Y/P_Y_ORF2_1.nwk", 
  "P_Y/P_Y_ORF2_2.nwk"
  )

clades = c(
  51, 76, 49, 63 
)

normalize_id = function(x) {
  x = gsub("'", "", x)
  x = sub("/.*$", "", x)
  x = gsub("_", "-", x)
  return(x)
}

extract_bootstrap = function(tree_file) {
  
  lines = readLines(
    tree_file,
    warn = FALSE
  )
  
  tree_text = paste(
    lines,
    collapse = ""
  )
  
  bootstrap_values = regmatches(
    tree_text,
    gregexpr(
      "\\[&label=[^\\]]+\\]",
      tree_text,
      perl = TRUE
    )
  )[[1]]
  
  if (length(bootstrap_values) == 0) {
    return(numeric(0))
  }
  
  bootstrap_values = gsub(
    "^\\[&label=",
    "",
    bootstrap_values
  )
  
  bootstrap_values = gsub(
    "\\]$",
    "",
    bootstrap_values
  )
  
  bootstrap_values = as.numeric(
    bootstrap_values
  )
  
  bootstrap_values = bootstrap_values[
    !is.na(bootstrap_values)
  ]
  
  return(bootstrap_values)
}

reference_tree = read.nexus(
  reference_tree_file
)

reference_tree$tip.label = normalize_id(
  reference_tree$tip.label
)

outgroup_index = which(
  reference_tree$tip.label == outgroup_id
)

if (length(outgroup_index) == 0) {
  stop(
    paste0(
      "Аутгруппа не найдена: ",
      outgroup_id
    )
  )
}

reference_tree = root(
  reference_tree,
  outgroup = outgroup_index,
  resolve.root = TRUE
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
  
  order_list[[as.character(node)]] =
    normalize_id(taxa)
}

base_colors = distinctColorPalette(
  length(order_list)
)

color_tables = list()

for (i in seq_along(order_list)) {
  
  taxa = order_list[[i]]
  
  colors_clade = colorRampPalette(
    c(
      lighten(base_colors[i], 0.4),
      darken(base_colors[i], 0.4)
    )
  )(
    max(1, length(taxa))
  )
  
  color_tables[[i]] =
    data.frame(
      ID = taxa,
      color = colors_clade
    )
}

color_table =
  bind_rows(color_tables) %>%
  distinct(
    ID,
    .keep_all = TRUE
  )

info = read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

info$ID = normalize_id(info$ID)

info_upd =
  info %>%
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
  
  tree = read.nexus(tree_file)
  
  tree$tip.label =
    normalize_id(tree$tip.label)
  
  outgroup_index = which(
    tree$tip.label == outgroup_id
  )
  
  if (length(outgroup_index) == 0) {
    stop(
      paste0(
        "Аутгруппа ",
        outgroup_id,
        " не найдена в ",
        tree_file
      )
    )
  }
  
  tree = root(
    tree,
    outgroup = outgroup_index,
    resolve.root = TRUE
  )
  
  bootstrap_values =
    extract_bootstrap(tree_file)
  
  p = ggtree(
    tree,
    size = 0.75
  )
  
  p$data$ID = NA_character_
  p$data$bootstrap = NA_real_
  
  tip_rows = which(
    p$data$isTip
  )
  
  p$data$ID[tip_rows] =
    tree$tip.label
  
  color_info =
    metadata %>%
    select(
      ID,
      color
    ) %>%
    filter(
      !is.na(color),
      color != ""
    ) %>%
    distinct(
      ID,
      .keep_all = TRUE
    )
  
  p$data =
    p$data %>%
    left_join(
      color_info,
      by = "ID"
    )
  
  internal_rows = which(
    !p$data$isTip
  )
  
  n_assign = min(
    length(internal_rows),
    length(bootstrap_values)
  )
  
  if (n_assign > 0) {
    
    p$data$bootstrap[
      internal_rows[
        seq_len(n_assign)
      ]
    ] =
      bootstrap_values[
        seq_len(n_assign)
      ]
  }
  
  tips_colored =
    p$data %>%
    filter(
      isTip,
      !is.na(color)
    )
  
  tips_uncolored =
    p$data %>%
    filter(
      isTip,
      is.na(color)
    )
  
  p =
    p +
    
    geom_tiplab(
      data = tips_uncolored,
      aes(
        label = ID
      ),
      color = "black",
      size = 4,
      hjust = 0
    ) +
    
    geom_tiplab(
      data = tips_colored,
      aes(
        label = ID,
        color = color
      ),
      size = 4,
      hjust = 0
    ) +
    
    scale_color_identity() +
    
    geom_treescale(
      fontsize = 5
    )
  
  high_bootstrap =
    p$data %>%
    filter(
      !isTip,
      !is.na(bootstrap),
      bootstrap > 90
    )
  
  if (nrow(high_bootstrap) > 0) {
    
    p =
      p +
      geom_nodepoint(
        data = high_bootstrap,
        aes(
          x = x,
          y = y
        ),
        colour = "black",
        size = 3.5,
        inherit.aes = FALSE
      )
  }
  
  p =
    p +
    theme(
      legend.position = "none",
      plot.margin = margin(
        10,
        20,
        10,
        10
      )
    )
  
  return(p)
}

for (tree_file in tree_files) {
  
  p =
    plot_gradient_tree(
      tree_file,
      info_upd
    )
  
  tree_name =
    tools::file_path_sans_ext(
      basename(tree_file)
    )
  
  p =
    p +
    ggtitle(tree_name)
  
  ggsave(
    paste0(
      tree_name,
      "_gradient.png"
    ),
    p,
    width = 20,
    height = 15,
    dpi = 300
  )
  
  ggsave(
    paste0(
      tree_name,
      "_gradient.pdf"
    ),
    p,
    width = 20,
    height = 15
  )
  
  ggsave(
    paste0(
      tree_name,
      "_gradient.svg"
    ),
    p,
    width = 20,
    height = 15
  )
}

cat(
  "\nВсе деревья обработаны.\n"
)