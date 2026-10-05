library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(randomcoloR)
library(dplyr)


# ============================================================
# ПАРАМЕТРЫ
# ============================================================

tree_file =
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/Aves_Amphibia_trees/Aves_1B_97_aa_2206_after_trim_iqtree_clr.nex"

metadata_file =
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/snakemake/constants/Astroviridae_Aves_Amphibia_14072026_host.csv"

cluster_file =
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/snakemake/constants/Astroviridae_Aves_30042026_clusters_shorter.csv"

hosts_colors_file =
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/snakemake/constants/host_colors.R"

outgroup_id =
  "OM480533"

output_dir =
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/snakemake/results/visualization/Aves_1B_97"


# ============================================================
# ЦВЕТА CLADE
# ============================================================

clade_color_map = c(
  G = "#2E8B57",
  B = "#2878B5",
  P = "#8E44AD",
  Y = "#E5B700",
  R = "#C0392B"
)


# ============================================================
# ПАРАМЕТРЫ ВИЗУАЛИЗАЦИИ
# ============================================================

tip_label_size =
  3.5

tree_line_size =
  0.4

heatmap_offset =
  0.8

heatmap_gap =
  0.13

heatmap_width =
  0.045


# ============================================================
# ПРОВЕРКА ФАЙЛОВ
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "ПРОВЕРКА ФАЙЛОВ\n"
)

cat(
  "========================================\n"
)

cat(
  "Дерево:",
  tree_file,
  "\n"
)

cat(
  "Metadata:",
  metadata_file,
  "\n"
)

cat(
  "Clusters:",
  cluster_file,
  "\n"
)

cat(
  "Host colors:",
  hosts_colors_file,
  "\n"
)


if (!file.exists(tree_file)) {
  
  stop(
    paste0(
      "\nНе найден файл дерева:\n",
      tree_file,
      "\n"
    )
  )
}


if (!file.exists(metadata_file)) {
  
  stop(
    paste0(
      "\nНе найден metadata файл:\n",
      metadata_file,
      "\n"
    )
  )
}


if (!file.exists(cluster_file)) {
  
  stop(
    paste0(
      "\nНе найден cluster файл:\n",
      cluster_file,
      "\n"
    )
  )
}


if (!file.exists(hosts_colors_file)) {
  
  stop(
    paste0(
      "\nНе найден файл hosts_colors.R:\n",
      hosts_colors_file,
      "\n"
    )
  )
}


# ============================================================
# ЗАГРУЗКА HOST COLORS
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "ЗАГРУЗКА HOST COLORS\n"
)

cat(
  "========================================\n"
)


source(
  hosts_colors_file
)


host_colors =
  manual_host_colors[
    host_order
  ]


cat(
  "Загружено цветов Host:",
  length(
    host_colors
  ),
  "\n"
)


cat(
  "Категорий Host:",
  length(
    host_order
  ),
  "\n"
)


# ============================================================
# ФУНКЦИИ
# ============================================================

normalize_id = function(x) {
  
  x =
    trimws(
      x
    )
  
  x =
    gsub(
      "'",
      "",
      x
    )
  
  x =
    gsub(
      '"',
      "",
      x
    )
  
  sub(
    "/.*$",
    "",
    x
  )
}


fix_na = function(x) {
  
  x =
    as.character(
      x
    )
  
  x[
    is.na(x) |
      trimws(x) == ""
  ] =
    "NA"
  
  x
}


get_value = function(x, position) {
  
  parts =
    strsplit(
      x,
      "/",
      fixed = TRUE
    )[[1]]
  
  if (
    length(parts) < position
  ) {
    
    return(
      "NA"
    )
  }
  
  value =
    trimws(
      parts[position]
    )
  
  if (
    is.na(value) ||
    value == ""
  ) {
    
    return(
      "NA"
    )
  }
  
  value
}


# ============================================================
# METADATA
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "ЗАГРУЗКА METADATA\n"
)

cat(
  "========================================\n"
)


metadata =
  read.csv(
    metadata_file,
    sep = ",",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


required_metadata_columns =
  c(
    "ID",
    "Host",
    "class",
    "Host_species",
    "clade_color"
  )


missing_metadata_columns =
  setdiff(
    required_metadata_columns,
    names(metadata)
  )


if (
  length(
    missing_metadata_columns
  ) > 0
) {
  
  stop(
    paste0(
      "\nВ metadata отсутствуют колонки: ",
      paste(
        missing_metadata_columns,
        collapse = ", "
      ),
      "\n"
    )
  )
}


metadata$ID_norm =
  normalize_id(
    metadata$ID
  )


metadata$Host =
  fix_na(
    metadata$Host
  )


metadata$class =
  fix_na(
    metadata$class
  )


metadata$Host_species =
  fix_na(
    metadata$Host_species
  )


metadata$clade_color =
  fix_na(
    metadata$clade_color
  )


cat(
  "Строк metadata:",
  nrow(metadata),
  "\n"
)


# ============================================================
# CLUSTER METADATA
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "ЗАГРУЗКА CLUSTER METADATA\n"
)

cat(
  "========================================\n"
)


cluster_metadata =
  read.csv(
    cluster_file,
    sep = ";",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


required_cluster_columns =
  c(
    "id",
    "aa83_1B",
    "aa90_1B",
    "nt86_1B",
    "nt80_1B",
    "Host_shorter"
  )


missing_cluster_columns =
  setdiff(
    required_cluster_columns,
    names(cluster_metadata)
  )


if (
  length(
    missing_cluster_columns
  ) > 0
) {
  
  stop(
    paste0(
      "\nВ cluster metadata отсутствуют колонки: ",
      paste(
        missing_cluster_columns,
        collapse = ", "
      ),
      "\n"
    )
  )
}


cluster_metadata$id =
  trimws(
    as.character(
      cluster_metadata$id
    )
  )


cluster_metadata$aa83_1B =
  fix_na(
    cluster_metadata$aa83_1B
  )


cluster_metadata$aa90_1B =
  fix_na(
    cluster_metadata$aa90_1B
  )


cluster_metadata$nt86_1B =
  fix_na(
    cluster_metadata$nt86_1B
  )


cluster_metadata$nt80_1B =
  fix_na(
    cluster_metadata$nt80_1B
  )


cluster_metadata$Host_shorter =
  fix_na(
    cluster_metadata$Host_shorter
  )


# ============================================================
# ПРОВЕРКА HOST
# ============================================================

unknown_hosts =
  setdiff(
    unique(
      cluster_metadata$Host_shorter
    ),
    host_order
  )


if (
  length(
    unknown_hosts
  ) > 0
) {
  
  cat(
    "\n========================================\n"
  )
  
  cat(
    "НЕИЗВЕСТНЫЕ HOST\n"
  )
  
  cat(
    "========================================\n"
  )
  
  cat(
    paste(
      unknown_hosts,
      collapse = "\n"
    ),
    "\n"
  )
  
  cat(
    "\nОни будут заменены на Other.\n"
  )
  
  cluster_metadata$Host_shorter[
    cluster_metadata$Host_shorter %in% unknown_hosts
  ] =
    "Other"
}


# ============================================================
# ЧТЕНИЕ NEXUS
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "ЧТЕНИЕ NEXUS\n"
)

cat(
  "========================================\n"
)


nexus_text =
  paste(
    readLines(
      tree_file,
      warn = FALSE
    ),
    collapse = "\n"
  )


if (
  !grepl(
    "#NEXUS",
    nexus_text,
    ignore.case = TRUE
  )
) {
  
  stop(
    "\nФайл не содержит #NEXUS.\n"
  )
}


# ============================================================
# ИЗВЛЕКАЕМ TREE
# ============================================================

tree_match =
  regexec(
    "(?is)tree\\s+[^=]+\\s*=\\s*(?:\\[&R\\]\\s*)?(.*?);",
    nexus_text,
    perl = TRUE
  )


tree_match_result =
  regmatches(
    nexus_text,
    tree_match
  )


if (
  length(tree_match_result) == 0 ||
  length(tree_match_result[[1]]) < 2
) {
  
  stop(
    "\nНе удалось найти tree внутри NEXUS.\n"
  )
}


tree_newick =
  tree_match_result[[1]][2]


cat(
  "Tree из NEXUS найден.\n"
)


# ============================================================
# УДАЛЯЕМ TIP COLOR ANNOTATIONS
# ============================================================

tree_newick =
  gsub(
    "\\[&!color=#[0-9A-Fa-f]+\\]",
    "",
    tree_newick,
    perl = TRUE
  )


# ============================================================
# ПРЕОБРАЗУЕМ BOOTSTRAP
# ============================================================

tree_newick =
  gsub(
    "\\[&label=([0-9.]+)\\]",
    "\\1",
    tree_newick,
    perl = TRUE
  )


# ============================================================
# ЧИТАЕМ NEWICK
# ============================================================

tree =
  read.tree(
    text =
      paste0(
        tree_newick,
        ";"
      )
  )


tree$edge.length[
  is.na(
    tree$edge.length
  )
] =
  0


cat(
  "Количество tips:",
  length(
    tree$tip.label
  ),
  "\n"
)

cat(
  "Количество внутренних узлов:",
  tree$Nnode,
  "\n"
)


# ============================================================
# BOOTSTRAP
# ============================================================

bootstrap_values =
  suppressWarnings(
    as.numeric(
      tree$node.label
    )
  )


cat(
  "Внутренних узлов с bootstrap:",
  sum(
    !is.na(
      bootstrap_values
    )
  ),
  "\n"
)


cat(
  "Внутренних узлов без bootstrap:",
  sum(
    is.na(
      bootstrap_values
    )
  ),
  "\n"
)


cat(
  "Bootstrap > 95:",
  sum(
    bootstrap_values > 95,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# ORIGINAL LABEL
# ============================================================

original_labels =
  tree$tip.label


# ============================================================
# ПОИСК OUTGROUP
# ============================================================

cat(
  "\nПоиск outgroup:",
  outgroup_id,
  "\n"
)


outgroup_position =
  which(
    normalize_id(
      original_labels
    ) ==
      normalize_id(
        outgroup_id
      )
  )


if (
  length(
    outgroup_position
  ) != 1
) {
  
  stop(
    paste0(
      "\nOutgroup ",
      outgroup_id,
      " не найден или найден ",
      length(
        outgroup_position
      ),
      " раз.\n"
    )
  )
}


cat(
  "Outgroup найден.\n"
)


# ============================================================
# ДЕЛАЕМ TIP LABEL УНИКАЛЬНЫМИ
# ============================================================

tree$tip.label =
  paste0(
    original_labels,
    "__tip",
    seq_along(
      original_labels
    )
  )


outgroup_label =
  tree$tip.label[
    outgroup_position
  ]


# ============================================================
# ROOT ПО OUTGROUP
# ============================================================

cat(
  "Root по",
  outgroup_id,
  "...\n"
)


tree =
  root(
    tree,
    outgroup = outgroup_label,
    resolve.root = TRUE
  )


# ============================================================
# УДАЛЯЕМ OUTGROUP
# ============================================================

tree =
  drop.tip(
    tree,
    outgroup_label
  )


cat(
  "Outgroup удалён.\n"
)


# ============================================================
# GGTREE
# ============================================================

p =
  ggtree(
    tree,
    size = tree_line_size
  )


tip_rows =
  which(
    p$data$isTip
  )


internal_rows =
  which(
    !p$data$isTip
  )


# ============================================================
# BOOTSTRAP В GGTREE
# ============================================================

p$data$bootstrap =
  NA_real_


p$data$bootstrap[
  internal_rows
] =
  suppressWarnings(
    as.numeric(
      p$data$label[
        internal_rows
      ]
    )
  )


cat(
  "\n========================================\n"
)

cat(
  "BOOTSTRAP ПОСЛЕ ROOT\n"
)

cat(
  "========================================\n"
)


cat(
  "Внутренних узлов:",
  length(
    internal_rows
  ),
  "\n"
)


cat(
  "Узлов с bootstrap:",
  sum(
    !is.na(
      p$data$bootstrap[
        internal_rows
      ]
    )
  ),
  "\n"
)


cat(
  "Bootstrap > 95:",
  sum(
    p$data$bootstrap[
      internal_rows
    ] > 95,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# BOOTSTRAP POINTS
# ============================================================

bootstrap_data =
  p$data[
    !p$data$isTip &
      !is.na(
        p$data$bootstrap
      ) &
      p$data$bootstrap > 95,
  ]


if (
  nrow(
    bootstrap_data
  ) > 0
) {
  
  p =
    p +
    geom_point2(
      data = bootstrap_data,
      aes(
        subset = TRUE
      ),
      size = 2
    )
}


# ============================================================
# ORIGINAL LABEL
# ============================================================

p$data$original_label =
  NA_character_


p$data$original_label[
  tip_rows
] =
  sub(
    "__tip[0-9]+$",
    "",
    p$data$label[
      tip_rows
    ]
  )


p$data$original_label[
  tip_rows
] =
  gsub(
    '["\']',
    "",
    p$data$original_label[
      tip_rows
    ]
  )


# ============================================================
# NORMALIZED ID
# ============================================================

p$data$ID_norm =
  NA_character_


p$data$ID_norm[
  tip_rows
] =
  normalize_id(
    p$data$original_label[
      tip_rows
    ]
  )


# ============================================================
# METADATA MAPPING
# ============================================================

metadata_match =
  match(
    p$data$ID_norm[
      tip_rows
    ],
    metadata$ID_norm
  )


p$data$Host_species =
  NA_character_


p$data$Host =
  NA_character_


p$data$class =
  NA_character_


p$data$clade_color =
  NA_character_


p$data$Host_species[
  tip_rows
] =
  metadata$Host_species[
    metadata_match
  ]


p$data$Host[
  tip_rows
] =
  metadata$Host[
    metadata_match
  ]


p$data$class[
  tip_rows
] =
  metadata$class[
    metadata_match
  ]


p$data$clade_color[
  tip_rows
] =
  metadata$clade_color[
    metadata_match
  ]


# ============================================================
# HOST SPECIES FALLBACK
# ============================================================

p$data$Host_species[
  tip_rows
] =
  ifelse(
    is.na(
      p$data$Host_species[
        tip_rows
      ]
    ) |
      p$data$Host_species[
        tip_rows
      ] == "NA",
    
    p$data$Host[
      tip_rows
    ],
    
    p$data$Host_species[
      tip_rows
    ]
  )


# ============================================================
# ЦВЕТ TIP LABEL
# ============================================================

p$data$tip_color =
  clade_color_map[
    p$data$clade_color
  ]


p$data$tip_color[
  is.na(
    p$data$tip_color
  )
] =
  "black"


# ============================================================
# ПРОВЕРКА METADATA
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "METADATA MAPPING\n"
)

cat(
  "========================================\n"
)


cat(
  "Всего tips:",
  length(
    tip_rows
  ),
  "\n"
)


cat(
  "Сопоставлено:",
  sum(
    !is.na(
      metadata_match
    )
  ),
  "\n"
)


cat(
  "Не сопоставлено:",
  sum(
    is.na(
      metadata_match
    )
  ),
  "\n"
)


if (
  any(
    is.na(
      metadata_match
    )
  )
) {
  
  cat(
    "\nНе найдены в metadata:\n"
  )
  
  print(
    p$data$ID_norm[
      tip_rows[
        is.na(
          metadata_match
        )
      ]
    ]
  )
}


# ============================================================
# DISPLAY LABEL
# ============================================================

p$data$display_label =
  NA_character_


p$data$display_label[
  tip_rows
] =
  paste0(
    
    sapply(
      p$data$original_label[
        tip_rows
      ],
      get_value,
      position = 1
    ),
    
    "/",
    
    sapply(
      p$data$original_label[
        tip_rows
      ],
      get_value,
      position = 5
    ),
    
    "/",
    
    p$data$Host_species[
      tip_rows
    ],
    
    "/",
    
    p$data$class[
      tip_rows
    ]
  )


p$data$display_label[
  tip_rows
] =
  gsub(
    '["\']',
    "",
    p$data$display_label[
      tip_rows
    ]
  )


# ============================================================
# TIP LABELS
# ============================================================

p =
  p +
  geom_tiplab(
    aes(
      label = display_label,
      color = tip_color
    ),
    size = tip_label_size,
    offset = 0.02,
    align = FALSE
  ) +
  scale_color_identity()


# ============================================================
# CLUSTER MAPPING
# ============================================================

tree_tip_ids =
  p$data$ID_norm[
    tip_rows
  ]


cluster_data =
  cluster_metadata[
    match(
      tree_tip_ids,
      cluster_metadata$id
    ),
  ]


cat(
  "\n========================================\n"
)

cat(
  "CLUSTER MAPPING\n"
)

cat(
  "========================================\n"
)


cat(
  "Всего tips:",
  length(
    tree_tip_ids
  ),
  "\n"
)


cat(
  "Сопоставлено:",
  sum(
    !is.na(
      cluster_data$id
    )
  ),
  "\n"
)


cat(
  "Не сопоставлено:",
  sum(
    is.na(
      cluster_data$id
    )
  ),
  "\n"
)


if (
  any(
    is.na(
      cluster_data$id
    )
  )
) {
  
  cat(
    "\nНе найдены в cluster CSV:\n"
  )
  
  print(
    tree_tip_ids[
      is.na(
        cluster_data$id
      )
    ]
  )
}


# ============================================================
# LEVELS CLUSTERS
# ============================================================

levels_f1 =
  unique(
    cluster_metadata$aa83_1B
  )


levels_f2 =
  unique(
    cluster_metadata$aa90_1B
  )


levels_f3 =
  unique(
    cluster_metadata$nt86_1B
  )


levels_f4 =
  unique(
    cluster_metadata$nt80_1B
  )


# ============================================================
# COLORS CLUSTERS
# ============================================================

color_map_f1 =
  setNames(
    distinctColorPalette(
      length(
        levels_f1
      )
    ),
    levels_f1
  )


color_map_f2 =
  setNames(
    distinctColorPalette(
      length(
        levels_f2
      )
    ),
    levels_f2
  )


color_map_f3 =
  setNames(
    distinctColorPalette(
      length(
        levels_f3
      )
    ),
    levels_f3
  )


color_map_f4 =
  setNames(
    distinctColorPalette(
      length(
        levels_f4
      )
    ),
    levels_f4
  )


# ============================================================
# FIXED HOST COLORS
# ============================================================

host_colors =
  manual_host_colors[
    host_order
  ]


# ============================================================
# HEATMAP DATA
# ============================================================

cluster1 =
  data.frame(
    AA83 = factor(
      cluster_data$aa83_1B,
      levels = levels_f1
    )
  )


cluster2 =
  data.frame(
    AA90 = factor(
      cluster_data$aa90_1B,
      levels = levels_f2
    )
  )


cluster3 =
  data.frame(
    NT86 = factor(
      cluster_data$nt86_1B,
      levels = levels_f3
    )
  )


cluster4 =
  data.frame(
    NT80 = factor(
      cluster_data$nt80_1B,
      levels = levels_f4
    )
  )


host_data =
  data.frame(
    Host = factor(
      cluster_data$Host_shorter,
      levels = host_order
    )
  )


# ============================================================
# ROW NAMES
# ============================================================

heatmap_rows =
  tree$tip.label


rownames(cluster1) =
  heatmap_rows


rownames(cluster2) =
  heatmap_rows


rownames(cluster3) =
  heatmap_rows


rownames(cluster4) =
  heatmap_rows


rownames(host_data) =
  heatmap_rows


# ============================================================
# HEATMAP AA83
# ============================================================

p =
  gheatmap(
    p,
    cluster1,
    width = heatmap_width,
    offset = heatmap_offset,
    colnames = TRUE,
    font.size = 2,
    colnames_position = "top"
  ) +
  scale_fill_manual(
    values = color_map_f1,
    guide = "none",
    na.value = "grey90"
  ) +
  theme(
    axis.text.x = element_blank()
  ) +
  new_scale_fill()


# ============================================================
# HEATMAP AA90
# ============================================================

p =
  gheatmap(
    p,
    cluster2,
    width = heatmap_width,
    offset =
      heatmap_offset +
      heatmap_gap,
    colnames = TRUE,
    font.size = 2,
    colnames_position = "top"
  ) +
  scale_fill_manual(
    values = color_map_f2,
    guide = "none",
    na.value = "grey90"
  ) +
  theme(
    axis.text.x = element_blank()
  ) +
  new_scale_fill()


# ============================================================
# HEATMAP NT86
# ============================================================

p =
  gheatmap(
    p,
    cluster3,
    width = heatmap_width,
    offset =
      heatmap_offset +
      heatmap_gap * 2,
    colnames = TRUE,
    font.size = 2,
    colnames_position = "top"
  ) +
  scale_fill_manual(
    values = color_map_f3,
    guide = "none",
    na.value = "grey90"
  ) +
  theme(
    axis.text.x = element_blank()
  ) +
  new_scale_fill()


# ============================================================
# HEATMAP NT80
# ============================================================

p =
  gheatmap(
    p,
    cluster4,
    width = heatmap_width,
    offset =
      heatmap_offset +
      heatmap_gap * 3,
    colnames = TRUE,
    font.size = 2,
    colnames_position = "top"
  ) +
  scale_fill_manual(
    values = color_map_f4,
    guide = "none",
    na.value = "grey90"
  ) +
  theme(
    axis.text.x = element_blank()
  ) +
  new_scale_fill()


# ============================================================
# HEATMAP HOST
# ============================================================

p =
  gheatmap(
    p,
    host_data,
    width = heatmap_width,
    offset =
      heatmap_offset +
      heatmap_gap * 4,
    colnames = TRUE,
    font.size = 2,
    colnames_position = "top"
  ) +
  scale_fill_manual(
    values = host_colors,
    breaks = host_order,
    limits = host_order,
    name = "Host",
    drop = FALSE,
    na.value = "grey85"
  ) +
  theme(
    axis.text.x = element_blank()
  )


# ============================================================
# ОФОРМЛЕНИЕ
# ============================================================

p =
  p +
  theme(
    plot.margin =
      margin(
        10,
        100,
        10,
        10
      )
  )


# ============================================================
# OUTPUT
# ============================================================

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


output_name =
  basename(
    tree_file
  )


output_name =
  sub(
    "\\.(nex|nexus|nwk)$",
    "",
    output_name,
    ignore.case = TRUE
  )


# ============================================================
# PDF
# ============================================================

pdf_file =
  file.path(
    output_dir,
    paste0(
      output_name,
      "_clade_heatmap.pdf"
    )
  )


ggsave(
  filename = pdf_file,
  plot = p,
  width = 16,
  height = 14
)


# ============================================================
# SVG
# ============================================================

svg_file =
  file.path(
    output_dir,
    paste0(
      output_name,
      "_clade_heatmap.svg"
    )
  )


ggsave(
  filename = svg_file,
  plot = p,
  width = 16,
  height = 14
)


# ============================================================
# ГОТОВО
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "ГОТОВО\n"
)

cat(
  "========================================\n"
)

cat(
  "PDF:\n",
  pdf_file,
  "\n\n"
)

cat(
  "SVG:\n",
  svg_file,
  "\n"
)

cat(
  "========================================\n"
)