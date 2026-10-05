library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(phytools)
library(randomcoloR)
library(colorspace)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)

get_arg = function(name) {
  
  i = match(
    name,
    args
  )
  
  if (
    is.na(i) ||
    i == length(args)
  ) {
    stop(
      paste(
        "Не найден аргумент:",
        name
      )
    )
  }
  
  args[i + 1]
}


# ============================================================
# ПАРАМЕТРЫ ИЗ SNAKEMAKE
# ============================================================

reference_tree_file = get_arg(
  "--reference-tree"
)

outgroup_id = get_arg(
  "--outgroup"
)

clades = as.numeric(
  strsplit(
    get_arg("--clades"),
    ",",
    fixed = TRUE
  )[[1]]
)

tree_files = strsplit(
  get_arg("--tree-files"),
  ",",
  fixed = TRUE
)[[1]]

metadata_file = get_arg(
  "--metadata"
)

cluster_file = get_arg(
  "--clusters"
)

host_colors_file = get_arg(
  "--host-colors"
)

output_dir = get_arg(
  "--output-dir"
)

done_file = get_arg(
  "--done-file"
)


source(
  host_colors_file
)


# ============================================================
# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# ============================================================

normalize_id = function(x) {
  
  x = trimws(x)
  x = gsub("'", "", x)
  
  sub(
    "/.*$",
    "",
    x
  )
}


get_value = function(x, position) {
  
  parts = strsplit(
    x,
    "/",
    fixed = TRUE
  )[[1]]
  
  if (length(parts) < position) {
    return("NA")
  }
  
  value = trimws(
    parts[position]
  )
  
  if (is.na(value) || value == "") {
    return("NA")
  }
  
  value
}


fix_na = function(x) {
  
  x = as.character(x)
  
  x[
    is.na(x) |
      trimws(x) == ""
  ] = "NA"
  
  x
}


# ============================================================
# ЗАГРУЗКА METADATA
# ============================================================

metadata = read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cluster_metadata = read.csv(
  cluster_file,
  sep = ";",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


metadata$ID_norm = normalize_id(
  metadata$ID
)

cluster_metadata$id = trimws(
  as.character(
    cluster_metadata$id
  )
)


cluster_metadata$aa83_1B = fix_na(
  cluster_metadata$aa83_1B
)

cluster_metadata$aa90_1B = fix_na(
  cluster_metadata$aa90_1B
)

cluster_metadata$nt86_1B = fix_na(
  cluster_metadata$nt86_1B
)

cluster_metadata$nt80_1B = fix_na(
  cluster_metadata$nt80_1B
)

cluster_metadata$Host_shorter = fix_na(
  cluster_metadata$Host_shorter
)


# ============================================================
# РЕФЕРЕНСНОЕ ДЕРЕВО
# ============================================================

reference_tree = read.tree(
  reference_tree_file
)

reference_tree$edge.length[
  is.na(
    reference_tree$edge.length
  )
] = 0


reference_outgroup = reference_tree$tip.label[
  normalize_id(
    reference_tree$tip.label
  ) ==
    normalize_id(
      outgroup_id
    )
]


if (length(reference_outgroup) != 1) {
  
  stop(
    paste0(
      "Аутгруппа не найдена или найдена несколько раз: ",
      outgroup_id
    )
  )
}


reference_tree = root(
  reference_tree,
  outgroup = reference_outgroup,
  resolve.root = TRUE
)


# ============================================================
# ЦВЕТА КЛАД
# ============================================================

clade_colors = distinctColorPalette(
  length(clades)
)


# ============================================================
# ЛИСТЬЯ КАЖДОЙ КЛАДЫ
# ============================================================

reference_clade_ids = vector(
  "list",
  length(clades)
)


for (i in seq_along(clades)) {
  
  clade_tree = extract.clade(
    reference_tree,
    clades[i]
  )
  
  reference_clade_ids[[i]] =
    normalize_id(
      clade_tree$tip.label
    )
  
  cat(
    "Клада",
    clades[i],
    ":",
    length(
      reference_clade_ids[[i]]
    ),
    "листьев\n"
  )
}


# ============================================================
# СОЗДАНИЕ GGTREE РЕФЕРЕНСНОГО ДЕРЕВА
# ============================================================

reference_tree_plot = ggtree(
  reference_tree,
  size = 0.4
)


reference_tip_rows = which(
  reference_tree_plot$data$isTip
)


reference_tree_plot$data$ID_norm =
  NA_character_


reference_tree_plot$data$ID_norm[
  reference_tip_rows
] = normalize_id(
  reference_tree_plot$data$label[
    reference_tip_rows
  ]
)


# ============================================================
# СОЗДАНИЕ ФИКСИРОВАННЫХ ГРАДИЕНТОВ
# ============================================================

reference_tip_colors = data.frame(
  ID_norm = character(),
  tip_color = character(),
  stringsAsFactors = FALSE
)


for (i in seq_along(clades)) {
  
  clade_rows = reference_tip_rows[
    reference_tree_plot$data$ID_norm[
      reference_tip_rows
    ] %in%
      reference_clade_ids[[i]]
  ]
  
  
  if (!length(clade_rows)) {
    next
  }
  
  
  # Порядок листьев определяется ТОЛЬКО
  # по референсному дереву
  
  clade_rows = clade_rows[
    order(
      reference_tree_plot$data$y[
        clade_rows
      ]
    )
  ]
  
  
  # Градиент создаётся ОДИН РАЗ
  
  gradient = colorRampPalette(
    c(
      lighten(
        clade_colors[i],
        amount = 0.35
      ),
      darken(
        clade_colors[i],
        amount = 0.55
      )
    )
  )(
    length(clade_rows)
  )
  
  # Сохраняем соответствие:
  # конкретный ID -> конкретный цвет
  
  reference_tip_colors = rbind(
    reference_tip_colors,
    data.frame(
      ID_norm =
        reference_tree_plot$data$ID_norm[
          clade_rows
        ],
      tip_color = gradient,
      stringsAsFactors = FALSE
    )
  )
}


# Убираем возможные дубли ID

reference_tip_colors =
  reference_tip_colors[
    !duplicated(
      reference_tip_colors$ID_norm
    ),
  ]


# ============================================================
# ПРОВЕРКА ЦВЕТОВ РЕФЕРЕНСНОГО ДЕРЕВА
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "Создана фиксированная карта цветов\n"
)

cat(
  "Количество окрашенных листьев:",
  nrow(
    reference_tip_colors
  ),
  "\n"
)

cat(
  "========================================\n\n"
)


print(
  head(
    reference_tip_colors,
    20
  )
)


# ============================================================
# УРОВНИ HEATMAP
# ============================================================

levels_f1 = unique(
  cluster_metadata$aa83_1B
)

levels_f2 = unique(
  cluster_metadata$aa90_1B
)

levels_f3 = unique(
  cluster_metadata$nt86_1B
)

levels_f4 = unique(
  cluster_metadata$nt80_1B
)


# ============================================================
# ЦВЕТА HEATMAP
# ============================================================

color_map_f1 = setNames(
  distinctColorPalette(
    length(levels_f1)
  ),
  levels_f1
)

color_map_f2 = setNames(
  distinctColorPalette(
    length(levels_f2)
  ),
  levels_f2
)

color_map_f3 = setNames(
  distinctColorPalette(
    length(levels_f3)
  ),
  levels_f3
)

color_map_f4 = setNames(
  distinctColorPalette(
    length(levels_f4)
  ),
  levels_f4
)


# ============================================================
# HOST
# ============================================================

color_map_host = manual_host_colors

levels_host = host_order


# ============================================================
# ПАРАМЕТРЫ ВИЗУАЛИЗАЦИИ
# ============================================================

tip_label_size = 5

tree_line_size = 0.4

heatmap_offset = 0.8

heatmap_gap = 0.13

heatmap_width = 0.045


# ============================================================
# ОБРАБОТКА ДЕРЕВЬЕВ
# ============================================================

for (file in tree_files) {
  
  cat(
    "\n========================================\n"
  )
  
  cat(
    "Обрабатывается:",
    file,
    "\n"
  )
  
  cat(
    "========================================\n"
  )
  
  
  # ==========================================================
  # ЧТЕНИЕ ДЕРЕВА
  # ==========================================================
  
  tree = read.tree(
    file
  )
  
  
  tree$edge.length[
    is.na(
      tree$edge.length
    )
  ] = 0
  
  
  original_labels = tree$tip.label
  
  
  tree$tip.label = paste0(
    original_labels,
    "__tip",
    seq_along(
      original_labels
    )
  )
  
  
  # ==========================================================
  # АУТГРУППА
  # ==========================================================
  
  outgroup_position = which(
    normalize_id(
      original_labels
    ) ==
      normalize_id(
        outgroup_id
      )
  )
  
  
  if (length(outgroup_position) != 1) {
    
    cat(
      "Аутгруппа не найдена:",
      file,
      "\n"
    )
    
    next
  }
  
  
  outgroup_label = tree$tip.label[
    outgroup_position
  ]
  
  
  tree = root(
    tree,
    outgroup = outgroup_label,
    resolve.root = TRUE
  )
  
  
  tree = drop.tip(
    tree,
    outgroup_label
  )
  
  
  # ==========================================================
  # СОЗДАНИЕ GGTREE
  # ==========================================================
  
  p = ggtree(
    tree,
    size = tree_line_size
  )
  
  
  tip_rows = which(
    p$data$isTip
  )
  
  
  # ==========================================================
  # ИСХОДНЫЕ НАЗВАНИЯ ЛИСТЬЕВ
  # ==========================================================
  
  p$data$original_label =
    NA_character_
  
  
  p$data$original_label[
    tip_rows
  ] = sub(
    "__tip[0-9]+$",
    "",
    p$data$label[
      tip_rows
    ]
  )
  
  
  # ==========================================================
  # НОРМАЛИЗОВАННЫЙ ID
  # ==========================================================
  
  p$data$ID_norm =
    NA_character_
  
  
  p$data$ID_norm[
    tip_rows
  ] = normalize_id(
    p$data$original_label[
      tip_rows
    ]
  )
  
  
  tip_data = p$data[
    p$data$isTip,
  ]
  
  
  tree_tip_ids = tip_data$ID_norm
  
  
  # ==========================================================
  # MAPPING С CLUSTER CSV
  # ==========================================================
  
  cluster_data = cluster_metadata[
    match(
      tree_tip_ids,
      cluster_metadata$id
    ),
  ]
  
  
  cat(
    "\n=== ПРОВЕРКА MAPPING ===\n"
  )
  
  
  cat(
    "Всего листьев:",
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
  
  
  mapping_check = data.frame(
    tree_id = tree_tip_ids,
    cluster_id = cluster_data$id,
    aa83 = cluster_data$aa83_1B,
    aa90 = cluster_data$aa90_1B,
    nt86 = cluster_data$nt86_1B,
    nt80 = cluster_data$nt80_1B,
    host = cluster_data$Host_shorter
  )
  
  
  print(
    head(
      mapping_check,
      20
    )
  )
  
  
  # ==========================================================
  # ОСНОВНОЙ METADATA
  # ==========================================================
  
  p$data$Host_species =
    NA_character_
  
  p$data$class =
    NA_character_
  
  
  metadata_match = match(
    p$data$ID_norm[
      tip_rows
    ],
    metadata$ID_norm
  )
  
  
  p$data$Host_species[
    tip_rows
  ] = metadata$Host_species[
    metadata_match
  ]
  
  
  p$data$class[
    tip_rows
  ] = metadata$class[
    metadata_match
  ]
  
  
  p$data$Host_species =
    fix_na(
      p$data$Host_species
    )
  
  
  p$data$class =
    fix_na(
      p$data$class
    )
  
  
  # ==========================================================
  # ЦВЕТА ЛИСТЬЕВ
  # ==========================================================
  
  p$data$tip_color =
    "black"
  
  
  color_match = match(
    p$data$ID_norm[
      tip_rows
    ],
    reference_tip_colors$ID_norm
  )
  
  
  p$data$tip_color[
    tip_rows
  ] = reference_tip_colors$tip_color[
    color_match
  ]
  
  
  p$data$tip_color[
    tip_rows[
      is.na(
        color_match
      )
    ]
  ] = "black"
  
  
  # ==========================================================
  # НАЗВАНИЯ ЛИСТЬЕВ
  # ==========================================================
  
  p$data$display_label =
    NA_character_
  
  
  p$data$display_label[
    tip_rows
  ] = paste0(
    
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
  
  
  # ==========================================================
  # ПОДПИСИ ЛИСТЬЕВ
  # ==========================================================
  
  p = p +
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
  
  
  # ==========================================================
  # BOOTSTRAP > 95
  # ==========================================================
  
  bootstrap_values =
    suppressWarnings(
      as.numeric(
        p$data$label
      )
    )
  
  
  bootstrap_data = p$data[
    !p$data$isTip &
      !is.na(
        bootstrap_values
      ) &
      bootstrap_values > 95,
  ]
  
  
  if (
    nrow(
      bootstrap_data
    ) > 0
  ) {
    
    p = p +
      geom_point2(
        data = bootstrap_data,
        aes(
          subset = TRUE
        ),
        size = 1.5
      )
  }
  
  
  # ==========================================================
  # ДАННЫЕ ДЛЯ HEATMAP
  # ==========================================================
  
  cluster1 = data.frame(
    AA83 = factor(
      cluster_data$aa83_1B,
      levels = levels_f1
    )
  )
  
  
  cluster2 = data.frame(
    AA90 = factor(
      cluster_data$aa90_1B,
      levels = levels_f2
    )
  )
  
  
  cluster3 = data.frame(
    NT86 = factor(
      cluster_data$nt86_1B,
      levels = levels_f3
    )
  )
  
  
  cluster4 = data.frame(
    NT80 = factor(
      cluster_data$nt80_1B,
      levels = levels_f4
    )
  )
  
  
  # ==========================================================
  # HOST HEATMAP
  # ==========================================================
  
  host_values = cluster_data$Host_shorter
  
  
  host_values[
    !host_values %in%
      names(manual_host_colors)
  ] = "Other"
  
  
  host_data = data.frame(
    Host = factor(
      host_values,
      levels = levels_host
    )
  )
  
  
  # ==========================================================
  # ROW NAMES HEATMAP
  # ==========================================================
  
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
  
  
  # ==========================================================
  # HEATMAP AA83
  # ==========================================================
  
  p = gheatmap(
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
  
  
  # ==========================================================
  # HEATMAP AA90
  # ==========================================================
  
  p = gheatmap(
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
  
  
  # ==========================================================
  # HEATMAP NT86
  # ==========================================================
  
  p = gheatmap(
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
  
  
  # ==========================================================
  # HEATMAP NT80
  # ==========================================================
  
  p = gheatmap(
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
  
  
  # ==========================================================
  # HEATMAP HOST
  # ==========================================================
  
  p = gheatmap(
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
      values = manual_host_colors,
      breaks = host_order,
      name = "Host",
      drop = FALSE,
      na.value = "grey85"
    ) +
    theme(
      axis.text.x = element_blank()
    )
  
  
  # ==========================================================
  # ОФОРМЛЕНИЕ
  # ==========================================================
  
  p = p +
    theme(
      plot.margin = margin(
        10,
        100,
        10,
        10
      )
    )
  
  
  # ==========================================================
  # СОХРАНЕНИЕ
  # ==========================================================
  
  output_name = sub(
    "\\.treefile$",
    "",
    basename(file)
  )
  
  
  # PDF
  
  ggsave(
    file.path(
      output_dir,
      paste0(
        output_name,
        "_gradient_heatmap.pdf"
      )
    ),
    p,
    width = 16,
    height = 14
  )
  
  
  # SVG
  
  ggsave(
    file.path(
      output_dir,
      paste0(
        output_name,
        "_gradient_heatmap.svg"
      )
    ),
    p,
    width = 16,
    height = 14
  )
  
  
  cat(
    "\nГотово:",
    file,
    "\n"
  )
}


# ============================================================
# MARKER ДЛЯ SNAKEMAKE
# ============================================================

file.create(
  done_file
)


cat(
  "\nВсе деревья обработаны.\n"
)
