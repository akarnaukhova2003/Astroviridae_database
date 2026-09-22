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

# Приведение ID к единому виду для сопоставления дерева и metadata
normalize_id = function(x) {
  x = trimws(x)
  x = gsub("'", "", x)
  x = sub("/.*$", "", x)
  x = gsub("_", "-", x)
  x
}

# Получение нужного поля из исходного названия листа
# Если поля нет или оно пустое, возвращается NA
get_value = function(x, position) {
  parts = strsplit(x, "/", fixed = TRUE)[[1]]
  
  if (length(parts) < position) {
    return("NA")
  }
  
  value = trimws(parts[position])
  
  if (is.na(value) || value == "") {
    return("NA")
  }
  
  value
}

# Удаление пробелов из названий
clean_label = function(x) {
  if (is.na(x) || trimws(x) == "") {
    return("NA")
  }
  gsub("[[:space:]]+", "_", trimws(x))
}

# Загружаем metadata
metadata = read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata$ID_norm = normalize_id(metadata$ID)

# Загружаем референсное дерево и находим в нём outgroup
reference_tree = read.tree(reference_tree_file)

reference_outgroup = reference_tree$tip.label[
  normalize_id(reference_tree$tip.label) == outgroup_id
]

if (length(reference_outgroup) == 0) {
  stop(
    paste(
      "Outgroup",
      outgroup_id,
      "не найден в референсном дереве."
    )
  )
}

# Укореняем референсное дерево по outgroup
reference_tree_root = root(
  reference_tree,
  outgroup = reference_outgroup,
  resolve.root = TRUE
)

p_reference = ggtree(reference_tree_root)

# Получаем ID листьев, входящих в заданные клады
clade_taxa = lapply(
  clades,
  function(x) {
    normalize_id(
      get_taxa_name(
        p_reference,
        x
      )
    )
  }
)

# Создаём основной цвет для каждой клады
base_colors = distinctColorPalette(length(clades))
color_table = data.frame()

# Для каждой клады создаём градиент цвета от светлого к тёмному
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
      ID_norm = ids,
      color = colors
    )
  )
}

# Добавляем цвета к metadata
metadata = left_join(
  metadata,
  color_table,
  by = "ID_norm"
)

write.csv(
  metadata,
  "metadata_upd.csv",
  row.names = FALSE
)

# Обрабатываем каждое дерево
for (file in tree_files) {
  
  tree = read.tree(file)
  
  # Ищем outgroup в текущем дереве
  outgroup_label = tree$tip.label[
    normalize_id(tree$tip.label) == outgroup_id
  ]
  
  if (length(outgroup_label) == 0) {
    warning(
      paste(
        "Outgroup",
        outgroup_id,
        "не найден в",
        file
      )
    )
    next
  }
  
  # Укореняем дерево и удаляем outgroup
  tree_root = root(
    tree,
    outgroup = outgroup_label,
    resolve.root = TRUE
  )
  
  tree_final = drop.tip(
    tree_root,
    outgroup_label
  )
  
  p_final = ggtree(tree_final)
  
  # Сохраняем исходные названия листьев
  p_final$data$original_label = p_final$data$label
  
  # Создаём ID для сопоставления с metadata
  p_final$data$ID_norm = normalize_id(
    p_final$data$original_label
  )
  
  # Добавляем Host, class и цвет
  p_final$data = left_join(
    p_final$data,
    metadata %>%
      select(
        ID_norm,
        Host,
        class,
        color
      ),
    by = "ID_norm"
  )
  
  # Формируем новые названия только для листьев
  tip_rows = which(p_final$data$isTip)
  
  for (i in tip_rows) {
    
    original = p_final$data$original_label[i]
    
    # Берём 1-е, 3-е и 5-е значения из исходного названия
    value_1 = get_value(original, 1)
    value_3 = get_value(original, 3)
    value_5 = get_value(original, 5)
    
    # Добавляем Host и class из metadata
    host = p_final$data$Host[i]
    class = p_final$data$class[i]
    
    # Если данных нет, ставим NA и заменяем пробелы на _
    value_1 = clean_label(value_1)
    value_3 = clean_label(value_3)
    value_5 = clean_label(value_5)
    host = clean_label(host)
    class = clean_label(class)
    
    # Формат названия:
    # 1-е/3-е/5-е значение/Host/class
    p_final$data$label[i] = paste(
      value_1,
      value_3,
      value_5,
      host,
      class,
      sep = "/"
    )
  }
  
  # Листья без цвета окрашиваем в чёрный
  p_final$data$color[
    is.na(p_final$data$color)
  ] = "black"
  
  p_final = p_final +
    geom_tiplab(
      aes(
        label = label,
        color = color
      ),
      size = 4,
      hjust = 0
    ) +
    scale_color_identity()
  
  # Показываем bootstrap >95
  bootstrap_data = p_final$data[
    !p_final$data$isTip &
      !is.na(p_final$data$label) &
      suppressWarnings(
        as.numeric(p_final$data$label)
      ) > 95,
  ]
  
  p_final = p_final +
    geom_nodepoint(
      data = bootstrap_data,
      size = 3
    ) +
    geom_treescale()
  
  # Увеличиваем пространство справа для длинных названий
  max_x = max(
    p_final$data$x,
    na.rm = TRUE
  )
  
  p_final = p_final +
    xlim(
      0,
      max_x * 1.8
    ) +
    theme(
      plot.margin = margin(
        10,
        30,
        10,
        10
      )
    )
  
  # Имя выходного файла
  name = tools::file_path_sans_ext(
    basename(file)
  )
  
  # Сохраняем дерево в трёх форматах
  ggsave(
    paste0(name, "_gradient.png"),
    p_final,
    width = 15,
    height = 20,
    dpi = 300
  )
  
  ggsave(
    paste0(name, "_gradient.pdf"),
    p_final,
    width = 15,
    height = 20
  )
  
  ggsave(
    paste0(name, "_gradient.svg"),
    p_final,
    width = 15,
    height = 20
  )
}

cat("Все деревья обработаны.\n")
