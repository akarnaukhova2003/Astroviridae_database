library(randomcoloR)
library(colorspace)
library(dplyr)


add_colors2meta = function(order_files, metadata){
  
  # ----------------------------------------------------------
  # 1. Читаем metadata
  # ----------------------------------------------------------
  
  meta = read.csv(
    metadata,
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 2. Проверяем ID
  # ----------------------------------------------------------
  
  if (!"ID" %in% colnames(meta)) {
    stop("В metadata нет столбца 'ID'")
  }
  
  
  # ----------------------------------------------------------
  # 3. Создаём нормализованный ID для объединения
  # ----------------------------------------------------------
  
  meta$ID_join = sub(
    "/.*$",
    "",
    meta$ID
  )
  
  meta$ID_join = gsub(
    "_",
    "-",
    meta$ID_join
  )
  
  
  # ----------------------------------------------------------
  # 4. Читаем список файлов с порядком таксонов
  # ----------------------------------------------------------
  
  order_files = read.table(
    order_files,
    stringsAsFactors = FALSE
  )$V1
  
  list_taxa_df = list()
  
  
  # ----------------------------------------------------------
  # 5. Читаем последовательности из каждой клады
  # ----------------------------------------------------------
  
  for (i in seq_along(order_files)) {
    
    taxa = read.table(
      order_files[i],
      stringsAsFactors = FALSE
    )$V1
    
    # ID до "/"
    taxa = sub(
      "/.*$",
      "",
      taxa
    )
    
    # NC_005790 -> NC-005790
    taxa = gsub(
      "_",
      "-",
      taxa
    )
    
    list_taxa_df[[i]] = taxa
  }
  
  
  # ----------------------------------------------------------
  # 6. Создаём цвета для клад
  # ----------------------------------------------------------
  
  colors = distinctColorPalette(
    length(list_taxa_df)
  )
  
  print(colors)
  
  
  # ----------------------------------------------------------
  # 7. Создаём градиент для каждой клады
  # ----------------------------------------------------------
  
  for (i in seq_along(list_taxa_df)) {
    
    num_colors = length(
      list_taxa_df[[i]]
    )
    
    print(i)
    print(colors[i])
    
    color1 = lighten(
      colors[i],
      0.4
    )
    
    color2 = darken(
      colors[i],
      0.4
    )
    
    print(color1)
    print(color2)
    
    color_range = colorRampPalette(
      c(color1, color2)
    )
    
    colors_clade = color_range(
      num_colors
    )
    
    
    # Таблица ID + цвет
    df = data.frame(
      ID_join = list_taxa_df[[i]],
      color = colors_clade,
      stringsAsFactors = FALSE
    )
    
    list_taxa_df[[i]] = df
  }
  
  
  # ----------------------------------------------------------
  # 8. Объединяем все клады
  # ----------------------------------------------------------
  
  df_colors = bind_rows(
    list_taxa_df
  )
  
  
  # ----------------------------------------------------------
  # 9. Проверяем ID, которые не нашлись
  # ----------------------------------------------------------
  
  missing = setdiff(
    df_colors$ID_join,
    meta$ID_join
  )
  
  if (length(missing) > 0) {
    
    warning(
      paste(
        "Эти ID есть в деревьях, но отсутствуют в metadata:",
        paste(missing, collapse = ", ")
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # 10. Добавляем цвета
  # ----------------------------------------------------------
  
  meta_upd = full_join(
    meta,
    df_colors,
    by = "ID_join"
  )
  
  
  # ----------------------------------------------------------
  # 11. Удаляем технический столбец
  # ----------------------------------------------------------
  
  meta_upd$ID_join = NULL
  
  
  return(meta_upd)
}