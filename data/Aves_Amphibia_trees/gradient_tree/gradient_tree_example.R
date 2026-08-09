library(ggtree)
library(ggplot2)
library(phytools)
library(randomcoloR)
library(colorspace)
library(dplyr)

source("add_gradient_colors.R")

setwd(
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/Aves_Amphibia_trees/gradient_tree"
)

tree = read.nexus("B_ORF1a.nwk")
tree_rooted = midpoint.root(tree)
t = ggtree(tree_rooted, size = 0.75)


# Здесь указываешь номера внутренних узлов
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

order_files = c()

for (node in clades) {
  
  taxa = get_taxa_name(t, node)
  
  file_name = paste0(
    "tree_order_",
    node,
    ".csv"
  )
  
  write.table(
    taxa,
    file = file_name,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
  
  order_files = c(
    order_files,
    file_name
  )
}


writeLines(
  order_files,
  "order_files.txt"
)



info_upd = add_colors2meta(
  "order_files.txt",
  "Astroviridae_Aves_Amphibia_14072026.csv"
)


write.csv(
  info_upd,
  "metadata_upd.csv",
  row.names = FALSE
)

print("Готово: metadata_upd.csv создан")



