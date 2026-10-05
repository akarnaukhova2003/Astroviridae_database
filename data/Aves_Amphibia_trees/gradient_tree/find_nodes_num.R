library(ape)
library(ggtree)
library(ggplot2)

setwd(
  "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/snakemake/results/trees/"
)

reference_tree_file = "B/B_ORF1a.nex"
reference_tree_file
outgroup_id = "'PQ893528/NA/NSW-2021/Eolophus-roseicapilla/2021'"

reference_tree = read.nexus(
  reference_tree_file
)

reference_tree$tip.label = trimws(reference_tree$tip.label)

outgroup_index = which(
  reference_tree$tip.label == outgroup_id
)

if (length(outgroup_index) == 0) {
  stop("Outgroup не найдена в дереве")
}

reference_tree = root(
  reference_tree,
  outgroup = outgroup_index,
  resolve.root = TRUE
)

reference_tree$node.label = as.character(
  (Ntip(reference_tree) + 1):(Ntip(reference_tree) + Nnode(reference_tree))
)

reference_plot = ggtree(
  reference_tree,
  size = 0.75
)

node_plot =
  reference_plot +
  geom_text2(
    aes(
      subset = !isTip,
      label = label
    ),
    size = 4,
    color = "red"
  ) +
  geom_tiplab(
    size = 3
  )

ggsave(
  "B_ORF1a_nodes.pdf",
  node_plot,
  width = 40,
  height = 30,
  dpi = 300
)



