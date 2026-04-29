library(ape)
library(ggplot2)
library(dplyr)

get_lengths <- function(fasta_path) {
  aln <- read.dna(fasta_path, format = "fasta")
  seqs_list <- as.list(aln)
  
  lengths <- sapply(seqs_list, function(seq) {
    seq_char <- as.character(seq)
    sum(seq_char != "-" & !is.na(seq_char))
  })
  
  return(as.numeric(lengths))
}

len_1a <- get_lengths("/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/Astroviridae_15102025_1B.fasta")

df_lengths <- data.frame(
  length = len_1a,
  ORF = rep("ORF1A", length(len_1a))
)

p <- ggplot(df_lengths, aes(x = length)) +
  geom_histogram(
    bins = 50,
    fill = "#F8BBD0",
    color = "#AD1457",
    alpha = 0.8
  ) +
  theme_bw() +
  labs(
    title = "Распределение длин ORF клады птиц",
    x = "Длина (п.н.)",
    y = "Количество"
  ) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 2)),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "#F8BBD0", color = "black"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p)

ggsave(
  filename = "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/dist_lengths_ORF1b.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)