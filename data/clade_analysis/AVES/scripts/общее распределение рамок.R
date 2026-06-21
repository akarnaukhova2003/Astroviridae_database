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

len_1a <- get_lengths("/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Aves_full_seq_1a.fasta")
len_1b <- get_lengths("/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Aves_full_seq_1b.fasta")
len_2  <- get_lengths("/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Aves_full_seq_2.fasta")

df_lengths <- data.frame(
  length = c(len_1a, len_1b, len_2),
  ORF = c(
    rep("ORF1A", length(len_1a)),
    rep("ORF1B", length(len_1b)),
    rep("ORF2",  length(len_2))
  )
)


p <- ggplot(df_lengths, aes(x = length)) +
  geom_histogram(
    bins = 50,
    fill = "#F8BBD0",
    color = "#AD1457",
    alpha = 0.8
  ) +
  facet_wrap(~ORF, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Распределение длин ОРС клады птиц, рептилий и амфибий",
    x = "Длина (п.н.)",
    y = "Количество"
  ) +
  theme(
    axis.title.x = element_text(size=16, face="bold", margin=margin(t=8)),
    axis.title.y = element_text(size=16, face="bold", margin=margin(r=2)),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    strip.text = element_text(face="bold", size=12),
    strip.background = element_rect(fill = "#F8BBD0", color = "black"),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white"),
    plot.title = element_text(size=16,hjust = 0.5, face = "bold")
  )

print(p)

ggsave(
  filename = "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/dist_lengths_ORFs.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

