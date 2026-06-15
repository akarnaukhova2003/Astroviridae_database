library(ggplot2)

lines <- readLines("/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/Aves_full_seq_1a.fasta")

seqs <- list()
current_seq <- ""

for (line in lines) {
  if (startsWith(line, ">")) {
    if (nchar(current_seq) > 0) {
      seqs <- c(seqs, current_seq)
    }
    current_seq <- ""
  } else {
    current_seq <- paste0(current_seq, line)
  }
}
if (nchar(current_seq) > 0) {
  seqs <- c(seqs, current_seq)
}

seqs <- seqs[seqs != ""]

lengths <- nchar(seqs)
df <- data.frame(lengths = lengths)

# -------------------------
# 3. Параметры оси
# -------------------------
bin_width <- 100

#min_len <- min(df$lengths, na.rm = TRUE)
min_len <- 0
max_len <- 3850
#max_len <- max(df$lengths, na.rm = TRUE)

x_start <- floor(min_len / 1000) * 1000
x_end   <- ceiling(max_len / 1000) * 1000

# -------------------------
# 4. График
# -------------------------
p <- ggplot(df, aes(x = lengths)) +
  geom_histogram(
    binwidth = bin_width,
    fill = "#5DADE2",
    color = "black",
    boundary = 0,
    na.rm = TRUE
  ) +
  scale_x_continuous(
    limits = c(min_len, max_len),
    breaks = seq(x_start, x_end, by = 500),
    minor_breaks = seq(x_start, x_end, by = 500),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title = "Распределение длин ORF1a последовательностей клады птиц",
    x = "Длина (п.н.)",
    y = "Количество"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_line(color = "grey90"),
    axis.line.x = element_line(color = "grey90", linewidth = 0.2),
    plot.margin = margin(10, 10, 15, 10) 
  )
p
ggsave(
  filename = "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/clade_analysis/AVES/dist_ORF1a_AVES.png",
  plot = p,
  width = 2000,
  height = 1600,
  units = "px",
  dpi = 300
)

