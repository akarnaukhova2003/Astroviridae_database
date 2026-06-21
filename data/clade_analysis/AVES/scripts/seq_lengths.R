lines <- readLines("/Users/abagavetdinova/Desktop/lab/Astroviridae_database/data/Astroviridae_15102025_1B.fasta")
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
seqs <- c(seqs, current_seq)

lengths <- nchar(seqs)
hist(lengths,
     breaks = 200,
     xaxt = "n",
     col = "skyblue",
     main = "Распределение длин последовательностей ORF1b",
     xlab = "Длина (п.н)",
     ylab = "Количество", 
     cex.main = 1.6,
     cex.lab = 1.4,
     cex.axis = 1.2)

ticks <- seq(0, max(lengths), by = 500)

axis(1, at = ticks, labels = ticks)

