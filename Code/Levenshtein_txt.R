#Levenshtein.R
# For use with DNA amplicon sequencing data, input is txt file

suppressWarnings(library(stringdist, quietly = T))

args = commandArgs(trailingOnly = T)

input = args[1]
output_base = args[2]

start_time <- Sys.time()

var_seq = list()
n_reads = list()
n = list()
Final = NULL

if (file.size(input) > 0) {
  print("Reading Grep file")
  x = read.delim(input, head = FALSE)
  x <- as.data.frame(table(x))
  x <- x[order(x$Freq, decreasing = TRUE),]
  colnames(x) <- c("var_seq", "n_reads")
  print("Writing pre-Levenshtein csv")
  write.csv(x, paste0(output_base, "_pre-Levenshtein.csv"), col.names = TRUE, row.names = FALSE)
  print("Running stringdist")
  repeat {
    for (j in 1:nrow(x)) {
      a <- stringdist(x$var_seq[1], x$var_seq[j])
      n = rbind(n, a)
    }
    x$n <- n
    x_merge <- x[x$n <= 1,]
    var_seq <- append(var_seq, x_merge$var_seq[1])
    n_reads <- append(n_reads, sum(x_merge$n_reads))
    x <- x[x$n > 1,]
    n = list()
    if (nrow(x) == 0) {
      break
    }
  }
  print("Finished stringdist")
  final <- do.call(rbind, Map(data.frame, var_seq = var_seq, n_reads = n_reads))
  print("Writing final csv")
  write.csv(final, paste0(output_base, "_post-Levenshtein.csv"), col.names = TRUE, row.names = FALSE)
  print(paste0("--- ", (Sys.time() - start_time), " seconds ---" ))
} else {
  print("Input file was empty")
}