#Levenshtein.R
#Combines reads within 1 edit distance from each other

suppressWarnings(library(stringdist, quietly = T))

args = commandArgs(trailingOnly = T)

input = args[1]
output = args[2]

var_seq = list()
n_reads = list()
n = list()
Final = NULL
x = read.csv(input, head = FALSE)
x <- as.data.frame(table(x))
x <- x[order(x$Freq, decreasing = TRUE),]
colnames(x) <- c("var_seq", "n_reads")
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
final <- do.call(rbind, Map(data.frame, var_seq = var_seq, n_reads = n_reads))
write.csv(final, output, col.names = TRUE, row.names = FALSE)