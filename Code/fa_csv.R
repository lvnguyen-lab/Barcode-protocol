#To create .fa file from .csv file

args = commandArgs(trailingOnly = T)
if (length(args) != 3) {
  message("\n\tUsage: fa_csv.R <working directory> <PATH to CSV> <output path and filename>\n\n")
  stop(1)
}

wd = args[1]
input = args[2]
output = args[3]

setwd(wd)
full = read.csv(paste0(wd, "/", input), head = TRUE)
full <- full$var_seq
bcs <- paste0(substr(full, 5, 6), substr(full, 10, 11), substr(full, 15, 16), substr(full, 20, 21), substr(full, 25, 26), substr(full, 30, 31))
bc_pref_str = c("GFP_")
#Want to have the full barcode in the ids
ids = paste0(">", bc_pref_str, full)
full <- as.character(full)
ids <- as.character(ids)
o = c()
o[seq(1, len=length(ids), by=2)] = ids
o[seq(2, len=length(full), by=2)] = full
write.table(o, paste0(wd, "/", output), row.names=F, col.names=F, quote=F)