# Associate 10X cell barcodes with the whitelisted barcodes. Also report on the GFP expression of each cell.

require(Seurat, quietly = T)
require(Rsamtools, quietly = T)
require(data.table, quietly = T)
require(dplyr, quietly = T)

args = commandArgs(trailingOnly = T)
if (length(args) != 4) {
  message("\n\tUsage: associate_barcodes_with_10x_cells.R <sample name> <10x base dir> <barcodes with mapped reads file name> <output dir>\n\n")
  stop(1)
}

sample_name = args[1]
TENx_base_dir = args[2]
barcodes_sam_fn = args[3]
odir = args[4]

dir.create(odir, showWarnings = F, recursive = T)

scdb_init(sprintf("%s/scdb", TENx_base_dir))
m = scdb_mat(sample_name)

gfp = m@mat['eGFP', ]

lib_bc_reads = data.table::fread(barcodes_sam_fn, select=c(1,3), col.names=c('qname', 'lib_bc'))

filt_bam_fn = tempfile(fileext=".bam")
x = filterBam(sprintf("%s/cellranger/%s/possorted_genome_bam.bam", TENx_base_dir, sample_name), param=ScanBamParam(what='qname', tag=c('CB', 'UB')), filter=FilterRules(list(by_qname = function(x) x$qname %in% lib_bc_reads$qname)), destination=filt_bam_fn)
indexBam(filt_bam_fn)
xf = scanBam(filt_bam_fn, param=ScanBamParam(what='qname', tag=c('CB', 'UB')))
xf_df_u = data.frame(qname=xf[[1]]$qname, cb=xf[[1]]$tag$CB, ub=xf[[1]]$tag$UB) %>% 
  left_join(lib_bc_reads) %>%
  select(-qname) %>%
  distinct

cell_bc_t = table(xf_df_u$cb, xf_df_u$lib_bc)

n_gfp_with_barcodes = length(intersect(rownames(cell_bc_t), names(which(gfp > 0))))
message(sprintf("\n%s has %d GFP positive cell out of %d cells (%.2f).\nFound barcodes for %d of them (%.2f).\nGFP negative cells with detected barcodes: %d\n", sample_name, sum(gfp > 0), length(gfp), mean(gfp > 0), n_gfp_with_barcodes, n_gfp_with_barcodes / sum(gfp > 0), nrow(cell_bc_t) - n_gfp_with_barcodes))

write.csv(cell_bc_t, sprintf("%s/%s_cell_bc_umis_count.csv", odir, sample_name))
write.csv(as.data.frame(gfp), sprintf("%s/%s_cell_GFP_umis_count.csv", odir, sample_name))