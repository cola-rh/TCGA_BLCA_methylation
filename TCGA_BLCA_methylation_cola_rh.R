setwd("/omics/groups/OE0246/internal/guz/cola_hc/examples/TCGA_BLCA_methylation")
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.BLCA.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_BLCA_methylation_cola_rh.rds")

cola_report(rh, output = "TCGA_BLCA_methylation_cola_rh_report", title = "cola Report for Hierarchical Partitioning - 'TCGA_BLCA_methylation'")
