diff_tb = read.table("data/E-GEOD-101794-analytics.tsv", sep = "\t", quote = "", header = TRUE, row.names = 1)[, 1:3]
colnames(diff_tb) = c("Symbol", "p_value", "log2foldChange")
diff_tb$p_adjust = p.adjust(diff_tb$p_value, "BH")
diff_tb$p_adjust[is.na(diff_tb$p_adjust)] = Inf
load("data/E-GEOD-101794-atlasExperimentSummary.Rdata")
counts = experiment_summary@listData$rnaseq@assays$data@listData$counts
coldata = experiment_summary@listData$rnaseq@colData
l = coldata[["clinical_information"]] == "A1a Paris classification"
counts = counts[, l]
coldata = coldata[l, ]
saveRDS(lt, file = "data/EBI_E-GEOD-101794_expression.rds")
