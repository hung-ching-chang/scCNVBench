# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#             CONICSmat           # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
library(optparse)
parser = OptionParser(description = 'Execute CONICSmat')
parser = add_option(parser, '--rootDir', default = "/your/local/path/to/repo/", type = "character")
parser = add_option(parser, '--sampleID', default = "toy_data", type = "character")
parser = add_option(parser, '--refPath', default = "", type = "character")
args <- parse_args(parser)

library(CONICSmat)
library(Seurat)
root_Dir <- args$rootDir
sample_id <-args$sampleID
ref_Path <- args$refPath
input_Dir <- paste0(root_Dir, "Real_data/", sample_id, "/filtered_feature_bc_matrix")
Out_Dir <- paste0(root_Dir, "Result/", sample_id, "/CONICSmat/")

dir.create(Out_Dir, showWarnings = FALSE, recursive = TRUE)
setwd(Out_Dir)
message("Method: CONICSmat ...")
message("Process ", sample_id, " ...")
message("Output Directory: ", Out_Dir)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#           Input data            # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
input_expr <- as.matrix(Read10X(data.dir = input_Dir, strip.suffix = TRUE))
ref_expr <- as.matrix(read.csv(ref_Path, row.names = 1))
colnames(ref_expr) <- paste0("Ref_", colnames(ref_expr))

overlapped.genes <- intersect(rownames(input_expr), rownames(ref_expr))
input_expr <- cbind(input_expr[overlapped.genes,], ref_expr[overlapped.genes,])
normal_cell_names <- colnames(ref_expr)

# coordinates of the human chromosome arms
regions <- read.table(paste0(root_Dir, "resources/CONICSmat_ref/chromosome_arm_positions_grch38.txt"),
                            sep = "\t",row.names = 1,header = T)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#            Execute              # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# get gene position
# gene_pos <- getGenePositions(rownames(input_expr))
gene_pos <- read.table(paste0(root_Dir, "resources/CONICSmat_ref/gene_position.txt"))

# log2(CPM/10+1)
sdepth <- colSums(input_expr)
input_expr <- t(t(input_expr)/sdepth*1000000)
input_expr <- log2(input_expr/10+1)
input_expr[which(is.na(input_expr))] <- 0

#filter uninformative genes (expressed in few cells)
num.prefilter.cells <- nrow(input_expr)
input_expr <- filterMatrix(input_expr, gene_pos[,"hgnc_symbol"], minCells = 5)
message("filter ", num.prefilter.cells - nrow(input_expr), " out of ", num.prefilter.cells, " genes")

#calculate a normalization factor for each cell
normFactor <- calcNormFactors(input_expr)

# identify cancer cells
l <- plotAll(input_expr, normFactor, regions, gene_pos, 
             normal = which(colnames(input_expr) %in% normal_cell_names),
             fname = paste0(Out_Dir, "CONICSmat_CNV"))

# to obtain a final assignment as malignant or non-malignant cells
lrbic <- read.table(paste0(Out_Dir, "CONICSmat_CNV_BIC_LR.txt"),
                    sep = "\t", header = T, row.names = 1, check.names = F)
candRegions <- rownames(lrbic)[which(lrbic[,"BIC difference"] > 200 & lrbic[,"LRT adj. p-val"] < 0.01)]


# #generate a heatmap of posterior probabilities for candidate Regions
pdf(paste0(Out_Dir, "post_heatmap.pdf"))
hi=plotHistogram(l[,candRegions],input_expr,clusters=2,zscoreThreshold=2) # adjust threshold to obtain a clear cluster separation (e.g., threshold = 1 or 2)
dev.off()

#assign labels to subsets (normal cluster should have more reference cells)
ref_cell_group <- names(which.max(table(hi[normal_cell_names])))
normal <- which(hi==as.integer(ref_cell_group)) # CONICSmat wiki instructs to use "hi==1", but the cluster 1 is not always be normal group
tumor <- which(hi!=as.integer(ref_cell_group))
tumor_normal_df <- data.frame(cell.names = names(hi),
                              CONICSmat.pred = ifelse(hi == as.integer(ref_cell_group), "normal", "tumor"))
write.table(tumor_normal_df, file = paste0(Out_Dir, "tumor_cell_pred.txt"),
            row.names = F, quote = F)

# calculate posterior
post_res <- plotAll(input_expr, normFactor, regions[candRegions,], gene_pos,
                    paste0(Out_Dir, "SUVA_CNVs_with_info"), 
                    normal = normal, tumor = tumor)
saveRDS(post_res, file = paste0(Out_Dir, "final_posterior.RDS"))