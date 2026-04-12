# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#             Chloris             # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
library(optparse)
parser = OptionParser(description = 'Execute Chloris')
parser = add_option(parser, '--rootDir', default = "/your/local/path/to/repo/", type = "character")
parser = add_option(parser, '--sampleID', default = "toy_data", type = "character")
parser = add_option(parser, '--refPath', default = "", type = "character")
args <- parse_args(parser)

library(VariantAnnotation)
library(GenomicFeatures)
library(GenomicRanges)
library(Matrix)
library(Seurat)
library(Chloris)
library(Matrix.utils)
source(system.file("BAF_window.R", package = "Chloris"))
source(system.file("infercnv_utils.R", package = "Chloris"))

root_Dir <- args$rootDir
sample_id <-args$sampleID
ref_Path <- args$refPath
input_Dir <- paste0(root_Dir, "Real_data/", sample_id, "/filtered_feature_bc_matrix")
Out_Dir <- paste0(root_Dir, "Result/", sample_id, "/Chloris/")

dir.create(Out_Dir, showWarnings = FALSE, recursive = TRUE)
setwd(Out_Dir)
message("Method: Chloris ...")
message("Process ", sample_id, " ...")
message("Output Directory: ", Out_Dir)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#           Input data            # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# (1) preprocess SNP data
vcf <- readVcf(paste0(Out_Dir, "cellSNP.base.vcf"), "hg38")
snp_gr <- rowRanges(vcf)
gtf <- makeTxDbFromGFF(paste0(root_Dir, "resources/Chloris/genes.gtf.gz"), format="gtf")
gene_gr <- genes(gtf)

# Overlap SNPs with genes
hits <- findOverlaps(snp_gr, gene_gr)
snp2gene <- data.frame(snp_index = queryHits(hits),
                       gene_id = names(gene_gr[subjectHits(hits)]))
# gene id to symbol
geneinfo <- read.table(paste0(root_Dir, "resources/Chloris/geneInfo.tab"), 
                       sep = '\t', skip = 1, header = F)
geneinfo <- geneinfo[,1:2]
colnames(geneinfo) <- c("gene_id", "symbol")
snp2gene <- merge(snp2gene, geneinfo, sort = F)
# load snp matrix
AD.mtx <- Matrix::readMM(paste0(Out_Dir, "cellSNP.tag.AD.mtx")) # snp x cell
DP.mtx <- Matrix::readMM(paste0(Out_Dir, "cellSNP.tag.DP.mtx")) # snp x cell

# Create a mapping from SNP to gene
# AD: extract only the rows we need and duplicate
AD.mtx <- AD.mtx[snp2gene$snp_index, , drop = FALSE]
DP.mtx <- DP.mtx[snp2gene$snp_index, , drop = FALSE]

# Update rownames to gene symbols (used as groupings)
rownames(AD.mtx) <- snp2gene$symbol
rownames(DP.mtx) <- snp2gene$symbol
AD.gene.mtx <- aggregate.Matrix(AD.mtx, groupings = rownames(AD.mtx), fun = "sum")
DP.gene.mtx <- aggregate.Matrix(DP.mtx, groupings = rownames(DP.mtx), fun = "sum")
# update cell ID
cell_info <- read.table(paste0(Out_Dir, "cellSNP.samples.tsv"))
colnames(AD.gene.mtx) <- colnames(DP.gene.mtx) <- cell_info$V1
save(AD.gene.mtx, DP.gene.mtx, file = paste0(Out_Dir, '/BAF.RData'))

# (2) expression data
input_expr <- as.matrix(Read10X(data.dir = input_Dir))
input_expr <- input_expr[,colnames(DP.gene.mtx)] # keep cell order consistent w/ A & D matrix
ref_expr <- as.matrix(read.csv(ref_Path, row.names = 1))
colnames(ref_expr) <- paste0("Ref_", colnames(ref_expr))

overlapped.genes <- intersect(rownames(input_expr), rownames(ref_expr))
overlapped.genes <- intersect(overlapped.genes, rownames(AD.gene.mtx)) # overlapped with A & D matrix
input_expr <- input_expr[overlapped.genes,]
ref_expr <- ref_expr[overlapped.genes,]
normal_cell_names <- colnames(ref_expr)
input_expr <- cbind(input_expr, ref_expr)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#            Execute              # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
## (0) sort gene order
# sort gtf file to get gene order
gene_gr <- sort(gene_gr)
gene_order <- data.frame(gene_id = gene_gr$gene_id,
                         chr = as.character(gene_gr@seqnames))
chr_keep <- paste0("chr", 1:22)
gene_order <- gene_order[gene_order$chr %in% chr_keep,] # keep chr1-22
gene_order$chr <- as.numeric(gsub("^chr", "", gene_order$chr))
gene_order <- merge(gene_order, geneinfo, sort = FALSE)
gene_order <- gene_order[gene_order$symbol %in% overlapped.genes,] # keep overlapped genes
gene_order <- gene_order[!duplicated(gene_order$symbol),]
sorted_genes <- gene_order$symbol
# update gene order
AD.gene.mtx <- as.matrix(AD.gene.mtx[sorted_genes,]) # A & D matrix shouldn't be dgCMatrix
DP.gene.mtx <- as.matrix(DP.gene.mtx[sorted_genes,])
input_expr <- input_expr[sorted_genes,] 

## (1) Identify normal cells via SNP
## Window smooth BAF to obercome sparsity
zero.col <- colSums(DP.gene.mtx) == 0 # BAF_window function can't process D matrix with empty column (cell)
AD.gene.mtx <- AD.gene.mtx[,!zero.col] 
DP.gene.mtx <- DP.gene.mtx[,!zero.col]
input_expr <- input_expr[,!zero.col]

AD.gene.mtx[nrow(AD.gene.mtx),] <- 0; DP.gene.mtx[nrow(DP.gene.mtx),] <- 0 # last snp can only be zero
BAF <- BAF_window(AD.gene.mtx, DP.gene.mtx, window_size = 20) 
save(BAF, file = paste0(Out_Dir, 'smoothed_BAF.RData'))
get_norm <- Chloris(A = BAF$A, D = BAF$D, S = 2, init = "random", break_idx = NULL,  
                    burnin_tol = 100, Gibbs_tol = 100, 
                    cluster_shrink_tol = 10, min_cluster_size = 2, verbose = FALSE)

is_normal_cell <- get_norm$cluster_est == which.min(apply(get_norm$state_est, 1, function(x) sum(x == 1)))
tumor_normal_df <- data.frame(cell.names = colnames(AD.gene.mtx),
                              Chloris.pred = ifelse(is_normal_cell, "normal", "tumor"))
write.table(tumor_normal_df, file = paste0(Out_Dir, "tumor_cell_pred.txt"),
            row.names = F, quote = F)

## (2) Normalization using normal panel
cs <- colSums(input_expr)
tmp_data <- sweep(input_expr, STATS = cs, MARGIN = 2, FUN="/")  
normalize_factor = median(cs)  # normalize_factor
normalized <- tmp_data * normalize_factor

## Get RDR from gene expr
normalized <- log2(normalized + 1)
# use predefined normal reference
ref <- rowMeans(normalized[,colnames(normalized) %in% normal_cell_names])  
RDR <- normalized - ref
RDR <- RDR[,colnames(BAF$A)] # remove external normal cells (to make dim(RDR) == dim(A) == dim(D))
RDR <- apply(RDR, 2, smooth_helper, window_length = 51) ## Window smoothing

## Recentering
col_median <- apply(RDR, 2, function(x) { median(x, na.rm = TRUE) } )
RDR <- t(apply(RDR, 1, "-", col_median))


## Get index of the first gene in each chromosome
chrs <- gene_order$chr
break_idx <- which(chrs[-1]  - chrs[-length(chrs)] > 0) + 1

## Run full model
res <- Chloris(RDR, A = BAF$A, D = BAF$D, 
               break_idx = break_idx, burnin_tol = 300, Gibbs_tol = 300,
               cluster_shrink_tol = 10, min_cluster_size = 2, verbose = TRUE)
names(res$cluster_est) <- colnames(RDR) # add cell ID
colnames(res$state_est) <- rownames(RDR) # add gene symbol
save(res, file = paste0(Out_Dir, "Chloris_res.RData"))

# subclone: res$cluster_est
# CNV matrix (subclones x genes): res$cluster_est