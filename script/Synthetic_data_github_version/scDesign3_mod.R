# -----------------------------------------------------------------------------
# scDesign3_mod.R — CNV-aware synthetic gene-expression data generation
#
# Adapted from the scDesign3 workflow by Dongyuan Song and collaborators:
# https://github.com/SONGDONGYUAN1994/scDesign3
# Reference: Song D, et al. (2024). scDesign3 generates realistic in silico data
# for multimodal single-cell and spatial omics. Nature Biotechnology 42:247–252.
# https://doi.org/10.1038/s41587-023-01772-1
#
# Study-specific modifications: Seurat-based reference-cell clustering,
# subclone-specific CNV dosage/signal factors, CNV-adjusted simulation,
# output files, and heatmap visualization. See MODIFICATIONS.md for details.
#
# Usage: Set scDesign3_path and input/output paths below, then run this file
# with Rscript and the --alpha_amp, --beta_amp, --alpha_del, --beta_del options.
# Supporting scDesign3 R scripts are sourced automatically; do not run them
# individually. Requires the dependencies loaded below.
#
# Original scDesign3 code is distributed under the MIT license; retain the
# upstream LICENSE and copyright notices when redistributing derived code.
# -----------------------------------------------------------------------------

library(Seurat)
library(SingleCellExperiment)
library(parallel)
scDesign3_path <- "path/to/your/downloaded/Synthetic_data_github_version/folder"
source(paste0(scDesign3_path, "/construct_data.R"))
source(paste0(scDesign3_path, "/data.R"))
source(paste0(scDesign3_path, "/extract_para.R"))
source(paste0(scDesign3_path, "/fit_copula.R"))
source(paste0(scDesign3_path, "/fit_marginal.R"))
source(paste0(scDesign3_path, "/gamlss_fix.R"))
source(paste0(scDesign3_path, "/perform_lrt.R"))
source(paste0(scDesign3_path, "/plot_reduceddim.R"))
source(paste0(scDesign3_path, "/scdesign3.R"))
source(paste0(scDesign3_path, "/simu_new.R"))
source(paste0(scDesign3_path, "/sparse_cov.R"))
library(optparse)


#############################
# Define preprocessing and synthetic data generating functions
get_cluster_seurat <- function(data_mat, res = 0.5){
  submat <- data_mat
  if(is.null(rownames(submat)))
    rownames(submat) <- 1:nrow(submat)
  count_seurat <- CreateSeuratObject(counts = submat)
  ### normalization
  count_seurat <- NormalizeData(count_seurat)
  ### select highly variable genes
  count_seurat <- FindVariableFeatures(count_seurat, selection.method = "vst", nfeatures = 4000)
  ### scale the data
  count_seurat <- ScaleData(count_seurat)
  ### PCA
  count_seurat <- RunPCA(count_seurat,
                         features = VariableFeatures(object = count_seurat))
  ### clustering
  dims = 1:min(10, min(ncol(data_mat), nrow(data_mat))-1)
  count_seurat <- FindNeighbors(count_seurat, dims = dims)
  
  cat("Performing Louvain clustering with resolution 0.5...\n")        
  count_seurat <- FindClusters(count_seurat, resolution = res)
  ### results
  cluster_predicted <- as.integer(Idents(count_seurat))
  colnames(submat) <- as.character(cluster_predicted)
  cat(sprintf("%s cell clusters generated.\n", length(unique(cluster_predicted))))
  
  return(cluster_predicted)
}

scDesign3_SyntheticCount <- function(count_matrix, output_file, cell_label = NULL, 
                                     cnv_factor_mat = NULL) {
  # Add garbage collection
  gc()
  
  matrix_num <- data.matrix(count_matrix[,2:ncol(count_matrix)])
  matrix_num <- apply(matrix_num, 2, as.numeric)
  count_pergene_vec <- rowSums(matrix_num)
  matrix_num_filtered <- matrix_num[count_pergene_vec>0,]

  # Add progress messages
  cat("Processing matrix of size:", dim(matrix_num_filtered), "\n")
  
  # Annotate with given cell types; Use cluster labels as cell types if not given
  cluster_label <- get_cluster_seurat(matrix_num_filtered)
  if(is.null(cell_label)){
    cell_type <- cluster_label 
  }else{
    cell_type <- cell_label 
  }
  meta_tbl <- data.frame(cluster_label = as.factor(cluster_label),
                         cell_type = as.factor(cell_type),
                         corr_group_label = as.factor(cluster_label))

  # Clear unnecessary objects
  #rm(matrix_num_filtered)
  gc()
  
  sct_dat <- SingleCellExperiment(list(counts = matrix_num), colData = meta_tbl)
  rownames(sct_dat) <- 1:nrow(sct_dat)
  
  if(length(unique(cell_label))==1){
    mu_formula <- "1"
  }else if(is.null(cell_label)){
    mu_formula <- "cluster_label"
  }else{
    mu_formula <- "cell_type"
  }
  
  # Modify scdesign3 parameters to use less memory
  simu_res <- scdesign3(
    sce = sct_dat,
    assay_use = "counts",
    celltype = "cell_type",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = NULL,
    mu_formula = mu_formula,
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 5,  
    usebam = FALSE,
    corr_formula = "corr_group_label",
    copula = "gaussian",
    DT = TRUE,
    pseudo_obs = FALSE,
    return_model = FALSE,
    nonzerovar = FALSE,
    parallelization = "mcmapply",
    cnv_factor_mat = cnv_factor_mat)
  
  # Extract and save the simulated count matrix
  simu_counts <- as.matrix(simu_res$new_count)
  
  print("Simulated count matrix generated!")
  
  rownames(simu_counts) <- count_matrix[,1]
  #colnames(simu_counts) <- colnames(count_matrix)[-1]

  # Save the cell types/clusters  
  write.table(simu_res$new_covariate, 
              file = paste0(output_file, ".scDesign3Simulated.CellTypeLabel.txt"), 
              row.names = FALSE, 
              col.names = FALSE)

  # Save the simulated counts
  write.table(simu_counts, 
              file = paste0(output_file, ".scDesign3Simulated.txt"), 
              sep = "\t", 
              row.names = TRUE, 
              col.names = TRUE)
  
  cat(sprintf("Saved simulated counts to: %s\n", output_file))
  cat("[scReadSim] Done.\n")
  
  # Final cleanup
  rm(simu_counts)
  gc()
  
}

cnv_signal_strength_mat <- function(count_mat, threshold = 0.25, dna_dosage_mat,
                                    subclone_prop, 
                                    shape.amp, scale.amp, 
                                    shape.del, scale.del){
  ncell <- dim(count_mat)[2]-1
  ngene <- dim(count_mat)[1]
  nsubclone <- dim(dna_dosage_mat)[2]
  if (sum(subclone_prop) != 1) stop("subclone proportions must sum to 1!")
  if (nsubclone != length(subclone_prop)) stop("number of subclones must match!")
  subclone_prop_cumulative <- cumsum(subclone_prop)
  subclone_boudary_indices <- floor(ncell * subclone_prop_cumulative)
  # Initiate a CNV signal strength matrix with all ones
  res <- matrix(1, ngene, ncell)
  rownames(res) <- count_mat[,1]
  print("Begin preparing cnv factors")
  for (i in 1:nsubclone) {
    # Find genes that has expression large/small enough to be suspected as CNV behavior
    # in subclone i
    cnv_genes <- dna_dosage_mat[(dna_dosage_mat[,i] > 1+threshold | 
                                  dna_dosage_mat[,i] < 1-threshold), i, drop=FALSE]
    common_genes <- intersect(rownames(cnv_genes), count_mat[,1])
    if (is.null(common_genes)) stop(paste0("C", as.character(i), " no cnv genes found in count matrix"))
    if (i==1) {
      cell_indices <- 1 : subclone_boudary_indices[i]
    } else {
      cell_indices <- (subclone_boudary_indices[i-1]+1) : subclone_boudary_indices[i]
    }
    ncommongene <- length(common_genes)
    nsubclonecell <- length(cell_indices)
    # Crafting Dg matrix for subclone i
    common_genes_cnv_mat <- matrix(rep(cnv_genes[common_genes, ], nsubclonecell), 
                            nrow = ncommongene, 
                            ncol = nsubclonecell)
    # Initiate a 2 random matrices, one for amplification Gamma, another for deletion Gamma
  
    gamma_mat_amp <- matrix(rgamma(ncommongene * nsubclonecell, shape = shape.amp, scale = scale.amp), 
                        nrow = ncommongene, ncol = nsubclonecell)
    gamma_mat_del <- matrix(rgamma(ncommongene * nsubclonecell, shape = shape.del, scale = scale.del), 
                            nrow = ncommongene, ncol = nsubclonecell)
    # Use amplification Gamma if Dg is amplification there (>1), use deletion Gamma
    # if Dg is deletion there (<1)
    gamma_mat <- gamma_mat_amp
    gamma_mat[common_genes_cnv_mat < 1] <- gamma_mat_del[common_genes_cnv_mat < 1]
    # We found log(1+Cg) ~ gamma(shape, scale), transferring back to calculate Cg
    # Some Cg get extremely large because exponential function will enlarge the scale, 
    # add a cap to large values of Cg
    Cg_mat <- exp(gamma_mat) - 1
    Cg_mat[Cg_mat>10] <- 10
    Cg_mat[Cg_mat<0.001] <- 0.001
    # (No need) For amplification, signal strength is multiplicative; For loss, signal effect is division
    # Cg_mat[common_genes_cnv_mat < 1] <- 1 / Cg_mat[common_genes_cnv_mat < 1]
    
    # Modify the factor of genes with CNV issues as dosage effect(Dg, common_genes_cnv_mat) times
    # signal strength (Cg_mat, Cg)
    print(dim(common_genes_cnv_mat))
    print(dim(Cg_mat))
    res[common_genes, cell_indices] <- common_genes_cnv_mat * Cg_mat
    print(paste0("c", as.character(i), " cnv factors OK!"))
  }
  t(res)
}

#############################
##### 4.1
# Call functions with input paths
sampleID = "SNU638"
filename = "path/to/output/directory"

# Load the real count matrix, removes genes with all 0 expression in all cells
normal_gene_count <- Read10X("path/to/downloaded/P5931_normal_ref_matrix/")
gene_count_mat <- as.matrix(normal_gene_count)
# Make gene names as the first column of count matrix
gene_count_mat <- data.frame(gene = rownames(gene_count_mat),
                             gene_count_mat,     
                             check.names = FALSE)
matrix_num <- data.matrix(gene_count_mat[,2:ncol(gene_count_mat)])
matrix_num <- apply(matrix_num, 2, as.numeric)
count_pergene_vec <- rowSums(matrix_num)
gene_count_mat <- gene_count_mat[count_pergene_vec>0,]

# gene_count_mat <- gene_count_mat[sample(nrow(gene_count_mat), 1000), c(1, sample(2:ncol(gene_count_mat), 500))]

# Input cell types. If using normal reference, all cells should have the same cell type
cellType <- rep("epithelial", times = dim(gene_count_mat)[2]-1)

#intergene_count_mat <- read.table(paste0(filename, ".intergene.countmatrix.txt"), sep="\t", header = FALSE)

##### 4.2
# Load segment-level DNA CNV dosage effect file
library(GenomicRanges)
library(dplyr)
# Subclone segment information file for SNU638
subclone_seg_cnv <- readRDS("path/to/downloaded/SNU638subclone_seg_cnv.rds")
gene.info <- readRDS("path/to/downloaded/gene_info.rds")
bin.info <- subclone_seg_cnv$bin.info
subclone.bin.cnv <- subclone_seg_cnv$subclone.cnv
# Finding overlaps
olaps <- suppressWarnings(GenomicRanges::findOverlaps(gene.info,
                                                      bin.info,
                                                      ignore.strand = TRUE))
# Create a data frame that will contain the genes and positions (index) in the pipeline ranges.
# some genes might overlap more than one range (more than one bin), in this case only one will be kept
gene.meta <- data.frame(
  gene = as.character(data.frame(gene.info)$symbol[queryHits(olaps)]),
  chr = as.character(data.frame(gene.info)$seqnames[queryHits(olaps)]),
  start = as.character(data.frame(gene.info)$start[queryHits(olaps)]),
  end = as.character(data.frame(gene.info)$end[queryHits(olaps)]),
  bin_pos = subjectHits(olaps)
) %>%
  dplyr::distinct(gene, .keep_all = TRUE)
gene.meta <- gene.meta[order(gene.meta$bin_pos, gene.meta$start),]
gene.meta <- gene.meta[!is.na( gene.meta$gene),] # remove missing gene symbol
# head(gene.meta, 10)

# (segment)bin-level to gene-level
subclone.gene.cnv <- subclone.bin.cnv[gene.meta$bin_pos,]
rownames(subclone.gene.cnv) <- gene.meta$gene
# head(subclone.gene.cnv)

# (Optional) Enlarge the deletion segment on chr14
# Make chr14:67040000 -chr14:107043718 all having deletion CNV
subclone.gene.cnv <- cbind(subclone.gene.cnv, gene.meta)
subclone.gene.cnv[(subclone.gene.cnv$chr=="chr14" & 
                     as.numeric(subclone.gene.cnv$start) >= 67040000 &
                     as.numeric(subclone.gene.cnv$start) < 107043718), c("c1", "c2", "c3")] <- rep(0.5032941, 3)
subclone.gene.cnv <- subclone.gene.cnv[, c("c1", "c2", "c3")]

##### 4.3
# Construct a CNV factor matrix that has the same dimension as the transpose 
# of count matrix
parser <- OptionParser()
parser <- add_option(parser, '--beta_amp', default = NA, type = "numeric")
parser <- add_option(parser, '--alpha_amp', default = 1, type = "numeric")
parser <- add_option(parser, '--beta_del', default = NA, type = "numeric")
parser <- add_option(parser, '--alpha_del', default = 1, type = "numeric")
args <- parse_args(parser)

subcl_p <- c(0.3, 0.4, 0.3)
gamma_alpha_amp <- args$alpha_amp
gamma_beta_amp <- args$beta_amp
gamma_alpha_del <- args$alpha_del
gamma_beta_del <- args$beta_del
print(paste("Amp alpha =", gamma_alpha_amp))
print(paste("Amp beta =", gamma_beta_amp))
print(paste("Del alpha =", gamma_alpha_del))
print(paste("Del beta =", gamma_beta_del))
cnv_factor_mat <- cnv_signal_strength_mat(count_mat = gene_count_mat, 
                                          dna_dosage_mat = subclone.gene.cnv,
                                          subclone_prop = subcl_p,
                                          shape.amp = gamma_alpha_amp,
                                          scale.amp = gamma_beta_amp,
                                          shape.del = gamma_alpha_del,
                                          scale.del = gamma_beta_del)

##### 4.4
# Call functions with output paths
gene_output <- paste0(filename, ".gene.countmatrix", "_", as.character(gamma_alpha_amp), 
                      "_", as.character(gamma_beta_amp), "_", as.character(gamma_alpha_del),
                      "_", as.character(gamma_beta_del))
#intergene_output <- paste0(filename, ".intergene.countmatrix")
# Test our CNV region & magnitude modified scDesign3_SyntheticCount()!
tryCatch({
  gene_res <- scDesign3_SyntheticCount(gene_count_mat, gene_output, cnv_factor_mat = cnv_factor_mat, cell_label = cellType)
  # load cell label file (cluster from previous scDesign3 function)
  cell_label <- read.table(file = paste0(gene_output, ".scDesign3Simulated.CellTypeLabel.txt"), sep="\t",header = FALSE)
  #intergene_res <- scDesign3_SyntheticCount(intergene_count_mat, intergene_output, cell_label[,1])
}, 
  error = function(e) {
  cat("Error in scDesign3_SyntheticCount:", conditionMessage(e), "\\n")
})

##### 4.5
# Plot a heatmap of synthetic count matrix with CNV effects
synthetic_count <- read.table(paste0(gene_output, ".scDesign3Simulated.txt"), sep="\t",header = FALSE)
rownames(synthetic_count) <- synthetic_count[, 1]
synthetic_count <- synthetic_count[, 2:ncol(synthetic_count)]
library(matrixStats)
library(pheatmap)
# CP10k + log1p
lib <- colSums(synthetic_count); keep <- lib > 0
X <- log1p( t( t(synthetic_count[, keep, drop=FALSE]) / lib[keep] ) * 1e4 )

# keep genes with signal (≥ 5% nonzeros)
X <- X[rowMeans(synthetic_count[, keep, drop=FALSE] > 0) >= 0.05, , drop=FALSE]

# row z-score, then clip to reduce effect of outliers
Z <- t(scale(t(X)))
#Z <- pmin(pmax(Z, -3), 3)

options(bitmapType = "cairo")
pheatmap(Z, scale = "none",
         color = colorRampPalette(c("navy","white","firebrick3"))(100),
         breaks = seq(-3, 3, length.out = 101),
         show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE,
         filename = paste0(gene_output, ".count_heatmap.png"), width = 8, height = 6)


