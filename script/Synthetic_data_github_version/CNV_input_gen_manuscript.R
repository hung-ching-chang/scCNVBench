# This file is used to approximate the range of C_ij, the designed parameter 
# representing DNA to RNA CNV translation strength. This rough distribution of log(1+C_ij)
# will serve as user instruction of how to choose the gamma_alpha_amp, gamma_beta_amp,
# gamma_alpha_del, gamma_beta_del for strong signal CNV and weak signal CNV simulated.

# This code will learn the rough distribution from SNU638, user may change to other
# samples they want to create for synthetic data.

library(GenomicRanges)
library(dplyr)
###### Generate a subclone bin level and gene level table
subclone_bin_cnv <- readRDS("path/to/SNU638subclone_bin_cnv.rds")
gene.info <- readRDS("path/to/gene_info.rds")
bin.info <- subclone_bin_cnv$bin.info
subclone.bin.cnv <- subclone_bin_cnv$subclone.cnv

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
head(gene.meta, 10)

# bin-level to gene-level
subclone.gene.cnv <- subclone.bin.cnv[gene.meta$bin_pos,]
rownames(subclone.gene.cnv) <- gene.meta$gene
subclone.gene.cnv <- cbind(subclone.gene.cnv, gene.meta)
head(subclone.gene.cnv)

# Enlarge the deletion segment on chr14
# Make chr14:67040000 -chr14:107043718 all having deletion CNV
# subclone.gene.cnv[(subclone.gene.cnv$chr=="chr14" & 
#                      as.numeric(subclone.gene.cnv$start) >= 67040000 &
#                      as.numeric(subclone.gene.cnv$start) < 107043718), c("c1", "c2", "c3")] <- rep(0.5032941, 3)
# saveRDS(subclone.gene.cnv, file = "C:/University/2025_Summer/GSR/SNU638/subclone.gene.cnv_mod.rds")     

##### Modify this table so that segments will not lie on two chromosomes (from end of chrN to start of chrN+1)
##### We do this by adding a blank row after the end row of each chromosome
# Get the indices of the last row of each chromosome
end_rows <- tapply(1:nrow(subclone.gene.cnv), subclone.gene.cnv$chr, tail, 1)
levels_chr <- paste0("chr", c(1:22, "X", "Y", "M"))
end_rows <- end_rows[order(factor(names(end_rows), levels = levels_chr), na.last = TRUE)]

# Create a blank row template with same columns
blank_row <- subclone.gene.cnv[1, ]
blank_row[,] <- NA

# Insert blank rows
out <- subclone.gene.cnv
shift <- 0
for (idx in end_rows) {
  if (idx!=max(end_rows)){
    out <- rbind(
      out[1:(idx + shift), ],
      blank_row,
      out[(idx + shift + 1):nrow(out), ])
  } else{
    out <- rbind(
      out[1:(idx + shift), ],
      blank_row)
  }
  
  shift <- shift + 1
}

subclone.gene.cnv <- out

##### Find CNV segments
library(data.table)
# Exclude genes on chromosome X
# Find genes with all 3 subclones having CNV behavior
cnv_genes <- which(
  ((subclone.gene.cnv$c1 > 1.25 & subclone.gene.cnv$c2 > 1.25 & subclone.gene.cnv$c3 > 1.25) |
     (subclone.gene.cnv$c1 < 0.75 & subclone.gene.cnv$c2 < 0.75 & subclone.gene.cnv$c3 < 0.75)) &
    subclone.gene.cnv$chr != "chrX")
seg_start_end <- data.table(cnv_genes)[
  , .(start = min(cnv_genes), end = max(cnv_genes)), 
  by = rleid(cnv_genes - seq_along(cnv_genes))][, .(start, end)]
seg_start_end$segment <- 1:nrow(seg_start_end)
summary_table <- data.frame(
  start_pos = paste0(subclone.gene.cnv[seg_start_end$start, "chr"],":", 
                     subclone.gene.cnv[seg_start_end$start, "start"]),
  end_pos = paste0(subclone.gene.cnv[seg_start_end$end, "chr"], ":",
                   subclone.gene.cnv[seg_start_end$end, "end"]),
  n_gene = seg_start_end$end - seg_start_end$start + 1,
  mean_dosage = mapply(function(a,b) mean(unlist(subclone.gene.cnv[a:b, c("c1", "c2", "c3")])), 
                       seg_start_end$start, seg_start_end$end),
  sd_dosage = mapply(function(a,b) sd(unlist(subclone.gene.cnv[a:b, c("c1", "c2", "c3")])), 
                     seg_start_end$start, seg_start_end$end),
  min_dosage = mapply(function(a,b) min(unlist(subclone.gene.cnv[a:b, c("c1", "c2", "c3")])), 
                      seg_start_end$start, seg_start_end$end),
  max_dosage = mapply(function(a,b) max(unlist(subclone.gene.cnv[a:b, c("c1", "c2", "c3")])), 
                      seg_start_end$start, seg_start_end$end)
)
print(summary_table, digits=2)

# This is the Dg file we need!
allcell_seg_gene_dosage <- data.frame(gene = subclone.gene.cnv[cnv_genes, "gene"],
                                                   dosage = rowMeans(subclone.gene.cnv[cnv_genes, c("c1", "c2", "c3")]),
                                                   segment = cumsum(c(1, diff(cnv_genes) > 1)))
print(dim(allcell_seg_gene_dosage))
print(head(allcell_seg_gene_dosage))


##### Crafting for Eg file
library(Matrix)
library(Seurat)
# r_tumor
scRNA_expr <- Read10X("path/to/SNU638_filtered_feature_bc_matrix")
scRNA_expr <- CreateSeuratObject(counts = Matrix(as.matrix(scRNA_expr), sparse = TRUE), project = "simulation")
# r_normal
scRNA_normal <- read.csv("path/to/GSM4850590_Stomach_Counts.csv")
rna_normal_ref <- CreateSeuratObject(counts = Matrix(as.matrix(scRNA_normal), sparse = TRUE), project = "simulation")
# Convert both count matrix to cpm
scRNA_cpm <- edgeR::cpm(as.matrix(Seurat::GetAssayData(scRNA_expr, slot = "counts")))
rna_normal_ref_cpm <- edgeR::cpm(as.matrix(Seurat::GetAssayData(rna_normal_ref, slot = "counts")))
# Select only genes with averaged expr > 0
rna_normal_ref_cpm <- rna_normal_ref_cpm[rowSums(rna_normal_ref_cpm) > 0, ]

# Find overlapping genes
common_genes <- intersect(rownames(scRNA_cpm), rownames(rna_normal_ref_cpm))
scRNA_cpm <- scRNA_cpm[common_genes, , drop = FALSE]
rna_normal_ref_cpm <- rna_normal_ref_cpm[common_genes, , drop = FALSE]
# Get Eg by r_tumor/r_normal
rna_ref_gene_level <- rowMeans(rna_normal_ref_cpm, na.rm = TRUE)
rna_CNV_ratio <- sweep(scRNA_cpm, 1, rna_ref_gene_level, "/") 

library(patchwork)
library(ggplot2)
library(fitdistrplus)
library(goftest)   

# Use all cells from scRNA_cpm
subclone_gene_rna_ratio <- rowMeans(rna_CNV_ratio[, colnames(scRNA_cpm)], na.rm = TRUE)
subclone_gene_rna_ratio <- data.frame(gene = names(subclone_gene_rna_ratio), 
                                      cnv_ratio = unname(subclone_gene_rna_ratio), 
                                      stringsAsFactors = FALSE)
# Inner join with DNA dosage segment information
subclone_seg_dosage <- allcell_seg_gene_dosage
Dg_Eg_Cg <- inner_join(subclone_seg_dosage, subclone_gene_rna_ratio, by = "gene")
Dg_Eg_Cg$Cg <- Dg_Eg_Cg$cnv_ratio / Dg_Eg_Cg$dosage

Cg_hist_plots <- list()
Cg_log1p_hist_plots <- list()
Cg_gamma_diag_plots <- list()
for (seg in unique(Dg_Eg_Cg$segment)) {
  Cg_seg <- subset(Dg_Eg_Cg, segment == seg)
  if (dim(Cg_seg)[1]<5) next
  Cg_hist_plots[[length(Cg_hist_plots) + 1]] <-
    ggplot(Cg_seg, aes(x = Cg)) +
    geom_histogram(bins = 30) +
    stat_bin(
      bins = 30, geom = "text",
      aes(y = after_stat(count),
          label = after_stat(ifelse(count > 0, count, NA))),
      vjust = -0.5, size = 3) +
    labs(x = "Cg", y = "Frenquency", title = paste("segment", seg, "n_gene=", dim(Cg_seg)[1]))
  x <- Cg_seg$Cg
  x_trans <- log1p(x)
  x_trans <- x_trans[x_trans>0 & is.finite(x_trans)]
  if (length(x_trans)==0) next
  fit <- fitdistrplus::fitdist(x_trans, "gamma")
  shape <- unname(fit$estimate["shape"])
  scale <- 1 / unname(fit$estimate["rate"])
  ad_p <- goftest::ad.test(x_trans, null = pgamma, shape = shape, scale = scale)$p.value
  print(paste("segment", seg, "has GOF p-val", ad_p))
  if (ad_p<0.05) {
    par(mfrow=c(2,2))
    plot(fit)
    title(main = paste("segment", seg))
  }
  
  # Get chr position label from summary_table
  seg_row <- summary_table[seg, ]
  chr_label <- paste0(seg_row$start_pos, " : ", seg_row$end_pos)
  
  # Compute gamma density overlay
  x_seq <- seq(min(x_trans), max(x_trans), length.out = 200)
  gamma_df <- data.frame(
    x = x_seq,
    y = dgamma(x_seq, shape = shape, scale = scale)
  )
  
  # Scale density to histogram counts
  bin_width <- diff(range(x_trans)) / 30
  gamma_df$y_scaled <- gamma_df$y * length(x_trans) * bin_width
  
  # Annotation label for top right
  annot_label <- paste0("alpha = ", round(shape, 3),
                        "\nbeta = ", round(scale, 3),
                        "\nnum of genes = ", length(x_trans))
  
  Cg_log1p_hist_plots[[length(Cg_log1p_hist_plots) + 1]] <-
    ggplot(data.frame(x_trans = x_trans), aes(x = x_trans)) +
    geom_histogram(bins = 30, fill = "grey40", color = "white") +
    stat_bin(
      bins = 30, geom = "text",
      aes(y = after_stat(count),
          label = after_stat(ifelse(count > 0, count, NA))),
      vjust = -0.5, size = 3) +
    geom_line(data = gamma_df, aes(x = x, y = y_scaled), color = "red", linewidth = 1) +
    annotate("text", x = Inf, y = Inf,
             label = annot_label,
             hjust = 1.1, vjust = 1.3,
             size = 3.5) +
    labs(x = "log(1+Cg)",
         y = "Frequency",
         title = paste0(chr_label)) +
    theme_bw()
}
wrap_plots(Cg_hist_plots, ncol = 3) |> print()
wrap_plots(Cg_log1p_hist_plots, ncol = 3) |> print()