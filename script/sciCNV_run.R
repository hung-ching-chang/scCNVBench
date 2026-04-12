# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#             sciCNV              # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
library(optparse)
parser = OptionParser(description = 'Execute sciCNV')
parser = add_option(parser, '--rootDir', default = "/your/local/path/to/repo/", type = "character")
parser = add_option(parser, '--sampleID', default = "toy_data", type = "character")
parser = add_option(parser, '--refPath', default = "", type = "character")
args <- parse_args(parser)

root_Dir <- args$rootDir
sample_id <-args$sampleID
ref_Path <- args$refPath
sciCNV.code <- paste0(root_Dir, "script/util/sciCNV-Analysis/")
input_Dir <- paste0(root_Dir, "Real_data/", sample_id, "/filtered_feature_bc_matrix")
Out_Dir <- paste0(root_Dir, "Result/", sample_id, "/sciCNV/")

dir.create(Out_Dir, showWarnings = FALSE, recursive = TRUE)
setwd(Out_Dir)
message("Method: sciCNV ...")
message("Process ", sample_id, " ...")
message("Output Directory: ", Out_Dir)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#    load package and function    # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
library(base)
library(robustbase)
library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(qlcMatrix)
library(svd)
library(biomaRt)

source(file.path(sciCNV.code, "Mito_umi_gn.R"))
source(file.path(sciCNV.code, "RTAM_normalization.R"))
source(file.path(sciCNV.code, "Scaling_CNV.R"))
source(file.path(sciCNV.code, "CNV_score.R"))
source(file.path(sciCNV.code, "sciCNV.R"))
source(file.path(sciCNV.code, "Sketch_AveCNV.R"))
source(file.path(sciCNV.code, "CNV_htmp_glist.R"))
source(file.path(sciCNV.code, "CNV_htmp_gloc.R"))
source(file.path(sciCNV.code, "Opt_MeanSD_RTAM1.R"))
source(file.path(sciCNV.code, "Opt_MeanSD_RTAM2.R"))
source(file.path(sciCNV.code, "heatmap_break_glist.R"))
source(file.path(sciCNV.code, "heatmap_break_gloc.R"))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#           Input data            # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
## Reading raw data with a list of genes on the first column
input_expr <- as.matrix(Read10X(data.dir = input_Dir, strip.suffix = T))
ref_expr <- as.matrix(read.csv(ref_Path, row.names = 1))
merge_expr <- cbind(input_expr[intersect(rownames(input_expr),rownames(ref_expr)),],
                    ref_expr[intersect(rownames(input_expr),rownames(ref_expr)),])

No.test <- ncol(input_expr)   # Number of test cells
No.control <- ncol(ref_expr)   # Number of control cells
test.names <- colnames(input_expr)
control.names <- colnames(ref_expr)

nUMI <- t(as.numeric(colSums(merge_expr)))
colnames(nUMI) <- colnames(merge_expr)
nUMI <- t(as.matrix(nUMI))

#----------------------------------------------------
##  Sorting of cells by UMI (from largest to smallest)
#----------------------------------------------------
if(No.control > 0){
  raw.data <- merge_expr[,c(order(as.data.frame(nUMI)[1:No.test,1], decreasing = TRUE),
                          order(as.data.frame(nUMI)[(No.test+1):ncol(merge_expr),1], decreasing = TRUE) + No.test), 
                         drop=FALSE]
} else {
  raw.data <- merge_expr[, order(as.data.frame(nUMI)[, 1], decreasing=TRUE), drop=FALSE]
}
rownames(raw.data) <- rownames(merge_expr)

#----------------------------
##  RTAM1/2 Normalization
#----------------------------
norm.data <- RTAM_normalization(mat = raw.data,            
                                method = "RTAM2",      
                                Min_nGn  = 250,       
                                Optimizing = FALSE)
print("RTAM normalization OK!")
No.test <- sum(colnames(norm.data) %in% test.names)

# rownames(norm.data) <- rownames(raw.data) 
# colnames(norm.data) <- colnames(raw.data)

save(norm.data, file = paste0(Out_Dir, "NormalziedData.RData"))

#-----------------------------------------------------------------------------------------------------------
##  Running the sciCNV function to derive CNV profiles for single cells from RTAM-normalized scRNA-seq data
#-----------------------------------------------------------------------------------------------------------
tst.index  <- seq(1, No.test , 1)                      # No of test cells
ctrl.index <- seq(No.test+1, ncol(norm.data), 1)      # No of control cells

## generating infered-CNV data for (test and/or control) cells 
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                    dataset = 'hsapiens_gene_ensembl', 
                    mirror = "useast")

gene.info <- getBM(values=rownames(norm.data),
                   attributes = c("hgnc_symbol", "chromosome_name","start_position","end_position"),
                   filters = "hgnc_symbol",
                   mart = mart.obj)

gene.info <- gene.info[order(gene.info$chromosome_name, gene.info$start_position),]
colnames(gene.info) <- c("Gene.name", "Chromosome", "Start", "end")
gene.info <- gene.info[gene.info$Chromosome %in% c(as.character(1:22), c("X", "Y")), ]
gene.info <- gene.info[!duplicated(gene.info$Gene.name), ]

Gen.loc <- gene.info
rownames(Gen.loc) <- Gen.loc$Gene.name

norm.data <- norm.data[intersect(rownames(Gen.loc), rownames(norm.data)),]
Gen.loc.my <- Gen.loc[intersect(rownames(Gen.loc), rownames(norm.data)), ]
avg.ctl <- row_means(norm.data[, ctrl.index])

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#            Execute              #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
CNV.data <- sciCNV(norm.mat = norm.data,
                   gen.Loc = Gen.loc.my[, c(1,2)],  
                   ave.ctrl = avg.ctl,
                   No.test = No.test, 
                   sharpness  = 1, 
                   baseline_adj  = FALSE,  
                   baseline = 0)
print("CNV data OK!")
save(CNV.data, file = paste0(Out_Dir, "CNV.data.RData"))
save(Gen.loc.my, file = paste0(Out_Dir, "GenLoc.RData"))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#         Post-processing         #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#------------------------------------------------------------------
##  Scaling the sciCNV curves and setting a noise filter threshold
#------------------------------------------------------------------
## Scaling CNV-curves to adjust one copy number gain/losses to height +1/-1 if applicable
CNV.data <- CNV.data[intersect(Gen.loc.my$Gene.name, rownames(CNV.data)), ]
Gen.loc.my <- Gen.loc[intersect(Gen.loc.my$Gene.name, rownames(CNV.data)), ]
CNV.data.scaled <- Scaling_CNV(V7Alt = CNV.data, 
                               n.TestCells = No.test, 
                               Gen.Loc = Gen.loc.my,
                               scaling.factor = 1)

## Defining M_NF as Noise-Free Matrix of test cells (and control cells)
M_NF <- CNV.data.scaled

#######  Noise Filteration after scaling
noise.thr = 0.4   # Noise threshold
for(w in 1:ncol(M_NF) ){
  for(j in 1:nrow(M_NF)){
    if( (M_NF[j,w] > -noise.thr) && (M_NF[j,w] < noise.thr)  ){
      M_NF[j,w] <- 0
    }
  }
}
M_NF <- as.matrix(M_NF)

#### Taking Square Rate
for(w in 1:ncol(M_NF)){
  for(j in 1:nrow(M_NF)){
    if (M_NF[j,w] > 0){
      M_NF[j,w]<- sqrt( as.numeric(M_NF[j,w]))
    } else if (M_NF[j,w] < 0){
      M_NF[j,w]<- -sqrt( -as.numeric(M_NF[j,w]))
    }
  }
}

#------------------------------------------------------------------
##  Calculating tumor CNV scores to segregate normal vs tumor cells
#------------------------------------------------------------------
TotScore <- CNV_score(M_nf = M_NF)
# CNVscore <- data.frame(group = rep(c("Test", "Control"), c(No.test, ncol(M_NF) - 1 - No.test)),
#                        score = TotScore[1,])
save(TotScore, file = paste0(Out_Dir, "TotScore.RData"))


M_NF <- M_NF[,1:No.test] # remove control cells
save(M_NF, file = paste0(Out_Dir, "M_NF_rm_noise.RData"))
