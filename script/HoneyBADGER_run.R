# %%%%%%%%%%%%%%%%%%%%% #
#     HoneyBADGER       # 
# %%%%%%%%%%%%%%%%%%%%% #
library(optparse)
parser = OptionParser(description = 'Execute HoneyBADGER ')
parser = add_option(parser, '--rootDir', default = "/your/local/path/to/repo/", type = "character")
parser = add_option(parser, '--sampleID', default = "toy_data", type = "character")
parser = add_option(parser, '--refPath', default = "", type = "character")
parser = add_option(parser, '--chrPrefix', default = TRUE, type = "logical")
args <- parse_args(parser)

library(HoneyBADGER)
library(Seurat)
library(biomaRt)
library(parallel)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
root_Dir <- args$rootDir
sample_id <-args$sampleID
ref_Path <- args$refPath
chr_prefix <- args$chrPrefix
input_Dir <- paste0(root_Dir, "data/", sample_id, "/filtered_feature_bc_matrix/")
Out_Dir <- paste0(root_Dir, "Result/", sample_id, "/HoneyBADGER/")

dir.create(Out_Dir, showWarnings = FALSE, recursive = TRUE)
setwd(Out_Dir)
message("Method: HoneyBADGER ...")
message("Process ", sample_id, " ...")
message("Output Directory: ", Out_Dir)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#           Input data            # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
input_expr <- as.matrix(Read10X(data.dir = input_Dir))
ref_expr <- as.matrix(read.csv(ref_Path, row.names = 1))

overlapped.genes <- intersect(rownames(input_expr), rownames(ref_expr))
input_expr <- input_expr[overlapped.genes,]
ref_expr <- ref_expr[overlapped.genes,]

# scale for library size differences 
input_expr <- scale(input_expr)
ref_expr <- rowMeans(ref_expr)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#        Execute (expression mode)            # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                    dataset = 'hsapiens_gene_ensembl', 
                    mirror = "useast") ## current version 20210603- default version hg38

hb <- new('HoneyBADGER', name = sample_id)
hb$setGexpMats(input_expr, ref_expr, mart.obj, filter = FALSE, scale = FALSE, verbose = TRUE)
hb$setMvFit(verbose = TRUE) ## model variance
hb$setGexpDev(verbose = TRUE) ## model necessary expression deviation to identify CNVs
hb$calcGexpCnvBoundaries(init = TRUE, verbose = FALSE, min.num.genes = 3) ## HMM

# posterior probability of each CNV in each cell
hb$retestIdentifiedCnvs(retestBoundGenes = TRUE, retestBoundSnps = FALSE, verbose = FALSE)
# final results
results <- hb$summarizeResults(geneBased = TRUE, alleleBased = FALSE)

save(hb, file = paste0(Out_Dir, "HMM_expr_out.RData"))
save(results, file = paste0(Out_Dir, "CNV_expr_result.RData"))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#      Allele data preparation    # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# list of snps for allele count matrix
snps.out <<- get(load(paste0(root_Dir, "/resources/HoneyBADGER_ref/ExAC_hg38_chr1_22_snps_overlapped.RData")))
snps.out <- snps.out[unique(names(snps.out)),]
snps <- data.frame(chr = as.character(snps.out@seqnames),
                   pos = as.integer(as.character(snps.out@ranges)),
                   ref = as.character(snps.out$REF),
                   alt = as.character(unlist(snps.out$ALT)))
if(!chr_prefix){
  # to check chr prefix: samtools idxstats your.bam | cut -f1
  message("remove chr prefix ...")
  snps$chr <- gsub("chr", "", snps$chr) # do not remove chr if bam file contains chr
}
## filter out snps that are actually indels and warn
kv <- (sapply(snps[[3]],nchar) == 1 & sapply(snps[[4]], nchar) == 1)
snps <- snps[kv,]
tmp.snps.file = paste0(Out_Dir, "snps.txt")
write.table(snps, 
            tmp.snps.file, 
            quote = F, col.names = F, row.names = F)

# cell barcode for allele count matrix
cellBarcodes <- colnames(input_expr)
tmp.barcodes.file = paste0(Out_Dir, "barcodes.txt")
write.table(cellBarcodes,
            tmp.barcodes.file, 
            quote = F, col.names = F, row.names = F)

# apply scAlleleCount
scAlleleCountExec <- paste0(root_Dir, "script/util/HoneyBADGER/scAlleleCount.py")
tmp.output.prefix <- paste0(Out_Dir, "AlleleCount")
path <- paste0(input_Dir, "../")
files <- list.files(path = path)
# list of paths to bam files (sorted and indexed)
bamFiles <- files[grepl('.bam$', files)]
bamFiles <- paste0(path, bamFiles)

cmd <- paste0("python3.8 ", scAlleleCountExec, ' -v' ,
              ' --snps ', tmp.snps.file, ' --barcodes ', tmp.barcodes.file,
              ' --output-format mm --max-depth 9999999 --output-prefix ', tmp.output.prefix, ' ',
              ' --bamfile ', bamFiles)
cmd.ret <- system(cmd)

# save results for future uses since scAlleleCount takes long time to run
r <- t(Matrix::readMM(paste0(Out_Dir, 'AlleleCountrefmat.mtx')))
alt <- t(Matrix::readMM(paste0(Out_Dir, 'AlleleCountaltmat.mtx')))
rownames(r) <- rownames(alt) <- paste(snps$chr, snps$pos, sep = ":")
colnames(r) <- colnames(alt) <- cellBarcodes
r <- as.matrix(r)
alt <- as.matrix(alt)
cov.sc <- r+alt

save(r, file = paste0(Out_Dir, "ref_allele.RData"))  ## alternate allele
save(cov.sc, file = paste0(Out_Dir, "coverage.RData"))  ## total coverage

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#      Execute (allele mode)      # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hb$setAlleleMats(r.init = r, n.sc.init = cov.sc, het.deviance.threshold = 0.05, n.cores = 5)
hb$setGeneFactors(txdb)
hb$calcAlleleCnvBoundaries(init = TRUE, verbose = FALSE)
hb$retestIdentifiedCnvs(retestBoundGenes = FALSE, retestBoundSnps = TRUE, verbose = FALSE)

## final results
results <- hb$summarizeResults(geneBased = FALSE, alleleBased = TRUE)
save(hb, file = paste0(Out_Dir, "HMM_allele_out.RData"))
save(results, file = paste0(Out_Dir, "CNV_allele_result.RData"))
