# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#             CaSpER              # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
library(optparse)
parser = OptionParser(description = 'Execute CaSpeER')
parser = add_option(parser, '--rootDir', default = "/your/local/path/to/repo/", type = "character")
parser = add_option(parser, '--sampleID', default = "toy_data", type = "character")
parser = add_option(parser, '--refPath', default = "", type = "character")
args <- parse_args(parser)

root_Dir <- args$rootDir
sample_id <-args$sampleID
ref_Path <- args$refPath
input_Dir <- paste0(root_Dir, "Real_data/", sample_id, "/filtered_feature_bc_matrix")
Out_Dir <- paste0(root_Dir, "Result/", sample_id, "/CaSpER/")

dir.create(Out_Dir, showWarnings = FALSE, recursive = TRUE)
setwd(Out_Dir)
message("Method: CaSpER ...")
message("Process ", sample_id, " ...")
message("Output Directory: ", Out_Dir)

library(Seurat)
library(CaSpER)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#     cytoband preprocessing      # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
data(hg38_cytoband) # cytoband
cytoband_hg38 <- data.frame(V1 = gsub("chr", "", cytoband_hg38[,1]), 
                            V2 = cytoband_hg38[,2], 
                            V3 = cytoband_hg38[,3], 
                            V4 = substring(cytoband_hg38$V4, 1, 1), 
                            stringsAsFactors = F)
start <- do.call(rbind, lapply(split(cytoband_hg38$V2, paste0(cytoband_hg38$V1, cytoband_hg38$V4)), min))
end <- do.call(rbind, lapply(split(cytoband_hg38$V3, paste0(cytoband_hg38$V1, cytoband_hg38$V4)), max))
cytoband_hg38 <- data.frame(V1 = gsub("p", "", gsub("q", "", rownames(start))), 
                            V2 = start, 
                            V3 = end, 
                            V4 = rownames(start), 
                            stringsAsFactors = F)
cytoband_hg38 <- cytoband_hg38 [as.vector(unlist(sapply(c(1:22, "X"), function(x) which(cytoband_hg38$V1 %in% x)))), ]
cytoband_hg38$V4[grep("q", cytoband_hg38$V4)] <- "q"
cytoband_hg38$V4[grep("p", cytoband_hg38$V4)] <- "p"
rownames(cytoband_hg38) <- NULL
head(cytoband_hg38)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#           Input data            # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
input_expr <- Read10X(data.dir = input_Dir, strip.suffix = T) 
ref_expr <- as.matrix(read.csv(ref_Path, row.names = 1))
colnames(ref_expr) <- paste0("Ref_", colnames(ref_expr))

overlapped.genes <- intersect(rownames(input_expr), rownames(ref_expr))
merge_expr <- as.matrix(cbind(input_expr[overlapped.genes,], ref_expr[overlapped.genes,]))

annotation <- read.csv(paste0(root_Dir, "resources/CaSpER_ref/gene_annotation.csv"),
                       row.names = 1)
annotation <- annotation[annotation$Gene %in% rownames(merge_expr),]
merge_expr <- merge_expr[annotation$Gene,] # filter genes
# read baf output
loh <- readBAFExtractOutput(path = Out_Dir, sequencing.type = "single-cell", suffix = "baf")
names(loh) <- sample_id

#Create sample - barcode mapping for LOH annotation
loh.name.mapping <- data.frame(loh.name = sample_id, sample.name = colnames(merge_expr))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#            Execute              #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
object <- CreateCasperObject(raw.data = merge_expr, annotation = annotation, 
                             control.sample.ids = colnames(ref_expr),  
                             loh.name.mapping = loh.name.mapping, 
                             cytoband = cytoband_hg38, cnv.scale = 3, 
                             loh.scale = 3, method = "iterative", 
                             loh = loh, matrix.type = "raw",
                             sequencing.type = "single-cell", 
                             expr.cutoff = 4.5)

final.objects <- runCaSpER(object, removeCentromere = T, method = "iterative")

save(final.objects, file = paste0(Out_Dir, "CaSpER_final_out.RData"))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#   Gene based CNV Summarization  #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
gamma <- 6
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary(final.objects)
save(segment.summary, file = paste0(Out_Dir, "segment_summary.RData"))
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loh <- segment.summary$all.summary.loh
loss.final <- loss[loss$count>gamma, ]
gain.final <- gain[gain$count>gamma, ]
loh.final <- loh[loh$count>gamma, ]

all.summary <- rbind(loss.final, gain.final)
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), 
                IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation, 
                                   keep.extra.columns = TRUE, 
                                   seqnames.field = "Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "GeneSymbol")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation[,2])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
rna.matrix <- gene.matrix(seg = all.summary, all.genes = all.genes, 
                          all.samples = all.samples, genes.ann = genes.ann)
save(rna.matrix, file = paste0(Out_Dir, "rna_matrix_raw.RData"))
rna.matrix  <- rna.matrix[,colnames(input_expr)]# remove control cells
save(rna.matrix, file = paste0(Out_Dir, "rna_matrix.RData"))

