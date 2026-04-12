# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#       CONGAS with reference     # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
library(optparse)
parser = OptionParser(description = 'Execute CONGAS with reference')
parser = add_option(parser, '--rootDir', default = "/your/local/path/to/repo/", type = "character")
parser = add_option(parser, '--sampleID', default = "toy_data", type = "character")
parser = add_option(parser, '--refPath', default = "", type = "character")
args <- parse_args(parser)

library(Rcongas)
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(biomaRt)
root_Dir <- args$rootDir
sample_id <-args$sampleID
ref_Path <- args$refPath
input_Dir <- paste0(root_Dir, "data/", sample_id, "/filtered_feature_bc_matrix")
Out_Dir <- paste0(root_Dir, "Result/", sample_id, "/CONGAS/")

dir.create(Out_Dir, showWarnings = FALSE, recursive = TRUE)
setwd(Out_Dir)
message("Method: CONGAS ...")
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

# segmentation
data(human_chr_locations)
segments <- data.frame(human_chr_locations[,c("chr", "start", "end")],
                       copies = 2)
colnames(segments) <- c("chr", "from", "to", "copies")
segments$from[segments$from == 1] <- 0
segments$chr <- paste0("chr", segments$chr)
# Remove these chromosomes
segments = tibble(segments) %>% dplyr::filter(chr != 'chrX', chr != 'chrY')
segments$from <- as.integer(segments$from)
segments$copies <- as.integer(segments$copies)

# features
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
z <- getBM(c("hgnc_symbol", 'chromosome_name','start_position','end_position'),
           filters = "hgnc_symbol",
           rownames(input_expr),
           ensembl) %>%
  dplyr::rename(gene = hgnc_symbol,
                chr = chromosome_name,
                from = start_position,
                to = end_position) %>%
  filter(chr %in% c(seq(1:22), 'X', 'Y')) %>%
  mutate(chr = paste0('chr', chr))

z <- z[!duplicated(z$gene),]
gene_names <- rownames(input_expr)
features <- tibble(gene = gene_names) %>% left_join(z) 

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#            Execute              # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
rna <- Rcongas::create_congas_tibble(counts = input_expr, 
                                     modality = 'RNA', 
                                     save_dir = Out_Dir, 
                                     features = features)
# Compute normalization factors
norm_rna <- Rcongas:::auto_normalisation_factor(rna) %>%
  mutate(modality = 'RNA')

rna <- rna %>% filter_known_genes(what='r') # Ribosomal
all_genes <- rna$gene %>% unique
mito <- all_genes %>% str_starts(pattern = 'MT-') # Mitochondrial DNA
all_genes <- setdiff(all_genes, all_genes[mito])
rna <- rna %>% dplyr::filter(gene %in% all_genes)
# force all expression as integer
rna$value <- as.integer(rna$value)

# initialization
x <- init(
  rna = rna,
  atac = NULL,
  segmentation = segments, 
  rna_normalisation_factors = norm_rna,
  atac_normalisation_factors = NULL,
  rna_likelihood = 'NB',
  description = '')

filt <- segments_selector_congas(x)

k = 4:8
binom_limits = c(40,1000)
model = "BIC"
lr = 0.01
temperature = 20
steps = 5000
lambda = 0.5
# Estimate hyperparameters
hyperparams_filt <- auto_config_run(filt, k, 
                                    prior_cn = c(0.2, 0.6, 0.1, 0.05, 0.05),
                                    init_importance = 0.6, 
                                    CUDA = FALSE, 
                                    normal_cells = TRUE)  # set to FASLE if no normal cells in the sample

hyperparams_filt$binom_prior_limits <- binom_limits

fit_filt <- Rcongas:::fit_congas(filt,
                                 K = k,
                                 lambdas = 1, # set to 1 if using RNA data only
                                 learning_rate = lr, 
                                 steps = steps,
                                 model_parameters = hyperparams_filt, 
                                 model_selection = model,
                                 latent_variables = "G",
                                 CUDA = FALSE,
                                 temperature = temperature, 
                                 same_mixing = TRUE, 
                                 threshold = 0.001)
saveRDS(fit_filt, paste0(Out_Dir, "congas_wRef_res.rds"))
