

################################################################

# this is the code analyzing the RNA-seq data of in vivo AT2-Tfam KO in mice
# main tools: `DESeq2`
# author: Qianjiang Hu
# date: 2026-03-17
################################################################



## ===========================================================================
## 0) Packages & options
## ===========================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(msigdbr)
  library(GSEABase)
  library(dplyr)
  library(readxl)
})

set.seed(1)

## ===========================================================================
## 1) Build the ranked vector
## ===========================================================================
{
  ranks_weighted <- as.data.frame(res) %>%
    na.omit() %>%
    mutate(
      Gene = rownames(.),
      pvalue_safe = ifelse(pvalue == 0, .Machine$double.xmin, pvalue),
      Rank = -log10(pvalue_safe) * sign(log2FoldChange)
    ) %>%
    arrange(desc(Rank)) %>%
    select(Gene, Rank) %>%
    na.omit()
  head(ranks_weighted)
  
  ## --- convert the ID to gene symbol ---
  {
    listEnsembl(version = 108)  # check the release 108
    # create Mart object
    ensembl = useEnsembl(biomart = 'genes', dataset='mmusculus_gene_ensembl', verbose = TRUE,
                         host="https://oct2022.archive.ensembl.org")
    ranks_weighted <- as.data.frame(ranks_weighted)
    genes <- rownames(ranks_weighted)
    annot <- getBM(attributes = c('ensembl_gene_id','mgi_symbol'),
      filters = 'ensembl_gene_id',
      values = genes,
      verbose = TRUE,
      mart = ensembl)
    names(annot)
    ranks_weighted$ID <- rownames(ranks_weighted)
    ranks_weighted <- merge(x = ranks_weighted,
                            y =  annot,
                            by.y = 'ensembl_gene_id',
                            all.x = T,
                            by.x = 'ID')
    
    colnames(ranks_weighted)
    ranks_weighted <- ranks_weighted[, c(1, ncol(ranks_weighted), 2:(ncol(ranks_weighted)-1))]
    head(ranks_weighted)
  }
  ###  --- convert ranked table to vector ---
  {
    prepare_gsea_rank <- function(ranked_df) {
      gene_list <- setNames(ranked_df$Rank, ranked_df$mgi_symbol)
      gene_list <- sort(gene_list, decreasing = TRUE)  # Required for fgsea
      return(gene_list)
    }
    ranked_list <- prepare_gsea_rank(ranks_weighted)
    ranked_list <- ranked_list[!duplicated(names(ranked_list))]
  }
}


## ===========================================================================
## 2) Gene Set Collections
## ===========================================================================
#### input Senescence list

# Read specific sheet and define column headers
file_path <- './senescence_list.xlsx'

##### SenMayo_senescence
SenMayodata <- read_excel(file_path, sheet = 'SenMayo', col_names = TRUE)
SenMayo_sen <- SenMayodata$symbol

##### CellAge_senescence
CellAgedata <- read_excel(file_path, sheet = 'CellAge Senescence Genes', col_names = TRUE)
CellAge_sen <- CellAgedata$Symbol

#### convert human gene as mouse gene

# Read the homologous gene mapping table
homologues <- read.delim("./matched_homologues.txt", header = TRUE, stringsAsFactors = FALSE)
# Extract relevant columns
homologues_filtered <- homologues %>% select(human, Symbol)

# Function to convert human gene symbols to mouse gene symbols
convert_human_to_mouse <- function(gene_list) {
  converted_genes <- homologues_filtered %>%
    filter(human %in% gene_list) %>%
    pull(Symbol)
  return(converted_genes)
}
# Convert human gene list to mouse gene list
SenMayo_sen <- convert_human_to_mouse(SenMayo_sen)
CellAge_sen <- convert_human_to_mouse(CellAge_sen)

#####  Prepare the Custom Gene Sets
custom_sen_sets <- rbind(
  data.frame(GeneSet = "SenMayo_sen", gene = SenMayo_sen),
  data.frame(GeneSet = "CellAge_sen", gene = CellAge_sen))
  
## ===========================================================================
## 2) RUN clusterProfiler::GSEA()
## ===========================================================================

# Verify structure
ranked_list
head(custom_sen_sets)

## --- run GSEA on Senescence gene sets ---
gsea_sen <- GSEA(
  geneList     = ranked_list, 
  TERM2GENE    = custom_sen_sets, 
  pvalueCutoff = 1,
  verbose      = TRUE,
  eps          = 0,  
  minGSSize    = 15,
  maxGSSize    = 500,
  seed         = TRUE,
  by           = "fgsea" 
)

## --- save the table ---
out_dir <- file.path(rs_dir, "gsea")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(gsea_sen@result, file.path(out_dir, "GSEA_Senescence.csv"), row.names = FALSE)
  

## ===========================================================================
## ridgeplot 
##    
## ===========================================================================

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# GSEA Senescence 
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

library(viridis)  # Load the viridis package
p <- ridgeplot(gsea_sen, showCategory = 2) +
  scale_fill_gradientn(colors = viridis(10, option = "plasma"), limits = c(0, 0.05))

# save
out_dir <- file.path(pl_dir, "gsea_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
safe_name <- gsub("[^A-Za-z0-9._-]+", "_", "senescence")
png_file  <- file.path(out_dir, paste0(safe_name, "_ridgeplot.png"))
ggsave(png_file, plot = p, width = 5, height = 2.5, dpi = 300)
    
  
  
  


