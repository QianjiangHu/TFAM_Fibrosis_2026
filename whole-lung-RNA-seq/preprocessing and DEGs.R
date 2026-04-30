

################################################################

# this is the code analyzing the RNA-seq data of in vivo AT2-Tfam KO in mice
# main tools: `DESeq2`
# author: Qianjiang Hu
# date: 2026-03-17
################################################################


# install packages
list_of_packages <- c('DESeq2',
                      'ggplot2',
                      'ggfortify',
                      'biomaRt',
                      'org.Mm.eg.db',
                      'tidyverse',
                      'ComplexHeatmap',
                      'cowplot',
                      'circlize',
                      'RColorBrewer',
                      'clusterProfiler',
                      'enrichplot',
                      'msigdbr',
                      'DOSE',
                      'dplyr',
                      'EnhancedVolcano',
                      'pathview',
                      'fgsea',
                      'ggnewscale'
)

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,'Package'])]
if(length(new_packages) >0) BiocManager::install(new_packages, dependencies=TRUE)

# load packages
sapply(list_of_packages, library, character.only = T)


##### input data matrix

{
  ##### Set your working directory
  {
    # set and difine work folder path
    here::here()
    cm_dir <- here::here('Counts/featureCounts')
    rs_dir <- here::here('results')
    pl_dir <- here::here('plots')
    sp_dir <- here::here('scripts')
    me_dir <- here::here('metadata')
    ls_dir <- here::here('list')
    fd_ls  <- list(cm_dir,rs_dir,pl_dir,sp_dir,me_dir,ls_dir)
    
    # check the folder existance and create it if not exsiting
    for (i in fd_ls){
      if (!dir.exists(i)){
        dir.create(i)
        print(paste(basename(i), "---dir created"))
      }else{
        print(paste(basename(i), "---dir exists"))
      }
    }
    
  }

  
  ##### read in count matrix and merge samples to one matrix
  {
    # list the sample ID/file ID
    files <- dir(cm_dir, pattern = "*.txt$")
    file_ID <- unlist(strsplit(files, split = '*.featureCounts.txt'))
    coldata <- read.delim(paste(me_dir,'/', 'metadata.txt', sep = ''), header = TRUE, row.names=NULL, comment.char="#")
    head(coldata)
    
    # read in counts and combined together
    library(dplyr)
    
    data <- read.delim(paste(cm_dir,'/', file_ID[1], '.featureCounts.txt', sep = ''), header = TRUE, row.names=NULL, comment.char="#") %>% dplyr::select(c(1,6,7))
    head(data)
    colnames(data)[length(colnames(data))] <- file_ID[1]
    
    for (i in 2:length(file_ID)) {
      print(file_ID[i])
      print('...read dataframe')
      idata <- read.delim(paste(cm_dir,'/', file_ID[i], '.featureCounts.txt', sep = ''), header = TRUE, row.names=NULL, comment.char="#") %>% dplyr::select(c(1,7))  
      print('...rename colname')
      colnames(idata)[length(colnames(idata))]<- file_ID[i]
      print('...re-sort rows')
      print('...combine columns')
      data <- full_join(data, idata, by="Geneid")
      cat(file_ID[i]," count matrixs have merged!\n")
    }
    
    cat("All count matrixs have merged!\n")
    head(data)
    tail(data)
    dim(data)
    colnames(data)
  }
  
  
  ##### generate counts matrix and gene length table
  {
    cts <- data %>% dplyr::select(-2) %>% tibble::column_to_rownames(var = "Geneid")
    gene_length <- data %>% dplyr::select(1,2) %>% tibble::column_to_rownames(var = "Geneid")
    head(cts)
    head(gene_length)
    filter_cts <- cts[rowSums(cts)>0,]
    dim(cts)
    cat(nrow(cts)-nrow(filter_cts)," genes were not expressed in all samples!\n")
  }
  
  
  ##### Metadata for the samples
  {
    coldata$ID <- coldata$SampleID
    coldata <- coldata %>% tibble::column_to_rownames(var = "ID")
    all(rownames(coldata) %in% colnames(cts))
    all(rownames(coldata) == colnames(cts))
    cts <- cts[, rownames(coldata)]
    all(rownames(coldata) == colnames(cts))
    
    # set groups as factor
    coldata$SampleID <- factor(coldata$SampleID)
    coldata$Groups <- factor(coldata$Groups)
    coldata$Sample <- factor(coldata$Sample) 
  }
  
}


##### construct DESeqDataSet object `dds`

{
  ### put 'sex' as the covariate. 
  
  coldata$Sex    <- factor(coldata$Sex)
  coldata$Groups <- factor(coldata$Groups)
  coldata$Cohort <- factor(coldata$Cohort)
  
  # Check there is replication (no complete confounding)
  table(coldata$Groups, coldata$Cohort)
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData   = coldata,
                                design    = ~ Batch + Sex + Groups)
}



##### pre-filtering
{
  # check matrix
  head(assay(dds))
  nrow(dds)
  ncol(dds)
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep,]
  
  head(assay(dds))
  nrow(dds)
  
}


##### Make the DESeq object
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Groups","Homo_KO","WT"))

head(res) 
summary(res)



###################################
########  DEG analysis

dds$Groups <- relevel(dds$Groups, ref = "WT")
dds <- DESeq(dds)
resultsNames(dds)

### Generate a results table
res <- results(dds, alpha = 0.1, contrast = c("Groups", "Homo_KO", "WT"), 
               pAdjustMethod = 'fdr', independentFiltering = FALSE)
head(res) 
summary(res)


### convert the ID to gene symbol
{
  listEnsembl(version = 108)  # check the release 108
  ensembl = useEnsembl(biomart = 'genes', dataset='mmusculus_gene_ensembl', verbose = TRUE,
                       host="https://oct2022.archive.ensembl.org")  
  res <- as.data.frame(res)
  genes <- rownames(res)
  
  listFilters(ensembl)    
  listAttributes(ensembl)   
  annot <- getBM(attributes = c('ensembl_gene_id','mgi_symbol'),
                 filters = 'ensembl_gene_id',
                 values = genes,
                 verbose = TRUE,
                 mart = ensembl)
  names(annot)
  res$ID <- rownames(res)
  res <- merge(x = res,
               y =  annot,
               by.y = 'ensembl_gene_id',
               all.x = T,
               by.x = 'ID')
  
  colnames(res)
  res <- res[, c(1, ncol(res), 2:(ncol(res)-1))]
  head(res)
  
}

### define DEG threshold
deg <- subset(res, padj < 0.05 & abs(log2FoldChange) > 0.58)
deg <- deg[order(deg$padj), ]

write.csv(as.data.frame(res), paste0(rs_dir, "/DEG/DEGs_all.csv"))
write.csv(as.data.frame(deg), paste0(rs_dir, "/DEG/DEGs_sig.csv"))







