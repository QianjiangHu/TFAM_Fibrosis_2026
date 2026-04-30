
################################################################
# this is the code analyzing the RNA-seq data of in vitro TFAM-deletion AT2 cells
# main tools: `DESeq2`
# author: Qianjiang Hu
# date: 2026-03-26
################################################################

# Load the WGCNA library
library(WGCNA)
# Enable parallel processing
options(stringsAsFactors = FALSE)
disableWGCNAThreads()


#===============================================================================
#
# Step 1: Data input and cleaning
#
#===============================================================================

vsd <- vst(dds, blind = FALSE)
mtx_df <- assay(vsd) 

{
  head(mtx_df)
  dim(mtx_df)
  names(mtx_df)
  coldata
  datExpr <- as.data.frame(t(mtx_df))

  sampleInfo <- as.data.frame(coldata$Groups)
  rownames(sampleInfo) <- coldata$SampleID
  colnames(sampleInfo) <- "Genotype"
  
  # Convert the Genotype to a numeric factor for WGCNA
  sampleInfo$Genotype <- as.factor(sampleInfo$Genotype)
  levels(sampleInfo$Genotype) <- c(1,0) 
  sampleInfo$Genotype <- as.numeric(as.character(sampleInfo$Genotype))
  
  # Check the sample trait dataframe
  print(sampleInfo)
  
  # Ensure that the sample names in 'expr_mat' match the row names in 'sampleInfo'
  all(rownames(datExpr) == rownames(sampleInfo))
}


#===============================================================================
#
# Step 2: Step-by-step network construction
#
#===============================================================================

# -----------------------------------------------
#  step 2.1 Choosing the soft-thresholding power
# -----------------------------------------------

# Choose a set of soft threshold parameters
powers = c(c(1:50), seq(from = 22, to=30, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 

par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") 
abline(h=0.80,col="blue") 
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

abline(h=300,col="red") 
dev.off()


# ------------------------------------------------------------------------------
#  step 2.2 Co-expression similarity and adjacency
# ------------------------------------------------------------------------------

power = 30
adjacency = adjacency(datExpr, type = "signed", power = power)

# ------------------------------------------------------------------------------
#  step 2.3 Topological Overlap Matrix (TOM)
# ------------------------------------------------------------------------------

TOM = TOMsimilarity(adjacency, TOMType = "signed"); # Turn adjacency into topological overlap
dissTOM = 1-TOM

# ------------------------------------------------------------------------------
#  step 2.4 Construct clustering using TOM
# ------------------------------------------------------------------------------

# Plot e hierarchical clustering tree
geneTree = hclust(as.dist(dissTOM), method = "average");
#pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()


# Identify Modules
# Module identification using dynamic tree cut:
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# ------------------------------------------------------------------------------
#  step 2.5 Merging of modules whose expression profiles are very similar
# ------------------------------------------------------------------------------


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;


pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

table(mergedColors)

#===============================================================================
#
# Step 3: Relating modules to Phenotype
#
#===============================================================================
colnames(mergedMEs)
colnames(MEs0)

# ------------------------------------------------------------------------------
#  step 3.1 Quantifying module–trait associations
# ------------------------------------------------------------------------------


# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, sampleInfo, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


# display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

# Step 1: Combine correlation and p-values into a single data frame
moduleTrait_df <- as.data.frame(moduleTraitCor)
moduleTrait_df$pvalue <- moduleTraitPvalue[, 1]  # Assuming single trait: Genotype
moduleTrait_df$Module <- rownames(moduleTrait_df)

# Step 2: Sort by correlation (descending)
sorted_df <- moduleTrait_df[order(moduleTrait_df$Genotype), ]

# Step 3: Rebuild the correlation matrix and text labels
sorted_cor_matrix <- as.matrix(sorted_df$Genotype)
rownames(sorted_cor_matrix) <- sorted_df$Module

sorted_text_matrix <- paste0(signif(sorted_df$Genotype, 2), "\n(p=", signif(sorted_df$pvalue, 1), ")")
dim(sorted_text_matrix) <- dim(sorted_cor_matrix)

# Step 4: Save to PNG (or use pdf() for vector graphics)
save.path <- file.path(pl_dir, "WGCNA", "Sorted_ModuleTrait_Heatmap_ggplot.png")
png(save.path, width = 250, height = 1000)
par(mar = c(6, 8.5, 3, 3))
# Step 5: Draw the heatmap
labeledHeatmap(Matrix = sorted_cor_matrix,
               xLabels = "Genotype",
               yLabels = rownames(sorted_cor_matrix),
               ySymbols = rownames(sorted_cor_matrix),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = sorted_text_matrix,
               setStdMargins = TRUE,
               cex.text = 0.7,
               zlim = c(-1,1),
               verticalSeparator.ext = 0,
               main = "Module–Genotype \nCorrelations")

dev.off()

save.path <- file.path(pl_dir, "WGCNA", "Sorted_ModuleTrait_Heatmap_ggplot.pdf")
pdf(save.path, width = 3.5, height = 12)
par(mar = c(6, 9, 3, 3))              # bottom, left, top, right

labeledHeatmap(
  Matrix      = sorted_cor_matrix,
  xLabels     = "Genotype",
  yLabels     = rownames(sorted_cor_matrix),
  ySymbols    = rownames(sorted_cor_matrix),
  colorLabels = FALSE,
  colors      = blueWhiteRed(50),
  textMatrix  = sorted_text_matrix,
  setStdMargins = FALSE,                  # use our par(mar)
  cex.text    = 0.7,
  zlim        = c(-1, 1),
  verticalSeparator.ext = 0,
  main        = "Module-Genotype\ncorrelations"  # plain hyphen
)

dev.off()


# ------------------------------------------------------------------------------
#  step 3.2 Gene Significance (GS) and Module Memberships (MM):  Gene relationship to trait and important modules
# ------------------------------------------------------------------------------


sampleInfo
# Define variable weight containing the weight column of datTrait
Genotype = as.data.frame(sampleInfo$Genotype)
names(Genotype) = "Genotype"

geneModuleMembership <- signedKME(datExpr, MEs)
colnames(geneModuleMembership) <- paste("MM", gsub("^ME", "", names(MEs)), sep="")

### Caculate p-value of module membership
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(MMPvalue) <- paste("p.MM", gsub("^ME", "", names(MEs)), sep="")

### Caculate Gene Significance (GS) for Genotype and it's p-value
geneTraitSignificance = as.data.frame(cor(datExpr, Genotype, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Genotype), sep="")
names(GSPvalue) = paste("p.GS.", names(Genotype), sep="")

# Calculate gene connectivity (degree) - number of connections each gene has
geneConnectivity <- rowSums(adjacency)

# Create a gene-to-connectivity dataframe
gene_names <- rownames(mtx_df)
gene_connectivity_assignment <- as.data.frame(cbind(gene = gene_names,
                                                    connectivity = geneConnectivity))

# ------------------------------------------------------------------------------
#  step 3.3 Intramodular analysis: identifying genes with high GS and MM
# ------------------------------------------------------------------------------


module = "blue"
column = match(module, modNames)
moduleGenes = moduleColors==module;

save.path <- file.path(pl_dir, "WGCNA", paste0("MM_vs_GS_", module, ".png"))
png(save.path, width = 600, height = 600)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Genotype",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

dev.off()







