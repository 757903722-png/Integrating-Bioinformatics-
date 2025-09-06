# Function to create boxplots for data visualization
bioBoxplot <- function(inputFile = NULL, outFile = NULL, titleName = NULL) {
  # Read data file
  rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  data <- t(rt)
  
  # Extract project and sample information
  Project <- gsub("(.*?)\\_.*", "\\1", rownames(data))
  Sample <- gsub("(.+)\\_(.+)\\_(.+)", "\\2", rownames(data))
  data <- cbind(as.data.frame(data), Sample, Project)
  
  # Reshape data for ggplot2
  rt1 <- melt(data, id.vars = c("Project", "Sample"))
  colnames(rt1) <- c("Project", "Sample", "Gene", "Expression")
  
  # Create boxplot
  pdf(file = outFile, width = 10, height = 5)
  p <- ggplot(rt1, mapping = aes(x = Sample, y = Expression)) +
    geom_boxplot(aes(fill = Project), notch = TRUE, outlier.shape = NA) +
    ggtitle(titleName) + theme_bw() + theme(panel.grid = element_blank()) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 2), 
          plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}

# Create boxplots before and after batch correction
bioBoxplot(inputFile = "merge.preNorm.txt", outFile = "boxplot.preNorm.pdf", 
           titleName = "Before batch correction")
bioBoxplot(inputFile = "merge.normalize.txt", outFile = "boxplot.normalize.pdf", 
           titleName = "After batch correction")
# Function for PCA analysis
bioPCA <- function(inputFile = NULL, outFile = NULL, titleName = NULL) {
  # Read data file
  rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  data <- t(rt)
  
  # Extract project information
  Project <- gsub("(.*?)\\_.*", "\\1", rownames(data))
  
  # Perform PCA
  data.pca <- prcomp(data)
  pcaPredict <- predict(data.pca)
  PCA <- data.frame(PC1 = pcaPredict[, 1], PC2 = pcaPredict[, 2], Type = Project)
  
  # Create PCA plot
  pdf(file = outFile, width = 5.5, height = 4.25)
  p1 <- ggscatter(data = PCA, x = "PC1", y = "PC2", color = "Type", shape = "Type", 
                  ellipse = TRUE, ellipse.type = "norm", ellipse.border.remove = FALSE, 
                  ellipse.alpha = 0.1, size = 2, main = titleName, legend = "right") +
    theme(plot.margin = unit(rep(1.5, 4), 'lines'), plot.title = element_text(hjust = 0.5))
  print(p1)
  dev.off()
}

# Perform PCA before and after batch correction
bioPCA(inputFile = "merge.preNorm.txt", outFile = "PCA.preNorm.pdf", 
       titleName = "Before batch correction")
bioPCA(inputFile = "merge.normalize.txt", outFile = "PCA.normalize.pdf", 
       titleName = "After batch correction")

# Set parameters for differential expression analysis
logFCfilter <- 0.585
adj.P.Val.Filter <- 0.05
inputFile <- "merge.normalize.txt"

# Read and preprocess data
rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

# Extract sample information
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data <- data[, order(Type)]
Project <- gsub("(.+)\\_(.+)\\_(.+)", "\\1", colnames(data))
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
colnames(data) <- gsub("(.+)\\_(.+)\\_(.+)", "\\2", colnames(data))

# Design matrix for differential expression
design <- model.matrix(~0 + factor(Type))
colnames(design) <- c("Control", "Treat")
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(Treat - Control, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Extract differential expression results
allDiff <- topTable(fit2, adjust = 'fdr', number = 200000)
allDiffOut <- rbind(id = colnames(allDiff), allDiff)
write.table(allDiffOut, file = "all.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# Extract significant differential expression results
diffSig <- allDiff[with(allDiff, (abs(logFC) > logFCfilter & adj.P.Val < adj.P.Val.Filter)), ]
diffSigOut <- rbind(id = colnames(diffSig), diffSig)
write.table(diffSigOut, file = "diff.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# Extract expression of differentially expressed genes
diffGeneExp <- data[row.names(diffSig), ]
diffGeneExpOut <- rbind(id = paste0(colnames(diffGeneExp), "_", Type), diffGeneExp)
write.table(diffGeneExpOut, file = "diffGeneExp.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# Create heatmap of top differentially expressed genes
geneNum <- 50
diffUp <- diffSig[diffSig$logFC > 0, ]
diffDown <- diffSig[diffSig$logFC < 0, ]
geneUp <- row.names(diffUp)
geneDown <- row.names(diffDown)
if (nrow(diffUp) > geneNum) { geneUp <- row.names(diffUp)[1:geneNum] }
if (nrow(diffDown) > geneNum) { geneDown <- row.names(diffDown)[1:geneNum] }
hmExp <- data[c(geneUp, geneDown), ]

# Prepare annotation data
names(Type) <- colnames(data)
Type <- as.data.frame(Type)
Type <- cbind(Project, Type)

# Create heatmap
pdf(file = "heatmap.pdf", width = 10, height = 7)
pheatmap(hmExp, 
         annotation_col = Type, 
         color = colorRampPalette(c("blue2", "white", "red2"))(50),
         cluster_cols = FALSE,
         show_colnames = FALSE,
         scale = "row",
         fontsize = 8,
         fontsize_row = 5.5,
         fontsize_col = 8)
dev.off()

# WGCNA analysis function
expFile <- "merge.normalize.txt"

# Read and preprocess data
rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)
data <- data[apply(data, 1, sd) > 0.5, ]

# Extract sample information
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
conCount <- length(Type[Type == "Control"])
treatCount <- length(Type[Type == "Treat"])
datExpr0 <- t(data)

# Check for missing values
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Sample clustering to detect outliers
sampleTree <- hclust(dist(datExpr0), method = "average")
pdf(file = "1_sample_cluster.pdf", width = 9, height = 6)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 20000, col = "red")
dev.off()

# Remove outliers
clust <- cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)
table(clust)
keepSamples <- (clust == 1)
datExpr0 <- datExpr0[keepSamples, ]

# Prepare trait data
traitData <- data.frame(Control = c(rep(1, conCount), rep(0, treatCount)),
                        Treat = c(rep(0, conCount), rep(1, treatCount)))
row.names(traitData) <- colnames(data)
fpkmSamples <- rownames(datExpr0)
traitSamples <- rownames(traitData)
sameSample <- intersect(fpkmSamples, traitSamples)
datExpr0 <- datExpr0[sameSample, ]
datTraits <- traitData[sameSample, ]

# Sample clustering with trait heatmap
sampleTree2 <- hclust(dist(datExpr0), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)
pdf(file = "2_sample_heatmap.pdf", width = 9, height = 7)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Soft threshold selection
enableWGCNAThreads()
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file = "3_scale_independence.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.8, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# Determine soft threshold power
sft
softPower <- sft$powerEstimate
adjacency <- adjacency(datExpr0, power = softPower)
softPower

# TOM matrix
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Gene clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")
pdf(file = "4_gene_clustering.pdf", width = 8, height = 6)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# Module identification
minModuleSize <- 60
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
pdf(file = "5_Dynamic_Tree.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Module eigengenes
MEList <- moduleEigengenes(datExpr0, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
pdf(file = "6_Clustering_module.pdf", width = 7, height = 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres <- 0.25
abline(h = MEDissThres, col = "red")
dev.off()

# Merge similar modules
merge <- mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
pdf(file = "7_merged_dynamic.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, mergedColors, "Merged dynamic",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors <- mergedColors
table(moduleColors)
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs

# Module-trait relationships
nGenes <- ncol(datExpr0)
nSamples <- nrow(datExpr0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
pdf(file = "8_Module_trait.pdf", width = 5.5, height = 5.5)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.75,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))
dev.off()

# Calculate module membership and gene significance
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")
traitNames <- names(datTraits)
geneTraitSignificance <- as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", traitNames, sep = "")
names(GSPvalue) <- paste("p.GS.", traitNames, sep = "")

# Module significance
y <- datTraits[, 1]
GS1 <- as.numeric(cor(y, datExpr0, use = "p"))
GeneSignificance <- abs(GS1)
ModuleSignificance <- tapply(GeneSignificance, mergedColors, mean, na.rm = TRUE)
pdf(file = "9_GeneSignificance.pdf", width = 12.5, height = 7.5)
plotModuleSignificance(GeneSignificance, mergedColors)
dev.off()

# Scatter plots of module membership vs gene significance
trait <- "Treat"
traitColumn <- match(trait, traitNames)
for (module in modNames) {
  column <- match(module, modNames)
  moduleGenes <- moduleColors == module
  if (nrow(geneModuleMembership[moduleGenes, ]) > 1) {
    outPdf <- paste("10_", trait, "_", module, ".pdf", sep = "")
    pdf(file = outPdf, width = 7, height = 7)
    par(mfrow = c(1, 1))
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, traitColumn]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for ", trait),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    dev.off()
  }
}

# Create gene information table
probes <- colnames(datExpr0)
geneInfo0 <- data.frame(probes = probes,
                        moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance)) {
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneTraitSignificance[, Tra],
                          GSPvalue[, Tra])
  names(geneInfo0) <- c(oldNames, names(geneTraitSignificance)[Tra],
                        names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[, mod],
                          MMPvalue[, mod])
  names(geneInfo0) <- c(oldNames, names(geneModuleMembership)[mod],
                        names(MMPvalue)[mod])
}
geneOrder <- order(geneInfo0$moduleColor)
geneInfo <- geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.xls", sep = "\t", row.names = FALSE)

# Export genes in each module
for (mod in 1:nrow(table(moduleColors))) {
  modules <- names(table(moduleColors))[mod]
  probes <- colnames(datExpr0)
  inModule <- (moduleColors == modules)
  modGenes <- probes[inModule]
  write.table(modGenes, file = paste0("module_", modules, ".txt"), sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Venn diagram analysis
diffFile <- "diff.txt"
moduleFile <- "module_red.txt"

# Read differential expression results
rt <- read.table("diff.txt", header = TRUE, sep = "\t", check.names = FALSE)
geneNames <- as.vector(rt[, 1])
geneNames <- gsub("^ | $", "", geneNames)
uniqGene <- unique(geneNames)
geneList <- list()
geneList[["DEG"]] <- uniqGene

# Read module genes
rt <- read.table(moduleFile, header = FALSE, sep = "\t", check.names = FALSE)
geneNames <- as.vector(rt[, 1])
geneNames <- gsub("^ | $", "", geneNames)
uniqGene <- unique(geneNames)
geneList[["WGCNA"]] <- uniqGene

# Create Venn diagram
venn.plot <- venn.diagram(geneList, filename = NULL, fill = c("cornflowerblue", "darkorchid1"),
                          scaled = FALSE, cat.pos = c(-1, 1), 
                          cat.col = c("cornflowerblue", "darkorchid1"), cat.cex = 1.2)
pdf(file = "venn.pdf", width = 5, height = 5)
grid.draw(venn.plot)
dev.off()

# Extract intersecting genes
interGenes <- Reduce(intersect, geneList)
write.table(interGenes, file = "interGenes.txt", sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)

# GO enrichment analysis
pvalueFilter <- 0.05
adjPvalFilter <- 1

# Set color scheme
colorSel <- "p.adjust"
if (adjPvalFilter > 0.05) {
  colorSel <- "pvalue"
}

# Read gene list
rt <- read.table("interGenes1.txt", header = FALSE, sep = "\t", check.names = FALSE)

# Convert gene symbols to Entrez IDs
genes <- unique(as.vector(rt[, 1]))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
rt <- cbind(rt, entrezIDs)
rt <- rt[rt[, "entrezIDs"] != "NA", ]
gene <- rt$entrezID

# Perform GO enrichment
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, 
               qvalueCutoff = 1, ont = "all", readable = TRUE)
GO <- as.data.frame(kk)
GO <- GO[(GO$pvalue < pvalueFilter & GO$p.adjust < adjPvalFilter), ]
write.table(GO, file = "GO.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Create barplot
pdf(file = "barplot.pdf", width = 9, height = 7)
bar <- barplot(kk, drop = TRUE, showCategory = 10, label_format = 100, 
               split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bar)
dev.off()

# Create bubble plot
pdf(file = "bubble.pdf", width = 9, height = 7)
bub <- dotplot(kk, showCategory = 10, orderBy = "GeneRatio", label_format = 100, 
               split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bub)
dev.off()

# KEGG enrichment analysis
pvalueFilter <- 0.05
adjPvalFilter <- 1

# Set color scheme
colorSel <- "p.adjust"
if (adjPvalFilter > 0.05) {
  colorSel <- "pvalue"
}

# Read gene list
rt <- read.table("interGenes1.txt", header = FALSE, sep = "\t", check.names = FALSE)

# Convert gene symbols to Entrez IDs
genes <- unique(as.vector(rt[, 1]))
entrezIDs <- mget(genes, org.Mm.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
rt <- cbind(rt, entrezIDs)
colnames(rt)[1] <- "gene"
rt <- rt[rt[, "entrezIDs"] != "NA", ]
gene <- rt$entrezID

# Perform KEGG enrichment
kk <- enrichKEGG(gene = gene, organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1)
kk@result$Description <- gsub(" - Homo sapiens \\(human\\)", "", kk@result$Description)
KEGG <- as.data.frame(kk)
KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x) {
  paste(rt$gene[match(strsplit(x, "/")[[1]], as.character(rt$entrezID))], collapse = "/")
}))
KEGG <- KEGG[(KEGG$pvalue < pvalueFilter & KEGG$p.adjust < adjPvalFilter), ]
write.table(KEGG, file = "KEGG.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Determine number of pathways to show
showNum <- 30
if (nrow(KEGG) < showNum) {
  showNum <- nrow(KEGG)
}

# Create barplot
pdf(file = "barplot.pdf", width = 8.5, height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, label_format = 100, color = colorSel)
dev.off()

# Create bubble plot
pdf(file = "bubble.pdf", width = 8.5, height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", label_format = 100, color = colorSel)
dev.off()

# GSEA analysis
gene <- "FIBP"
expFile <- "merge.normalize.txt"
gmtFile <- "c2.cp.kegg.Hs.symbols.gmt"

# Read and preprocess data
rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)
data <- data[rowMeans(data) > 0, ]

# Filter for treatment samples
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data <- data[, Type == "Treat", drop = FALSE]

# Calculate logFC for target gene
dataL <- data[, data[gene, ] < median(data[gene, ]), drop = FALSE]
dataH <- data[, data[gene, ] >= median(data[gene, ]), drop = FALSE]
meanL <- rowMeans(dataL)
meanH <- rowMeans(dataH)
meanL[meanL < 0.00001] <- 0.00001
meanH[meanH < 0.00001] <- 0.00001
logFC <- meanH - meanL
logFC <- sort(logFC, decreasing = TRUE)
genes <- names(logFC)

# Read gene set file
gmt <- read.gmt(gmtFile)

# Perform GSEA
kk <- GSEA(logFC, TERM2GENE = gmt, pvalueCutoff = 1)
kkTab <- as.data.frame(kk)
kkTab <- kkTab[kkTab$pvalue < 0.05, ]
write.table(kkTab, file = "GSEA.result.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Create enrichment plots for high expression group
termNum <- 5
kkUp <- kkTab[kkTab$NES > 0, ]
if (nrow(kkUp) >= termNum) {
  showTerm <- row.names(kkUp)[1:termNum]
  gseaplot <- gseaplot2(kk, showTerm, base_size = 8, title = "Enriched in high expression group")
  pdf(file = "GSEA.highExp.pdf", width = 6.5, height = 5.5)
  print(gseaplot)
  dev.off()
}

# Create enrichment plots for low expression group
termNum <- 5
kkDown <- kkTab[kkTab$NES < 0, ]
if (nrow(kkDown) >= termNum) {
  showTerm <- row.names(kkDown)[1:termNum]
  gseaplot <- gseaplot2(kk, showTerm, base_size = 8, title = "Enriched in low expression group")
  pdf(file = "GSEA.lowExp.pdf", width = 6.5, height = 5.5)
  print(gseaplot)
  dev.off()
}

# Machine learning modeling function
RunRF <- function(Train_set, train_label, mode = "XXX") {
  # Random forest implementation
  # [Implementation details would go here]
}

# Load required libraries
library(openxlsx)

# Read method list
methods <- read.xlsx("./data/Machine_learning_methods.xlsx")$Model
methods <- gsub("-| ", "", methods)

# Pre-training feature extraction
min.selected.var <- 2
Variable <- colnames(Train_set)
preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method <- unique(unlist(preTrain.method))

# Feature selection using various methods
set.seed(seed = 123)
preTrain.var <- list()
for (method in preTrain.method) {
  preTrain.var[[method]] <- RunML(method = method,
                                  Train_expr = Train_set,
                                  Train_surv = train_label,
                                  mode = "Variable")
}
preTrain.var[["simple"]] <- colnames(Train_set)

# Model training
model <- list()
set.seed(seed = 123)
for (method in methods) {
  cat(match(method, methods), ":", method, "\n")
  method_name <- method
  method <- strsplit(method, "\\+")[[1]]
  
  if (length(method) == 1) method <- c("simple", method)
  selected.var <- preTrain.var[[method[1]]]
  
  if (length(selected.var) <= min.selected.var) {
    model[[method_name]] <- NULL
  } else {
    model[[method_name]] <- RunML(method = method[2],
                                  Train_expr = Train_set[, selected.var],
                                  Train_surv = train_label,
                                  mode = "Model")
  }
}

# Extract features from models
fea_df <- lapply(model, function(fit) {
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, "./result/fea_df.xls", sep = "\t", row.names = FALSE, 
            col.names = TRUE, quote = FALSE)

