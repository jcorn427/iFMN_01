######### --Libraries--#########
setwd("/vol08/ngs/P51/iFMN/iFMN_01_Tissues/Cornelius_analysis/de_work/")
library(stringr)
library(limma)
library(edgeR)
library(PCAtools)
library(plyr)
library(dplyr)
library(tidyr)
library(WebGestaltR)
source("/vol08/ngs/P51/iFMN/iFMN_01_Tissues/Cornelius_analysis/de_work/heatmap3LW_function.R")

# Helper file
rhesus2human <- read.csv(
  file = "/vol08/ngs/P51/iFMN/iFMN_01_Tissues/Cornelius_analysis/de_work/Macaca_mulatta_Mmul_10.109/rhesus2human109v2.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# --Read in target files
message("STATUS: Load tables")
cm <- read.table("./count_matrix.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
target <- read.csv("./targetfile_slim.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

# Rename samples
newsampleIDs <- c()
for (i in colnames(cm)) {
  i <- str_remove(i, "_Lib\\d+\\S*$")
  i <- str_replace_all(i, "-", "_")
  newsampleIDs <- c(newsampleIDs, i)
}
colnames(cm) <- newsampleIDs

# --generate figure of all counts
png("de_intensities_raw_counts.png", res = 100)
par(xpd = TRUE)
if (length(rownames(target)) > 10) {
  plotDensities(log2(cm + 0.1), legend = FALSE)
} else {
  plotDensities(log2(cm + 0.1),
                legend = "topright",
                inset = c(-0.2, 0), levels(rownames(target))
  )
}
dev.off()

## Generate model matrix ##

target$bioreps <- paste(target$Tissue, target$Time_Point, sep = "_")
target$bioreps <- str_replace_all(target$bioreps, "-", "_")
bioreps <- factor(target$bioreps)
mm <- model.matrix(~ 0 + bioreps)

# id <- factor(target[, "Animal_ID"])
# tp <- factor(target[, "Time_Point"])
# age <- factor(target[, "Aged"])
# tr <- factor(target[, "Treatement"])
# mm <- model.matrix(~ 0 + id:tp:age)

rownames(mm) <- rownames(target)
colnames(mm) <- make.names(colnames(mm))
mm <- mm[, colnames(mm)[order(tolower(colnames(mm[, ])))]]
mm <- mm[, colSums(mm) > 0]

excludeAll <- nonEstimable(mm)
if (length(excludeAll) > 0) {
  message("WARNING: These samples are nonEstimatable, design matrix ", excludeAll)
}

if ("ti" %in% excludeAll) {
  return("interactions term non estimable")
}
mm <- mm[, !colnames(mm) %in% excludeAll]
if (!is.fullrank(mm)) {
  return("not full rank")
}


# Save model matrix for future analysis
saveRDS(mm, file = "model_matrix.RDS")
mm <- readRDS("model_matrix.RDS")

# -- Normalize data
# order target and count matrix so they are the same (THIS IS IMPORTANT)
cm <- cm[, rownames(target)]

# CHECK IF ORDER IS THE SAME
if (all.equal(colnames(cm), rownames(target)) != TRUE) {
  print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
  print(rownames(target))
  print(colnames(cm))
}

# normalize
cm2 <- DGEList(counts = cm)
cm2 <- calcNormFactors(cm2, method = "TMM") # TMM normalization
corfit <- duplicateCorrelation(cm2$counts, block = factor(target$Animal_ID))
png("mean_variance_norm.png")
Pi.CPM <- voom(counts = cm2, design = mm, correlation = corfit, normalize.method = "none", plot = T, span = 0.1)
dev.off()
write.csv(Pi.CPM$E, "Pi.CPM$E_all.csv")


# box plot of unfiltered data
png("boxplot_vnorm_all.png", width = 10, height = 8, units = "in", res = 100)
# par(mar=c(1,1,1,1))
minvalue <- min(Pi.CPM$E)
maxvalue <- max(Pi.CPM$E)
boxplot(Pi.CPM$E,
        labels = target$GaleID, ylim = c(minvalue - 1, maxvalue + 1),
        ylab = "voom expression", main = "Count matrix", cex.axis = .6, las = 2,
        frame = FALSE
)
dev.off()



#-- filter out genes from each group that are below mean count of 3 across samples 
#Iteratively adjusted thresholds and decided that 3 was the best cutoff to get rid of the 
# mean-variance hook shown in the initial mean-variance plot (iterative results not saved just ran in R)
A <- rowMeans(cm)
isexpr <- A >= 3
cmfl_counts <- cm[isexpr, ]
write.csv(cmfl_counts, "count_matrix_renamed_fl.csv")

#Normalize again with filtering
cm2 <- DGEList(counts = cmfl_counts)
cm2 <- calcNormFactors(cm2, method = "TMM") # TMM normalization
corfit <- duplicateCorrelation(cm2$counts, block = factor(target$Animal_ID))
png("mean_variance_norm_fl.png")
Pi.CPM <- voom(counts = cm2, design = mm, correlation = corfit, normalize.method = "none", plot = T, span = 0.1)
dev.off()
write.csv(Pi.CPM$E, "Pi.CPM$E_all_fl.csv")

# Save voom normalized object to resume analysis
saveRDS(Pi.CPM, file = "Pi.CPM.rds")
Pi.CPM <- readRDS("Pi.CPM.rds")

# Get sig genes with gene names
sig_HGNC <- merge(rhesus2human, Pi.CPM$E,
                  by.x = "Gene.stable.ID",
                  by.y = "row.names",
                  all.X = T, all.Y = T
)

sig_HGNC <- sig_HGNC[, !(names(sig_HGNC) %in% c("Gene.stable.ID"))]
sig_HGNC <- avereps(sig_HGNC,
                    ID = sig_HGNC$HGNC.symbol
)
rownames(sig_HGNC) <- sig_HGNC[, "HGNC.symbol"]
sig_HGNC <- sig_HGNC[, !(colnames(sig_HGNC) %in% c("HGNC.symbol"))]
sig_HGNC <- as.matrix(data.frame(sig_HGNC))
write.csv(sig_HGNC, "norm_matrix_HGNC_Aging_02.csv", quote = FALSE)


# box plot of filtered data
png("boxplot_vnorm_all_fl.png", width = 10, height = 8, units = "in", res = 100)
# par(mar=c(1,1,1,1))
minvalue <- min(Pi.CPM$E)
maxvalue <- max(Pi.CPM$E)
boxplot(Pi.CPM$E,
        labels = target$GaleID, ylim = c(minvalue - 1, maxvalue + 1),
        ylab = "voom normalized expression", main = "Normalized count matrix", cex.axis = .6, las = 2,
        frame = FALSE
)
dev.off()




#Feature reduction
p <- PCAtools::pca(Pi.CPM$E, 
                   metadata = target)

PCAtools::biplot(p, x='PC1', y='PC2', 
                 lab=NULL, 
                 shape='Time_Point',
                 colby='Shorthand_Animal_ID',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 10,
                 legendLabSize = 10,
                 legendTitleSize = 10,
                 showLoadings = FALSE,
                 legendIconSize = 3,
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)

screeplot(p, xlab = "Principal component", title = "SCREE plot", ylim = c(0,100), components = getComponents(p, 1:10),
          colBar = "dodgerblue")

PCAtools::biplot(p, x='PC1', y='PC2', 
                 lab=NULL, 
                 shape='Time_Point',
                 colby='Tissue',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 10,
                 legendLabSize = 10,
                 legendTitleSize = 10,
                 showLoadings = TRUE,
                 legendIconSize = 3,
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)



## Generate lmfit object ##

Pi.lmfit <- lmFit(Pi.CPM, design = mm, block = target$Animal_ID, correlation = corfit$consensus)
# saveRDS(Pi.lmfit, file = "lmfitobject.RDS")
#Pi.lmfit <- readRDS("lmfitobject.RDS")

contrastsmatrix <- c(
  "biorepsOB_lt_D03-biorepsOB_D_11",
  "biorepsOB_rt_D03-biorepsOB_D_11",
  "biorepsOB_lt_D45-biorepsOB_D_11",
  "biorepsOB_rt_D45-biorepsOB_D_11"
)

contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- eBayes(fit, robust = TRUE, trend = TRUE)
#fit2 <- treat(fit, lfc = 0.58, robust = TRUE, trend = TRUE)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.05)
sum_results  <- summary(results)
write.csv(sum_results,file = "sum_results.csv")

n <- 0
for (i in contrastsmatrix) {
  n <- n+1
  topTable_data <- topTable(fit2, coef=n, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
  
  topTable_data$id <- 1:nrow(topTable_data)
  
  topTable_data_genes <- merge(rhesus2human, topTable_data,
                               by.x = "Gene.stable.ID",
                               by.y = "row.names",
                               all.X = T, all.Y = T
  )
  
  topTable_data_genes <- topTable_data_genes[order(topTable_data_genes$id),]
  write.csv(topTable_data_genes, file = paste0(i, "_gene_table.csv"))
}




dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

ExpressMatrixde <- merge(rhesus2human, ExpressMatrixde,
                         by.x = "Gene.stable.ID",
                         by.y = "row.names",
                         all.X = T, all.Y = T
)

ExpressMatrixde$HGNC.symbol[ExpressMatrixde$HGNC.symbol==''] = ExpressMatrixde$Gene.stable.ID[ExpressMatrixde$HGNC.symbol ==""]
rownames(ExpressMatrixde) <- ExpressMatrixde[, "HGNC.symbol"]
ExpressMatrixde <- ExpressMatrixde[,c(3,4,5,6)]

ExpressMatrixde <- as.matrix(data.frame(ExpressMatrixde, check.names = FALSE))
ExpressMatrixde

class(ExpressMatrixde) <- "numeric"
collabels <- c(
  "Left Cheek D03 - D-11",
  "Right Cheek D03 - D-11",
  "Left Cheek D45 - D-11",
  "Right Cheek D45 - D-11"
)

colnames(ExpressMatrixde) <- collabels

hmap <- heatmap.L.4(ExpressMatrixde,
            figmargins = c(20, 10),
            cutoff = 1, distmethod = "euclidean", cexcol = 2,
            clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9, labRow = TRUE
)

###WebGestaltR analysis####
for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  #print(genes)
  WebGestaltR(enrichMethod="ORA",
              organism="hsapiens",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = rhesus2human[,2],
              referenceGeneType = "genesymbol",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  # write.csv(names(genes), paste0(cluster,"_genes_first_two.csv"))
  
  temp <- ExpressMatrixde[row.names(ExpressMatrixde)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs.csv"))
}




#### Look at loadings ####
loadings <- p$loadings

loadings <- merge(rhesus2human, loadings,
                  by.x = "Gene.stable.ID",
                  by.y = "row.names",
                  all.X = T, all.Y = T
)

write.csv(loadings, file = "loadings.csv")

# Generate PCA with gene names on the loadings
E <- Pi.CPM$E

E <- merge(rhesus2human, E,
                         by.x = "Gene.stable.ID",
                         by.y = "row.names",
                         all.X = T, all.Y = T
)

E$HGNC.symbol[E$HGNC.symbol==''] = E$Gene.stable.ID[E$HGNC.symbol ==""]
E$HGNC.symbol <- make.unique(E$HGNC.symbol)
rownames(E) <- E[, "HGNC.symbol"]
E <- subset(E, select = -c(Gene.stable.ID, HGNC.symbol))


q <- PCAtools::pca(E, 
                   metadata = target)

PCAtools::biplot(q, x='PC1', y='PC2', 
                 lab=NULL, 
                 shape='Time_Point',
                 colby='Shorthand_Animal_ID',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 10,
                 legendLabSize = 10,
                 legendTitleSize = 10,
                 showLoadings = FALSE,
                 legendIconSize = 3,
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)


PCAtools::biplot(q, x='PC1', y='PC2', 
                 lab=NULL, 
                 shape='Time_Point',
                 colby='Tissue',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 10,
                 legendLabSize = 10,
                 legendTitleSize = 10,
                 showLoadings = FALSE,
                 legendIconSize = 3,
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)

## Cibersort
cibersort <- Pi.CPM$E
cibersort <- merge(rhesus2human, cibersort,
                   by.x = "Gene.stable.ID",
                   by.y = "row.names",
                   all.X = T, all.Y = T
)
cibersort <- cibersort[, !(names(cibersort) %in% c("Gene.stable.ID"))]
cibersort <- avereps(cibersort,
                     ID = cibersort$HGNC.symbol
)
rownames(cibersort) <- cibersort[, "HGNC.symbol"]
cibersort <- cibersort[, !(colnames(cibersort) %in% c("HGNC.symbol"))]
cibersort <- as.matrix(data.frame(cibersort))
# cibersort <- make.unique(rownames(cibersort))
class(cibersort) <- "numeric"

cibersortnotlog <- 2^cibersort
cibersortnotlog <- cibersortnotlog[!(rownames(cibersortnotlog) %in% c("")), ]
write.table(cibersortnotlog, file.path("HGNC_cibersort_allnotlog.txt"), sep="\t", quote=F)
write.table(cibersort, file.path("HGNC_cibersort_all.txt"), sep="\t", quote=F)

## Generate tables for IPA analysis ##
topTable_data <- topTable(fit2, coef=1, resort.by=NULL, number=length(fit2$coefficients))
dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

topTable_data_sort <- topTable_data[match(rownames(dataMatrixde), rownames(topTable_data)),]
topTable_data_sig <- subset(topTable_data_sort, rowSums(sigMask) != 0)
topTable_data_sig_genes <- merge(rhesus2human, topTable_data_sig,
                                 by.x = "Gene.stable.ID",
                                 by.y = "row.names",
                                 all.X = T, all.Y = T
)

topTable_data_sig_genes[duplicated(topTable_data_sig_genes$Gene.stable.ID),]


n <- 0
top_big <- data.frame(matrix(NA, nrow = 38, ncol = 1))
for (i in contrastsmatrix) {
  n <- n+1
  topTable_data <- topTable(fit2, coef=n, resort.by=NULL, number=length(fit2$coefficients))
  
  topTable_data_sort <- topTable_data[match(rownames(dataMatrixde), rownames(topTable_data)),]
  topTable_data_sig <- subset(topTable_data_sort, rowSums(sigMask) != 0)
  topTable_data_sig_genes <- merge(rhesus2human, topTable_data_sig,
                                   by.x = "Gene.stable.ID",
                                   by.y = "row.names",
                                   all.X = T, all.Y = T
  )
  top_big <- cbind(top_big, topTable_data_sig_genes)
}

write.csv(top_big, file = "IPA_DEG_all.csv")



# Further analysis looking at individual animals
# Using EdgeR to get this. NOT STATISTICALLY SIGNIFICANT
group <- factor(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5))
group2 <- factor(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30))
y <- DGEList(counts = cm, group = group2)
keep <- filterByExpr(y)
y <- y[keep,]
y <- normLibSizes(y)
design <- model.matrix(~group2)
y <- estimateGLMCommonDisp(y, method = "deviance", robust = TRUE, subset = NULL)


fit_edgeR <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit_edgeR, coef = 2)
topTags(qlf)

fit_edgeR2 <- glmFit(y, design)
lrt <- glmLRT(fit_edgeR2, coef = 2)
topTags(lrt)
