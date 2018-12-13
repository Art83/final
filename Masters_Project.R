# Final project script 

library(UpSetR)
library(GEOquery)
library(limma)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(ReactomePA)
library(reactome.db)
library(fgsea)

# Downloading GSE57465 to get entrez gene ids
gse57465_norm <- getGEO("GSE57465")[[1]]
# sometimes doesn't work, I think it's accounted for some technical issues with the internet connection and server response

# creating a data frame containing all gene IDs for this experiment
geneIDS <- data.frame(fData(gse57465_norm)$Entrez_Gene_ID)
rownames(geneIDS) <-  fData(gse57465_norm)$ID
colnames(geneIDS) <- "ENTREZ_GENE_ID"

# Downloading the whole dataset
gse57465 <- getGEO("GSE57465", GSEMatrix = F)

#creating a vector containing gene IDs 
probesets <- Table(GSMList(gse57465)[[6]])

# Organizing the GSE data frame for further analysis
data.matrix <- do.call('cbind',lapply(GSMList(gse57465),function(x)
{tab <- Table(x)
mymatch <- match(probesets[,1],tab$ID_REF)
return(tab$VALUE[mymatch])
}))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
rownames(data.matrix) <- probesets[,1]
colnames(data.matrix) <- names(GSMList(gse57465))


pdata <- data.frame(samples=names(GSMList(gse57465)))
rownames(pdata) <- names(GSMList(gse57465))
pheno <- as(pdata,"AnnotatedDataFrame")
es <- new('ExpressionSet',exprs=data.matrix,phenoData=pheno)
pData(es)$condition <- c(rep('two_hrs', 4), rep("half_hour", 4), rep("Home_cage", 4))


# working with the expression counts
# Normalization 
exprs(es) <- normalizeBetweenArrays(log2(exprs(es)+1), method="quantile")
pca <- prcomp(t(exprs(es)))
# PCA plot, doesn't look very neat, 3 samples from diferent groups look shifted,  apparently it has something to do with batch effects???
dataPCA <- data.frame(PC1 = pca$x[,1],PC2 = pca$x[,2], condition = pData(es)$condition)
qplot(PC1, PC2, data = dataPCA, color = condition) + geom_text(aes(label=row.names(dataPCA)))


# Stats
f = factor(pData(es)$condition,levels=c("Home_cage", "half_hour",'two_hrs'))
design = model.matrix(~ 0 + f)
colnames(design) <- c("Home_cage", "half_hour",'two_hrs')
data.fit = lmFit(exprs(es),design)
contrast.matrix = makeContrasts(half_hour-Home_cage,levels=design)
contrast.matrix.2h = makeContrasts(two_hrs-Home_cage, levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.con.2h = contrasts.fit(data.fit, contrast.matrix.2h)
data.fit.eb = eBayes(data.fit.con)
data.fit.eb.2h = eBayes(data.fit.con.2h)
tab = topTable(data.fit.eb, adjust.method = "BH", number = Inf)
tab2h = topTable(data.fit.eb.2h, adjust.method = "BH", number = Inf)

# Adding to the final result tables a column with the names of the genes
tab$symbol <- geneIDS[match(rownames(tab), rownames(geneIDS)), 1]
tab$name <- as.character(mapIds(org.Mm.eg.db, keys = as.character(tab$symbol), column = "SYMBOL",keytype = "ENTREZID"))

tab2h$symbol <- geneIDS[match(rownames(tab2h), rownames(geneIDS)), 1]
tab2h$name <- as.character(mapIds(org.Mm.eg.db, keys = as.character(tab2h$symbol), column = "SYMBOL",keytype = "ENTREZID"))

# Getting rid of all NAs in the final tables
tab <- tab[complete.cases(tab),]
tab2h <- tab2h[complete.cases(tab2h),]

# detecting up- and downregulated genes 153 overall 
tab %>% filter(logFC > 0 & adj.P.Val < 0.05) # 114
tab %>% filter(logFC < 0 & adj.P.Val < 0.05) # 39

#ranking according to foldcnange or t-statistics
tab <- tab[order(tab$logFC, decreasing = T),]
tab2h <- tab2h[order(tab2h$logFC, decreasing = T),]
#sorting by t
tab <- tab[order(tab$t, decreasing = T),]
tab2h <- tab2h[order(tab2h$t, decreasing = T),]

ranks30t <- setNames(tab$t, tab$name)
ranks2ht <- setNames(tab2h$t, tab2h$name) 

ranks30FC <- setNames(tab$logFC, tab$symbol)
ranks2hFC <- setNames(tab$logFC, tab$symbol)




# reactome enrichment
genes30m <- tab$symbol[tab$adj.P.Val<0.01]
xEnrich30 <- enrichPathway(gene=genes30m,pvalueCutoff=0.05, readable=T, organism = "mouse")
xEnrich30frame <- as.data.frame(xEnrich30)
barplot(xEnrich30, showCategory = 8)
enrichMap(xEnrich30, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
cnetplot(xEnrich30, showCategory = 5, categorySize="pvalue", foldChange=ranks30FC)



genes2h <- tab2h$symbol[tab2h$adj.P.Val<0.05]
xEnrich2h <- enrichPathway(gene=genes2h,pvalueCutoff=0.05, readable=T, organism = "mouse")
xEnrich2frame <- as.data.frame(xEnrich2h)
barplot(xEnrich2h, showCategory = 8)
enrichMap(xEnrich2h, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
cnetplot(xEnrich2h, showCategory = 5, categorySize="pvalue", foldChange=ranks30FC)



# looking for overlapping my results with the results from the paper

# extracting DEGs from my results
namesGenes30 <- tab %>%
  filter(adj.P.Val < 0.05) %>%
  select(name)
namesGenes30 <- as.vector(unlist(namesGenes30)) 

namesGenes2h <- tab2h %>%
  filter(adj.P.Val < 0.05) %>%
  select(name)
namesGenes2h <- as.vector(unlist(namesGenes2h))  

# reading the list of DEGs from the paper

ref_genes30_up <- scan("ref_up_30.txt", character(), quote = "")
ref_genes30_down <- scan("ref_down_30.txt", character(), quote = "")
ref_genes30 <- c(ref_genes30_up, ref_genes30_down)

ref_genes2h_up <- scan("ref_up_2h.txt", character(), quote = "")
ref_genes2h_down <- scan("ref_down_2h.txt", character(), quote = "")
ref_genes2h <- c(ref_genes2h_up, ref_genes2h_down)


# reading the data from GEO2r service

geo2r.30 <- read.table("geo2r_30.txt", header = T)
geo2r.2h <- read.table("geo2r_2hrs.txt", header = T)

# a lot of blank values, getting rid of them
geo2r.30 <- geo2r.30[!geo2r.30$Gene.symbol == '',]
geo2r.2h <- geo2r.2h[!geo2r.2h$Gene.symbol == '',]
length(geo2r.30[geo2r.30$adj.P.Val < 0.05,8]) # 320 DEGs!!
length(geo2r.2h[geo2r.2h$adj.P.Val < 0.05,8]) # 30 DEGs!!

# extracting the names of the DEGs in geo2r results
namesGenes30_geo2r <- geo2r.30 %>%
  filter(adj.P.Val < 0.05) %>%
  select(Gene.symbol)
namesGenes30_geo2r <- as.vector(unlist(namesGenes30_geo2r)) 

namesGenes2h_geo2r <- geo2r.2h %>%
  filter(adj.P.Val < 0.05) %>%
  select(Gene.symbol)
namesGenes2h_geo2r <- as.vector(unlist(namesGenes2h_geo2r))


# Looking for overlaps 

# Preparing the dataframe with all genes
Totalgenes <- c(namesGenes2h, namesGenes30,namesGenes2h_geo2r,namesGenes30_geo2r, ref_genes30, ref_genes2h)
Totalgenes <- unique(Totalgenes)
dataFr <- data.frame(Totalgenes, HalfH = 0, TwoH = 0, paper_HalfH = 0, paper_TwoH=0, geo2r_HalfH=0, geo2r_TwoH=0)
dataFr$HalfH <- ifelse(dataFr$Totalgenes %in% namesGenes30, 1, 0)
dataFr$TwoH <- ifelse(dataFr$Totalgenes %in% namesGenes2h, 1, 0)
dataFr$ref_30 <- ifelse(dataFr$Totalgenes %in% ref_genes30, 1, 0)
dataFr$ref_2h <- ifelse(dataFr$Totalgenes %in% ref_genes2h, 1, 0)
dataFr$geo2r_HalfH <- ifelse(dataFr$Totalgenes %in% namesGenes30_geo2r, 1, 0)
dataFr$geo2r_TwoH <- ifelse(dataFr$Totalgenes %in% namesGenes2h_geo2r, 1, 0)


# building a plot with overlaps
upset(dataFr)

# common genes
length(intersect(ref_genes30, namesGenes30)) # 23
length(intersect(ref_genes2h, namesGenes2h)) # 34
length(intersect(ref_genes30, namesGenes30_geo2r)) # 46
length(intersect(ref_genes2h, namesGenes2h_geo2r)) # 23
length(intersect(namesGenes2h, namesGenes2h_geo2r)) # 22
length(intersect(namesGenes30, namesGenes30_geo2r)) # 132


23/(length(ref_genes30) + length(namesGenes30))*100 # 9.38%
34/(length(ref_genes2h) + length(namesGenes2h))*100 # 18.27%
46/(length(ref_genes2h) + length(namesGenes2h_geo2r))*100 # 32.39%
23/(length(ref_genes30) + length(namesGenes30_geo2r))*100 # 5.58%
132 / (length(namesGenes30) + length(namesGenes30_geo2r))*100 # 27.91
22 / (length(namesGenes2h) + length(namesGenes2h_geo2r))*100 # 21.1538

# Common for 3 datasets
# 30
length(Reduce(intersect, list(ref_genes30,namesGenes30, namesGenes30_geo2r))) # 23 genes
# "Egr4"    "Arc"     "Ddit4"   "Dnajb1"  "Gad1"    "Per1"    "Hspa8"   "Fos"     "Nfkbia"  "Pdcd4"   "Sgk1"    "Fosb"   
#"Cox6c"   "Tsc22d3" "Plekhf1" "Il33"    "Acsl6"   "Mat2a"   "Sparc"   "Setd7"   "Dgkg"    "Lypd1"   "Gpr22" 

length(Reduce(intersect, list(ref_genes30,namesGenes30, namesGenes30_geo2r)))/(length(ref_genes30) +length(namesGenes30) + length(namesGenes30_geo2r))*100
# 4.07%

length(Reduce(intersect, list(ref_genes2h,namesGenes2h, namesGenes2h_geo2r))) # 18

# "Cox6c"   "Tsc22d3" "Plekhf1" "Il33"    "Acsl6"   "Mat2a"   "Sparc"   "Setd7"   "Dgkg"    "Lypd1"   "Gpr22" 

length(Reduce(intersect, list(ref_genes2h,namesGenes2h, namesGenes2h_geo2r))) / (length(ref_genes2h) + length(namesGenes2h) + length(namesGenes2h_geo2r))*100
# 8.33%
