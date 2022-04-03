library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)
library(goseq)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(goseq)

#Import count matrix and round the expected counts to integers
Counts <- read.table("output.genes.counts.matrix.CTRLBM1KD.txt", stringsAsFactors = FALSE)
Counts <- as.matrix(Counts)
Counts <- round(Counts)
Conditions <- c("Control", "Control", "BMI1KD", "BMI1KD")

#Perform differential expression analysis using EdgeR
DGE <- DGEList(Counts, group = Conditions)
keep <- filterByExpr(DGE) # Determine low count genes
summary(keep)
y <- DGE[keep, , keep.lib.sizes=FALSE] # Filter low count genes and recompute the library size
y <- calcNormFactors(y) # Calculate Normalization Factors
y <- estimateDisp(y)
plotBCV(y)
plotMDS(y, labels=rownames(DGE$samples$group), col=c("blue", "red") [factor(DGE$samples$group)]) # Plot MDS for sample distances

et <- exactTest(y, pair = c("Control", "BMI1KD"))
is.de <- decideTestsDGE(et, p.value = 0.01, lfc = 1) #Identify which genes are significantly differentially expressed
summary(is.de)
ResAll <- topTags(et, n=Inf, sort.by="none")$table
DE_genes <- rbind(ResAll[is.de==-1L,], ResAll[is.de==1L,])

#Converting ensemble gene ids to Gene symbols using biomaRt package
attributes <- c('ensembl_gene_id','gene_biotype','hgnc_symbol', 'external_gene_name', 'description','ensembl_transcript_id') #Attributes we want to retrieve from bioMart
mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
iv <- match(rownames(DE_genes),anno$ensembl_gene_id) # match the rownames of your resulting DEs with ensembl gene ids
DE_genes$gene <- anno[iv,'external_gene_name'] # Add a column to your dataframe with the gene SYMBOL, external_gene_id column in the anno object
write.csv(DE_genes, file = "DE_genesBMI1KDvsCTRL.csv", row.names = TRUE)

#Plotting
plotMD(et, status=is.de, values=c(1,-1), col=c("red","blue"), legend=FALSE); abline(h=c(-1, 1), col="blue") # MA plot with DE highlighted

nc = cpm(y, prior.count = 2, normalized.lib.sizes = TRUE, log = TRUE) # depth adjusted normalized counts per million (CPM) in log2 for HEATMAPS
ncDE = subset(nc, rownames(nc) %in% rownames(DE_genes))
pheatmap(ncDE, clustering_method = "ward.D2", show_colnames = TRUE, show_rownames = FALSE, drop_levels = TRUE, fontsize = 14, scale = "row") # Heatmap of ALL genes post-filtering

#Volcano plotting with Enhanced Volcano
EnhancedVolcano(ResAll,
                lab = rownames(ResAll),
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12, 12),
                pCutoff = 0.01,
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                legendPosition = 'bottom',
                selectLab = FALSE,
                title = "",
                subtitle = "",
                caption = "",
                legendVisible = FALSE,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1,
                FCcutoff = 1,)

#GoAnalysis using an R package called GOseq, which was developed for RNA-seq data
#First step is to retrieve the gene lenghts that will be used with goseq, obtain them from any RSEM results.
Control1est <- read.table("Control1.genes.results", stringsAsFactors = FALSE, header = TRUE)
GeneLength <- Control1est[,c(1,3)]
AllGenMeas <- rownames(Counts)
gene.vector <- as.integer(AllGenMeas%in%rownames(DE_genes)) #allmeasuredgenes from unfiltered counts chr vector and deg chr vector
names(gene.vector) = AllGenMeas
head(gene.vector)
pwf <- nullp(gene.vector,bias.data=GeneLength$length) # fitting probability weighting function. Weight based on gene length
GO.Wall <- goseq(pwf, "hg38", "ensGene") # Goseq command for the pwf in hg38 genome, gene names are ensGene
enriched.GO <- GO.Wall$category[p.adjust(GO.Wall$over_represented_pvalue, method="BH")<.05] # Subset based on FDR <0.05
SignTerms <- subset(GO.Wall, GO.Wall$category %in% enriched.GO) # Create data frame with enriched GO terms






