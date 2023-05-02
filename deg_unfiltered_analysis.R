# 201110_A00558_0095_BHVNWGDMXX

library (openxlsx)
library (DESeq2)
library (ggplot2)
library (ggrepel)
library (pheatmap)


anno <- read.delim ("gencode.v43.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

a <- read.delim ("subread.counts.txt", skip=1)
a <- a[ ,grep ("ene|bam", colnames (a))]
a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("_S[0-9]+.*", "", colnames (a)) 
colnames (a) <- gsub ("star.IIT_RNM_", "", colnames (a))

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 
a <- a[ ,grep ("gene_type.y|gene_name.y", colnames (a), invert=TRUE)]
colnames (a) [colnames (a) == "gene_name.x"] <- "gene_name"
colnames (a) [colnames (a) == "gene_type.x"] <- "gene_type"

#write.xlsx (a, "star_gene_raw_counts.xlsx", rowNames=F)


annot <- a
annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "hgnc_id", "external_gene_name", "description")]

torm <- c("gene_name", "gene_type", "hgnc_id", "external_gene_name", "description")
a <- a[ ,!colnames (a) %in% torm]
row.names (a) <- a[ ,1]
colnames (a) <- gsub ("star.", "", colnames (a))
a <- a[ ,-1]


pheno <- data.frame (matrix (nrow=dim (a)[2], ncol=2))
colnames (pheno) <- c("sample", "genotype")
row.names (pheno) <- pheno$sample <- colnames (a)
pheno$genotype <- gsub ("_.*", "", pheno$sample)
pheno$genotype [pheno$genotype == "PGBD5"] <- "shPGBD5" 



## OE contrast

#pheno.s <- pheno[grep ("DOX", pheno$genotype), ]
#pheno.s

#a.s <- a[ ,colnames (a) %in% pheno.s$sample]
#stopifnot (colnames (a.s) == pheno.s$sample)

#dds <- DESeqDataSetFromMatrix(countData = round (a.s), colData = pheno.s, design = ~ genotype)

#keep <- rowSums(counts(dds) >= 50) >= 3
#dds <- dds[keep,]
#dds

#dds <- DESeq(dds)
#res <- results(dds, contrast=c("genotype", "PGBD5OEplusDOX", "CONTROLplusDOX"))

#res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
#res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
#colnames (res)[1] <- "Geneid"
#res <- res[order (res$padj), ]

## Sanity check
#res[res$gene_name == "PGBD5", ] 
## padj= 0.3889514 !!! not significant



## OE paired contrast
pheno.s <- pheno[grep ("DOX", pheno$genotype), ]
pheno.s

pheno.s$batch <- factor (c(1,2,3,2,1,3)) 

a.s <- a[ ,colnames (a) %in% pheno.s$sample]
stopifnot (colnames (a.s) == pheno.s$sample)

dds <- DESeqDataSetFromMatrix(countData = round (a.s), colData = pheno.s, design = ~ batch + genotype)

keep <- rowSums(counts(dds) >= 50) >= 3
dds <- dds[keep,]
dds

dds <- DESeq(dds)
res <- results(dds)
## MA plot
plotMA(res, ylim=c(-5,5))

# Wald test p-value: genotype PGBD5OEplusDOX vs CONTROLplusDOX 

res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

# Sanity check
res[res$gene_name == "PGBD5", ] 
# padj=0.017. significant

write.xlsx (res, "piggybac_PB_overexpressionDOX_vs_CTRL_DOX_human_in-vitro.xlsx", rowNames=F)


## PCA plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("genotype", "sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype, label=sample)) +
  		geom_point(size=3) +
  		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  		geom_text_repel()  + 
		  coord_fixed () 

ggsave ("PCA plot.pdf")



## Sample to sample correlation, see https://rockefelleruniversity.github.io/RU_RNAseq/presentations/singlepage/RU_RNAseq_p3.html

sampleDists <- as.dist(1 - cor(log2 (counts(dds,normalized=TRUE)+1), method="pearson"))
sampleDistMatrix <- as.matrix(sampleDists)

library(RColorBrewer)
blueColours <- brewer.pal(9, "Blues")
colors <- colorRampPalette(rev(blueColours))(255)

df <- as.data.frame(colData(dds)[,c("genotype","sample")])

pdf ("Distance between samples plot.pdf")
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, color = colors, annotation_row=df)
dev.off()




## sh contrast

pheno.s <- pheno[grep ("DOX", pheno$genotype, invert=TRUE), ]
pheno.s

a.s <- a[ ,colnames (a) %in% pheno.s$sample]
stopifnot (colnames (a.s) == pheno.s$sample)

dds <- DESeqDataSetFromMatrix(countData = round (a.s), colData = pheno.s, design = ~ genotype)

keep <- rowSums(counts(dds) >= 50) >= 3
dds <- dds[keep,]
dds

dds <- DESeq(dds)
res <- results(dds, contrast=c("genotype", "shPGBD5", "WT"))

res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

## Sanity check
res[res$gene_name == "PGBD5", ] 
# not significant


## Sample to sample correlation, see https://rockefelleruniversity.github.io/RU_RNAseq/presentations/singlepage/RU_RNAseq_p3.html

sampleDists <- as.dist(1 - cor(log2 (counts(dds,normalized=TRUE)+1), method="pearson"))
sampleDistMatrix <- as.matrix(sampleDists)

library(RColorBrewer)
blueColours <- brewer.pal(9, "Blues")
colors <- colorRampPalette(rev(blueColours))(255)

df <- as.data.frame(colData(dds)[,c("genotype","sample")])

pdf ("Distance between samples plot.pdf")
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, color = colors, annotation_row=df)
dev.off()






















