library (limma)
library (edgeR)
library (openxlsx)


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



## OE paired contrast
pheno.s <- pheno[grep ("DOX", pheno$genotype), ]
pheno.s

pheno.s$batch <- factor (c(1,2,3,2,1,3)) 

counts <- a[ ,colnames (a) %in% pheno.s$sample]
stopifnot (colnames (counts) == pheno.s$sample)


## limma chunk

x <- DGEList(counts=counts) 

# filter for cpm > 6 in at least 3 samples
isexpr <- rowSums(cpm(x) > 6) >= 3
x <- x[isexpr, ]
dim (x$counts)


## paired limma test (the paired factor is treated as a batch factor)

celltype <- factor (pheno.s$genotype)
#colnames (x$counts)
batch <- factor (c (1,2,3,2,1,3))

design <- model.matrix (~ batch + celltype) 

v <- voomWithQualityWeights(x, design=design, plot=TRUE)
## CONTROLplusDOX_D0_14_S5 and PGBD5OEplusDOX_D0_14_S6 is down-weigthed by limma
vfit <- lmFit(v, design)
efit <- eBayes(vfit, trend=TRUE)

res <- topTable(efit,coef="celltypePGBD5OEplusDOX", n="inf")

boxplot (res$logFC)
abline (h=0)
abline (h=1)
abline (h=-1)

anno <- read.delim ("gencode.v43.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

res <- merge (res, anno, by.x="row.names", by.y="gene_id")
colnames (res)[1] <- "Geneid"

res[res$gene_name == "PGBD5", ]
#                 Geneid     logFC  AveExpr        t      P.Value    adj.P.Val
#8571 ENSG00000177614.11 0.7560978 5.819697 10.92561 1.205759e-06 4.362904e-05


## see individual paired log2 fold changes. These are not normalized values
norm.exprs <- v$E
norm.dif1 <- data.frame (D014= norm.exprs[ ,"PGBD5OEplusDOX_D0_14"] -  norm.exprs[ ,"CONTROLplusDOX_D0_14"])
norm.dif2 <- data.frame (D001= norm.exprs[ ,"PGBD5OEplusDOX_D0.1"] -  norm.exprs[ ,"CONTROLplusDOX_D01_2"])
norm.dif3 <- data.frame (D002= norm.exprs[ ,"PGBD5OEplusDOX_D0.2"] -  norm.exprs[ ,"CONTROLplusDOX_D02_2"])
norm.dif <- cbind (norm.dif1, norm.dif2, norm.dif3)

norm.dif$consistent <- "No"
norm.dif$consistent [apply (norm.dif[ ,1:3], 1, function (x) all (x > 0) )] <- "Up"
norm.dif$consistent [apply (norm.dif[ ,1:3], 1, function (x) all (x < 0) )] <- "Down"
table (norm.dif$consistent)
#Down   No   Up 
#2042 6381 2794 

res <- merge (res, norm.dif, by.x="Geneid", by.y="row.names")
res <- res[ ,c(1:11, 13:16, 12)]
res <- res[order (res$adj.P.Val), ]
write.xlsx (res, "deg_unfiltered_piggybac_overexpression_limma_new_pipeline.xlsx", rowNames=F)



## Sanity check (with the old pipeline: the old pipeline and the new IIT pipeline are highly correlated)

prev <- read.xlsx ("/Volumes/texas/iit_projects/devide/deg_unfiltered_piggybac_overexpression_limma.xlsx")
prev <- merge (res, prev, by.x="gene_name", by.y="Geneid")
plot (prev$logFC.x, prev$logFC.y, col=ifelse (prev$adj.P.Val.x < 0.05 & prev$adj.P.Val.y < 0.05, "darkblue", "black"), xlab="new_pipe_limma", ylab="prev_pipe_limma")
abline (h=0)
abline (v=0)
cor (prev$logFC.x, prev$logFC.y, method="pearson")
# 0.99












