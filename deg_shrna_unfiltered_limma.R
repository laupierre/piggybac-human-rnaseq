## not working experiment

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



## shRNA paired contrast
pheno.s <- pheno[grep ("DOX", pheno$genotype, invert=TRUE), ]
pheno.s$genotype[pheno.s$genotype == "WT"] <- "CONTROL"

a <- a[ ,colnames (a) %in% pheno.s$sample]
stopifnot (colnames (a) == pheno.s$samples)

x <- DGEList(counts=a) 

# filter for cpm > 6 in at least 3 samples
isexpr <- rowSums(cpm(x) > 6) >= 3
x <- x[isexpr, ]
dim (x$counts)
# 11310     9

celltype <- factor (pheno.s$genotype)
batch <- factor (c (5,6,9,5,6,9,5,6,9))

design <- model.matrix (~ 0 + celltype + batch) 
colnames (design) <- gsub ("celltype", "", colnames (design))

contr.matrix <- makeContrasts (shpgbvsctrl = shPGBD5-CONTROL, shpgbdvsshctrl = shPGBD5-shCONTROL, shctrlvsctrl= shCONTROL-CONTROL, levels = colnames(design))

v <- voomWithQualityWeights(x, design=design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, trend=TRUE)

res1 <- topTable(efit,coef=1,sort.by="P", n= Inf)
res2 <- topTable(efit,coef=2,sort.by="P", n= Inf)
res3 <- topTable(efit,coef=3,sort.by="P", n= Inf)

res1[row.names (res1) == anno[anno$gene_name == "PGBD5", ]$gene_id, ]
#                        logFC  AveExpr         t    P.Value  adj.P.Val
#ENSG00000177614.11 -0.2473902 5.302327 -3.427104 0.00657502 0.01958996
res2[row.names (res2) == anno[anno$gene_name == "PGBD5", ]$gene_id, ]
#                        logFC  AveExpr         t     P.Value  adj.P.Val
#ENSG00000177614.11 -0.2639496 5.302327 -4.289607 0.001628369 0.03022514
res3[row.names (res3) == anno[anno$gene_name == "PGBD5", ]$gene_id, ]
#                        logFC  AveExpr         t  P.Value adj.P.Val
#ENSG00000177614.11 0.01655938 5.302327 0.2379178 0.816802 0.8655751


boxplot (res1$logFC)
abline (h=0)
abline (h=1)
abline (h=-1)

boxplot (res2$logFC)
abline (h=0)
abline (h=1)
abline (h=-1)


colnames (res1) <- paste (colnames (res1), "shpgbvsctrl", sep=".")
colnames (res2) <- paste (colnames (res2), "shpgbvsshctrl", sep=".")
colnames (res3) <- paste (colnames (res3), "shctrolvsctrl", sep=".")

resall <- merge (res1, res2, by="row.names")
colnames (resall)[1] <- "gene_id"
resall <- merge (resall, res3, by.x="gene_id", by.y="row.names")
colnames (resall)[1] <- "gene_id"
res <- merge (resall, anno, by="gene_id")
res <- res[order (res$adj.P.Val.shpgbvsctrl), ]

#res[res$gene_name == "PGBD5", ]

## There are more genes in shPGBD5 vs CTRL (than vs shCTRL) and as much in the control contrast !
table (res$adj.P.Val.shpgbvsctrl < 0.05)
#FALSE  TRUE 
# 6044  5266 
table (res$adj.P.Val.shpgbvsshctrl < 0.05)
#FALSE  TRUE 
#10278  1032 
table (res$adj.P.Val.shctrolvsctrl < 0.05)
#FALSE  TRUE 
# 6974  4336

write.xlsx (res, "deg_unfiltered_piggybac_shrna_hesc_limma_new_pipeline.xlsx", rowNames=F)



## Comparison of contrasts
par (mfrow=c(2,1))
res.s <- res[res$adj.P.Val.shpgbvsctrl < 0.05, ]

plot (res.s$logFC.shpgbvsctrl, res.s$logFC.shpgbvsshctrl, xlab="shPGBD5 vs CTRL", ylab="shPGBD5 vs shCTRL", main="significant shPGBD5 vs CTRL genes", xlim=c(-3,3), ylim=c(-3,3))
abline (h=0)
abline (v=0)
abline (0,1, col="red")

plot (res.s$logFC.shpgbvsctrl, res.s$logFC.shctrolvsctrl, xlab="shPGBD5 vs CTRL", ylab="shCTRL vs CTRL", main="significant shPGBD5 vs CTRL genes", xlim=c(-3,3), ylim=c(-3,3))
abline (h=0)
abline (v=0)
abline (0,1, col="red")


par (mfrow=c(2,1))
res.s <- res[res$adj.P.Val.shpgbvsshctrl < 0.05, ]

plot (res.s$logFC.shpgbvsctrl, res.s$logFC.shpgbvsshctrl, xlab="shPGBD5 vs CTRL", ylab="shPGBD5 vs shCTRL", main="significant shPGBD5 vs shCTRL genes", xlim=c(-3,3), ylim=c(-3,3))
abline (h=0)
abline (v=0)
abline (0,1, col="red")

plot (res.s$logFC.shpgbvsctrl, res.s$logFC.shctrolvsctrl, xlab="shPGBD5 vs CTRL", ylab="shCTRL vs CTRL", main="significant shPGBD5 vs shCTRL genes", xlim=c(-3,3), ylim=c(-3,3))
abline (h=0)
abline (v=0)
abline (0,1, col="red")











