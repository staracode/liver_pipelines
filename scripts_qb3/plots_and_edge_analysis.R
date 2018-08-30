# Goal: 
#compare affected to unaffected

library(pheatmap)
library(edgeR)

#scatterplot
#pca highest variance genes
#dendrogram (why does it algn with tsne clusters????)
#DE

filename = "~/Downloads/readcounts.txt"
counts1 = read.table(filename, header = T, stringsAsFactors = F) 
counts2 = counts1[,7:ncol(counts1)]
#keep = which(apply(counts2, 1, function(x) sum(which(x>0))==ncol(counts2)))
keep = rowMeans(counts2) > 5 & rowMeans(counts2) < 5000 

counts3 = counts2[keep,]
counts3 = counts3 + 1
y <- DGEList(counts=counts3)
y <- calcNormFactors(y)

counts4 = rpkm (y, gene.length=counts1[keep,"Length"], log = TRUE)

#scatterplot pairwise

#dendrogram 
hc <- hclust(dist(t(counts4)))
dend1 <- as.dendrogram(hc)
plot(dend1)
library(DESeq2)

cold = colnames(counts4)
#condition=
mmm = counts4[1:400,]
ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData=mmm, colData=cold);
#ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData=mmm, colData=cold, design=~condition);
colnames(ddsMat) <- colnames(mmm)
rld <- DESeq2::rlog(ddsMat, blind=blind_status) # LOG TRANSFORMED.

HEATMAP_HEIGHT <- 45 # inches
HEATMAP_WIDTH  <- 15
bitmap("test.png")
#most highly variable genes ???
pheatmap(counts4[1:500,],
         , width=HEATMAP_WIDTH, height=HEATMAP_HEIGHT
         , scale="row"  # actually scale each row
         , clustering_distance_rows="correlation" # or: euclidean
         , clustering_distance_cols="correlation"
         , clustering_method="complete"
         , display_numbers=FALSE
         , border_color="#00000000" # note: the last two digits here are the alpha channel/ transparency
         , show_rownames=FALSE, show_colnames = FALSE
         , fontsize_row=3.5)
dev.off()

#pca = prcomp(t(counts4), center = T)
pca  <- prcomp(t(assay(rld)[select,]))

#percent variance explained
bitmap('biplot.png')
biplot(pca)
dev.off()

#differential expression 

y = estimateDisp(y, design=design) 
plotBCV(y)
fit <- glmFit(y, design, dispersion=y$trended.dispersion)
lrt1 <- glmLRT(fit, contrast=c(0,1,0,0,0, 0,0))

smoothScatter(x=lrt1$table$logFC, y=-log10(lrt1$table$PValue))