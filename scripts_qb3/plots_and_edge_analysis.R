library(pheatmap)
library(edgeR)
library(DESeq2)

filename = "~/Downloads/readcounts.txt"
counts1 = read.table(filename, header = T, stringsAsFactors = F) 
rownames(counts1) = counts1$Geneid
counts2 = counts1[,7:ncol(counts1)]
keep = rowMeans(counts2) > 5 & rowMeans(counts2) < 5000 
colnames(counts2) = sub("_1_algn_Aligned", "", matrix(unlist(strsplit(colnames(counts2), "\\.")), byrow =T, ncol=6)[,3])
counts3 = counts2[keep,]
counts3 = counts3 + 1

y <- DGEList(counts=counts3)
fpkm_log_matrix = rpkm (y, gene.length=counts1[keep,"Length"], log = TRUE)

#scatterplot pairwise

#dendrogram 
par(cex=1,font=3)
hc2 <- hclust(stats::dist(t(fpkm_log_matrix), method="minkowski"), "ward.D2")
par(mar=c(10, 4, 4, 10), pty='s')
plot(hc2, lwd=4, lty=1, col="black", col.lab="red", xlab="Samples", main="Samples clustered by log2(FPKM+1) using Ward's clustering & Euclidean distance")

groupNamesPerSample = matrix (unlist(strsplit(colnames(mmm), "_")), byrow=T, ncol=4) [,1]
groupNamesPerSample = unlist( lapply (groupNamesPerSample, function(x) sub("S475", "", sub("S449", "", sub("PLC", "", sub("Focus", "", x))))) )
groupNamesPerSample 
conditionFactor  <- factor(groupNamesPerSample)
subjectgroup = matrix (unlist(strsplit(colnames(mmm), "_")), byrow=T, ncol=4) [,1]
subjectgroup = unlist( lapply (subjectgroup, function(x) sub("SC", "", sub("YT", "", sub("T", "", sub("Y", "", x))))) )
subjectgroup
Subject <- factor(subjectgroup)

#pca
mmm = counts3
cold <- data.frame("conditionFactor"=conditionFactor, "Subject"=Subject) #"sample_name"=colnames(mmm)
design_matrix <- model.matrix( ~Subject + conditionFactor)
ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData=mmm, colData=cold, design=~conditionFactor);
colnames(ddsMat) <- colnames(mmm)
rld <- DESeq2::rlog(ddsMat) # LOG TRANSFORMED.
#rld <- DESeq2::rlog(ddsMat, blind=blind_status) 

pca              <- prcomp(t(assay(rld)))
allPer.vec       <- 100 * summary(pca)$importance[2,]
allPer.text      <- paste(format(allPer.vec,digits=0, nsmall=1, scientific=FALSE), "%", sep='')
allPerLabels.vec <- paste(names(allPer.vec), "\n", allPer.text, sep='')

barplot(allPer.vec, main=paste0("PCA: % variance explained by each component\nCalculated with the ", length(select), " highest-variance genes", additional_main), ylim=c(0,100), names.arg=allPerLabels.vec)

d <- data.frame(pca$x, group=conditionFactor, Subject, name=colnames(rld)) # returns PC1, PC2, ... PC8... etc
ppp <- ggplot(d, aes_string(x="PC1", y="PC2", shape="Subject", color="group"))
ppp <- ppp + geom_point(size=6)
plot(ppp)

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

#differential expression 
y <- calcNormFactors(y)
y = estimateDisp(y, design=design_matrix) 
plotBCV(y)
fit <- glmFit(y, design=design_matrix, dispersion=y$trended.dispersion)
lrt <- glmLRT(fit,coef = "conditionFactorT")
topTags(lrt)
smoothScatter(x=lrt$table$logFC, y=-log10(lrt$table$PValue))

#concatenate fastq files 
# general gene names
#scatter plot of replicates