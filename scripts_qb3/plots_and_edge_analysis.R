# Goal: 
#compare affected to unaffected

library(pheatmap)
library(edgeR)
#scatterplot
#pca highest variance genes
#dendrogram (why does it algn with tsne clusters????)
#DE

filename = "Downloads/readcounts.txt"
counts1 = read.table(filename, header = T) 
counts2 = counts[,7:ncol(counts1)]
#keep = which(apply(counts2, 1, function(x) sum(which(x>0))==ncol(counts2)))
keep = rowMeans(counts2) > 5 & rowMeans(counts2) < 5000 

counts3 = counts2[keep,]
counts3 = counts3 + 1
y <- DGEList(counts=counts3)
y <- calcNormFactors(y)

counts4 = rpkm (y, gene.length=counts[keep,"Length"], log = TRUE)

#scatterplot pairwise

#dendrogram 


bitmap("test.png")
#most highly variable genes ???
pheatmap(counts4[1:500,], show_rownames=FALSE, show_colnames = FALSE)
dev.off()

pca = prcomp(t(counts4), center = T)
#percent variance explained
bitmap('biplot.png')
biplot(pca)
dev.off()

#differential expression 