#human cell lines 
# si RNA knockdown experiment

library(pheatmap)
library(edgeR)
library(DESeq2)

# read in files
filename = "~/Downloads/readcounts.txt"
counts1 = read.table(filename, header = T, stringsAsFactors = F) 
rownames(counts1) = counts1$Geneid

# remove metadata
counts2 = counts1[,7:ncol(counts1)]
colnames(counts2) = sub("_1_algn_Aligned", "", matrix(unlist(strsplit(colnames(counts2), "\\.")), byrow =T, ncol=6)[,3])

# filter out genes with low and high counts ( high counts might be mitochondrial genes)
keep = rowMeans(counts2) > 5 & rowMeans(counts2) < 5000 
counts3 = counts2[keep,]

# add count of 1 so that I can take the log later
counts4 = counts3 + 1

# take the log (RPKM)
y <- DGEList(counts=counts4)
fpkm_log_matrix = rpkm (y, gene.length=counts1[keep,"Length"], log = TRUE)

# scatterplot pairwise

# dendrogram 
par(cex=1,font=3)
hc2 <- hclust(stats::dist(t(fpkm_log_matrix), method="minkowski"), "ward.D2")
par(mar=c(10, 4, 4, 10), pty='s')
plot(hc2, lwd=4, lty=1, col="black", col.lab="red", xlab="Samples", main="Samples clustered by log2(FPKM+1) using Ward's clustering & Euclidean distance")

# setup up design matrix for pca plot
groupNamesPerSample = matrix (unlist(strsplit(colnames(counts4), "_")), byrow=T, ncol=4) [,1]
groupNamesPerSample = unlist( lapply (groupNamesPerSample, function(x) sub("S475", "", sub("S449", "", sub("PLC", "", sub("Focus", "", x))))) )
groupNamesPerSample 
conditionFactor  <- factor(groupNamesPerSample)
subjectgroup = matrix (unlist(strsplit(colnames(counts4), "_")), byrow=T, ncol=4) [,1]
subjectgroup = unlist( lapply (subjectgroup, function(x) sub("SC", "", sub("YT", "", sub("T", "", sub("Y", "", x))))) )
subjectgroup
Subject <- factor(subjectgroup)

# pca
cold <- data.frame("conditionFactor"=conditionFactor, "Subject"=Subject) #"sample_name"=colnames(counts4)
design_matrix <- model.matrix( ~Subject + conditionFactor)
ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData=counts4, colData=cold, design=~conditionFactor);
colnames(ddsMat) <- colnames(counts4)
rld <- DESeq2::rlog(ddsMat) # LOG TRANSFORMED.
#rld <- DESeq2::rlog(ddsMat, blind=blind_status) 

pca              <- prcomp(t(assay(rld)))
allPer.vec       <- 100 * summary(pca)$importance[2,]
allPer.text      <- paste(format(allPer.vec,digits=0, nsmall=1, scientific=FALSE), "%", sep='')
allPerLabels.vec <- paste(names(allPer.vec), "\n", allPer.text, sep='')

barplot(allPer.vec, main="PCA: % variance explained by each component\n" , ylim=c(0,100), names.arg=allPerLabels.vec)

d <- data.frame(pca$x, group=conditionFactor, Subject, name=colnames(rld)) # returns PC1, PC2, ... PC8... etc
ppp <- ggplot(d, aes_string(x="PC1", y="PC2", shape="Subject", color="group"))
ppp <- ppp + geom_point(size=6)
plot(ppp)

#concatenate duplicates 
#done manually using counts3 because I don't want to use matrix counts with 1 added to the value
counts5 = data.frame(
            counts3$FocusSC_USR17001380L_HJHCJBBXX_L6
          , counts3$FocusT_USR17001382L_HJHCJBBXX_L7 + counts3$FocusT_USR17001382L_HKF5JBBXX_L1
          , counts3$FocusYT_USR17001383L_hkff7BBXX_L1 + counts3$FocusYT_USR17001383L_HKY5KBBXX_L1
          , counts3$FocusY_USR17001381L_HJHCJBBXX_L6
          
          , counts3$PLCSC_USR17002043L_HJHCJBBXX_L7
          , counts3$PLCT_USR17002044L_HJHCJBBXX_L7
          , counts3$PLCYT_USR17001379L_HJHCJBBXX_L6
          , counts3$PLCY_USR17001377L_HJHCJBBXX_L6
          
          , counts3$S449SC_USR17001384L_HJHCJBBXX_L7
          , counts3$S449T_USR17001386L_HJHCJBBXX_L7 
          , counts3$S449YT_USR17002046L_HJHCJBBXX_L7 + counts3$S449YT_USR17002046L_HKF5JBBXX_L3 
          , counts3$S449Y_USR17001385L_HJHCJBBXX_L7 + counts3$S449Y_USR17001385L_HKF5JBBXX_L3
          
          , counts3$S475SC_USR17001372L_HJHCJBBXX_L6
          , counts3$S475T_USR17001374L_HJHCJBBXX_L6
          , counts3$S475YT_USR17001375L_HJHCJBBXX_L6
          , counts3$S475Y_USR17001373L_HJHCJBBXX_L6 
          )

colnames(counts5) = unlist(lapply (strsplit(sub("counts3.", "" ,colnames(counts5)) , "_"), `[[`, 1))



HEATMAP_HEIGHT <- 45 # inchescol
HEATMAP_WIDTH  <- 15

#most highly variable genes ???
pheatmap(counts5[1:500,],
         , width=HEATMAP_WIDTH, height=HEATMAP_HEIGHT
         , scale="row"  # actually scale each row
         , clustering_distance_rows="correlation" # or: euclidean
         , clustering_distance_cols="correlation"
         , clustering_method="complete"
         , display_numbers=FALSE
         , border_color="#00000000" # note: the last two digits here are the alpha channel/ transparency
         , show_rownames=FALSE, show_colnames = FALSE
         , fontsize_row=3.5)


# setup up design matrix for DE
Condition = unlist( lapply (colnames(counts5), function(x) sub("S475", "", sub("S449", "", sub("PLC", "", sub("Focus", "", x))))) )
Condition  <- factor(Condition)
Subject2 = unlist( lapply (colnames(counts5), function(x) sub("SC", "", sub("YT", "", sub("T", "", sub("Y", "", x))))) )
Subject2 <- factor(Subject2)
cold <- data.frame("Condition"=Condition, "Subject2"=Subject2) #"sample_name"=colnames(counts4)
design_matrix2 <- model.matrix( ~Subject2 + Condition)

#differential expression
y_for_DE <- DGEList(counts=counts5)
y_for_DE <- calcNormFactors(y_for_DE)
y_for_DE = estimateDisp(y_for_DE, design=design_matrix2) 
plotBCV(y_for_DE)
fit <- glmFit(y_for_DE, design=design_matrix2, dispersion=y_for_DE$trended.dispersion)

lrt_T <- glmLRT(fit,coef = "ConditionT")
topTags(lrt_T)
smoothScatter(x=lrt_T$table$logFC, y=-log10(lrt_T$table$PValue))

lrt_Y <- glmLRT(fit,coef = "ConditionY")
topTags(lrt_Y)
smoothScatter(x=lrt_Y$table$logFC, y=-log10(lrt_Y$table$PValue))

lrt_YT <- glmLRT(fit,coef = "ConditionYT")
topTags(lrt_YT)
smoothScatter(x=lrt_YT$table$logFC, y=-log10(lrt_YT$table$PValue))


#concatenate fastq files 
# general gene names
#scatter plot of replicates