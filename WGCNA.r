library(data.table)
library(WGCNA);
options(stringsAsFactors = FALSE);

# load in all files
expression<-fread("")
expression<-as.data.frame(expression)
pheno<-read.table("")
covs<-read.table("")
# take away excessive information from Counts, leaving just values
 names(expression)[5:length(expression)]=covs[,1]

expression_info<-expression[,1:4]
expression<-expression[,-1:-4]
expression[,length(expression)]<-NULL

pheno<-pheno[pheno[,1] %in% covs[,1],]
covs<-covs[covs[,1] %in% pheno[,1],]

# set column names to sampleIDs
names(expression) = covs$[,1]
#transformed so that the samples are row.names and exons are column names
expression<-as.data.frame(t(expression))
#set column names to exon ID
names(expression)=expression_info[,1]

#set phenotype row names to SampleIDs
rownames(pheno)=pheno[,1]
pheno[,1]<-NULL
#add BMI and age to the phenotype file
pheno$BMI<-covs$BMI
pheno$age<-covs$age

# Cluster samples by Euclidian distance and plot phenotype heatmap below
sampleTree = flashClust(dist(expression), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(pheno, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree, traitColors,
groupLabels = names(pheno),
main = "Dendrogram and phenotype heatmap, Adipose")

# Work out network connectivity and power value to use
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(expression, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



#######################################################################################################
#Automated Blockwise module detection - Change power to suitable value based on plot above





