#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

# Display the current working directory
getwd()

# If necessary, change the path to the directory where the data files are stored. 
setwd("C:/Users/feder/OneDrive/Desktop/WGCNA/")
workingDir = "./"

path_in <- paste(workingDir, "data/", sep = "")
path_out <- paste(workingDir, "script/", sep = "")

# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#Read in the female liver data set
Data = read.table(paste0(path_in, "matrix.txt"), header = T, sep = "\t", check.names = F, quote = "", row.names = 1)

# Take a quick look at what is in the data set:
dim(Data)
names(Data)

####
datExpr0<- apply(Data,2, as.numeric) 
datExpr0 <- log2(datExpr0+1)

row.names(datExpr0) = row.names(Data)

rm(Data)
#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

datExpr0 = as.data.frame(t(datExpr0))
names(datExpr0) 
rownames(datExpr0) 

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
d <- 1 - cor(t(datExpr0), method = "spearman")
mydist <- as.dist(d)
sampleTree = hclust(mydist, method = "average")
#sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================
# If you WANT remove outliers
# Plot a line to show the cut
# abline(h = 240000, col = "red");
# Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = 240000, minSize = 10)
# table(clust)
# clust 1 contains the samples we want to keep.
# keepSamples = (clust==1)
# datExpr = datExpr0[keepSamples, ]

# If you DO NOT WANT remove outliers
datExpr <- datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


rm(datExpr0,Data,gsg,sampleTree,clust,keepSamples,nGenes,nSamples)
#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================
traitData = read.table(paste0(path_in, "ClinicalTraits.txt"), header = T, sep = "\t", check.names = F, quote = "", row.names = 1)
dim(traitData)
names(traitData)
 
# Order datTraits as datExpr
traitData = traitData[rownames(datExpr),]
# 
# # Transform in numeric some columns of datTraits
traitData$Condition <- ifelse((traitData$Condition=="CTRL"),0,1)
# traitData$Sesso <- ifelse((traitData$Sesso=="M"),0,
#                                      ifelse((traitData$Sesso=="F"),1,NA))
#                              
datTraits <- data.frame(sapply(traitData,as.numeric))
rownames(datTraits) <- rownames(traitData)
colnames(datTraits) <- colnames(traitData)
# 
rm(traitData)

#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================

# Re-cluster samples
d <- 1 - cor(t(datExpr), method = "spearman")
mydist <- as.dist(d)
sampleTree = hclust(mydist, method = "average")
#sampleTree = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation:
traitColors = numbers2colors(datTraits, signed = FALSE)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree, traitColors,
                    cex.colorLabels = 0.5, 
                    cex.dendroLabels = 0.5,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

rm(sampleTree,traitColors)
#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================

dirRes <- paste(path_out, "Rdata/", sep = "")
if(!file.exists(dirRes)){
  dir.create(dirRes)
}

filename = paste(dirRes, "dataInput.RData" , sep = "")
save(datExpr, datTraits, file = filename)
