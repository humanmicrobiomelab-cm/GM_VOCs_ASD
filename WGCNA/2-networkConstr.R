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

path <- paste(workingDir, "script/Rdata/", sep = "")

# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Load the data saved in the first part
filename = paste0(path, "dataInput.RData")
load(file = filename)

rm(filename,workingDir)
#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers = c(1:10, seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
par(mfrow = c(1,2), mar = c(5,5,2,2));
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = 1, 
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       saveTOMFileBase = paste(path, "TOM", sep = ""),
                       verbose = 3, maxBlockSize = 30000)

# Use the function "recutBlockwiseTreeswould" to change some of the tree cut, module membership, and module merging criteria without having to recompute the network and the clustering dendrogram

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
#table(net$colors) #to know number and size of modules
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms;

#BarPlot of of number of genes for each module
t <- data.frame(table(moduleColors))
t <- t[order(t$Freq, decreasing = T),]
par(mar=c(3,3,3,1))
lab <- paste0(as.character(t$moduleColors),"(",as.character(t$Freq),")")
barplot(t$Freq, col = as.character(t$moduleColors), 
        ylab = "# Genes", names.arg = lab, las = 2, cex.names = 1.1)

filename = paste0(path, "networkConstruction.RData")
save(MEs, moduleLabels, moduleColors, geneTree, file = filename)

