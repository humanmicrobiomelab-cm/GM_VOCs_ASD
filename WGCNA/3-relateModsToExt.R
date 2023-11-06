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

# Load network data saved in the second part.
filename = paste0(path, "networkConstruction.RData")
load(file = filename);

rm(filename,workingDir)


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
 
 
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 1), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 2, 1, 1));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               #colors = greenWhiteRed(50),
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               cex.lab.x = 1,
               zlim = c(-1,1))
               #main = paste("Module-trait relationships"
#))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Define variable disease_status containing the disease_status column of datTrait
disease_status = as.data.frame(datTraits$Condition);
names(disease_status) = "Condition"


# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
 
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, disease_status, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(disease_status), sep="");
names(GSPvalue) = paste("p.GS.", names(disease_status), sep="");


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
# 
# module = "brown"
# column = match(module, modNames);
# moduleGenes = moduleColors==module;
# 
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(geneModuleMembership[moduleGenes, column],
#                    geneTraitSignificance[moduleGenes, 1],
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = "Gene significance",
#                    #main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,
#                    abline = T, abline.color = "red", abline.lty = 2)
# #=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


#names(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# names(datExpr)[moduleColors==module]


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# Create the starting data frame
geneInfo0 = data.frame(moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for disease_status
modOrder = order(-abs(cor(MEs, disease_status, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}


# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Condition));
geneInfo = geneInfo0[geneOrder, ]

#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================

write.table(geneInfo, "WGCNAInfo.txt", sep = "\t", row.names = T, col.names = NA, quote = F)

