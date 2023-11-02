
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mixOmics")

#Numeric matrix with the rows as individual observations
data_pc=dat_expr[-1]

#PCA
MyResult.pca <- pca(data_pc, scale = F) 

#score plot of PCA
plotIndiv(MyResult.pca, style = "graphics", group = as.factor(dat_expr$Group), legend=F, 
          col.per.group = c( "royalblue1","orange1"), ind.names = F,
          X.label = "PC1 (17.8%)",
          Y.label = "PC2 (10.6%)",
          cex=1.5,
          pch=19,
          legend.title = "Condition", title = "", point.lwd = 2, 
          size.axis = 1.5, size.xlabel = 1.5, size.ylabel = 1.5, size.title = 0.001, 
          ellipse = TRUE)

#Perform PLS-DA

#factor for the discrete outcome
Group= as.factor(dat_expr$Group)

MyResult.pls= plsda(data_pc, Group, scale = T)

#score plot of PLS-DA
plotIndiv(MyResult.pls, style = "graphics", group = as.factor(dat_expr$Group), legend=F, 
          col.per.group = c( "royalblue1","orange1"), ind.names = T,
          X.label = "C1 (12%)",
          Y.label = "C2 (8%)",
          cex=2,
          pch=16,
          legend.title = "Condition", title = "", point.lwd = 2, 
          size.axis = 1.5, size.xlabel = 1.5, size.ylabel = 1.5, size.title = 0.001, 
          ellipse = TRUE,
          ellipse.level = 0.95,
          xlim=c(-7,8.5))
