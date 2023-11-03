
#Numeric matrix with the rows as individual observations
data_pc=dat_expr[-1]

#PCA
PCA_cov=prcomp(dati_pc, scale. = F)
percentVar_cov <- round(100*PCA_cov$sdev^2/sum(PCA_cov$sdev^2), 1)

data_pca_result <- data.frame(PC1.v = PCA_cov$x[,1], PC2.v = PCA_cov$x[,2],
                     Condition = dat_expr$Group)
sd_ratio_cov <- sqrt(percentVar_cov[2] / percentVar_cov[1])

#score plot of PCA
ggplot(data_pca_result, aes(PC1.v, PC2.v, colour=Condition), label=rownames(PCA_cov$x)) + 
  geom_point(size=4) +
  guides(colour = guide_legend(override.aes=list(shape = 19)),
         fill= guide_legend(title="Condition"))+
  xlab(paste0("PC1 (", percentVar_cov[1], "%)")) +
  ylab(paste0("PC2 (", percentVar_cov[2], "%)")) +
  #xlim(c(-10,10))+
  #ylim(c(-10,10))+
  theme_classic()+
  theme(axis.title = element_text(size=18,colour = "black"),
        axis.text.x = element_text(size=18, colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"),
        axis.line.x = element_line(color="black", size = 1 ),
        axis.line.y = element_line(color="black", size = 1),
        legend.spacing.y = unit(2, "mm"),
        panel.border = element_rect(colour = "black",linewidth = 1,  fill=NA),
        legend.text=element_text(size=18),
        legend.title = element_text(size=18))+
  coord_fixed() +
  scale_color_manual(values = c("royalblue1","orange1"))+
  stat_ellipse(aes(PC1.v, PC2.v), linewidth=1, level=0.95, show.legend = F)+
  scale_x_continuous(breaks = seq(-30, 40, by = 10))

#install mixOmics
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mixOmics")

#Perform PLS-DA

#factor for the discrete outcome
Group= as.factor(dat_expr$Group)

MyResult.pls= plsda(data_pc, Group)

#score plot of PLS-DA
plotIndiv(MyResult.pls, style = "graphics", group = as.factor(dat_expr$Group), legend=F, 
          col.per.group = c( "royalblue1","orange1"), ind.names = T,
          cex=2,
          pch=16,
          legend.title = "Condition", title = "", point.lwd = 2, 
          size.axis = 1.5, size.xlabel = 1.5, size.ylabel = 1.5, size.title = 0.001, 
          ellipse = TRUE,
          ellipse.level = 0.95)
