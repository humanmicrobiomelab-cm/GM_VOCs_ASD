
#Import data
data0 <- read_xlsx("../metabolites.xlsx")
data_expr1= apply(data0[c(3:112)], 2, as.numeric)

#Log2 transformation
data1= as.data.frame(log2(data_expr1+1))

meta_CTRL= as.data.frame(data1[c(1:35),])
meta_ASD= as.data.frame(data1[c(36:76),])

#Calculate mean for each column
lmean_ASD<-colMeans(meta_ASD)
lmean_CTRL<-colMeans(meta_CTRL)

#Find log fold change value
logFC=lmean_ASD-lmean_CTRL

dat_expr=as.data.frame(data1)
dat_expr$Group= data0$Group
dat_expr= dat_expr%>% relocate(Group, .before=`1-butanol`)
rownames(dat_expr)=data0$ID


#Wilcoxon Mann-Whitney test
pvalue_multitest= sapply(dat_expr[-1], function(i) wilcox.test(i~dat_expr$Group, correction=FALSE)$p.value)
adjP= p.adjust(pvalue_multitest, "BH")

results=cbind(logFC, adjP)
results = as.data.frame(results)
results$probename <- rownames(results)

# if log2Foldchange > log2(2) and pvalue adjusted < 0.05, set as "UP" 
results$diffexpressed[results$logFC > 1 & results$adjP <= 0.05] <- "UP"
# if log2Foldchange < log2(2) and pvalue adjusted < 0.05, set as "DOWN"
results$diffexpressed[results$logFC < -1 & results$adjP <= 0.05] <- "DOWN"
