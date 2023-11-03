
#import metagenomics data
df=read.csv("OTU.csv", header=TRUE, sep=";")
rownames(df)=df$ID
otu=df[3:588]

#Wilcoxon Mann-Whitney test
pvalue_test_otu= sapply(otu[-1], function(i) wilcox.test(i~otu$Group, correction=FALSE)$p.value)

otu_rot=t(otu[-1])
otu_rot=as.data.frame(otu_rot)
otu_rot$pvalue= pvalue_test_otu

sign.otu= subset(otu_rot, pvalue<=0.05)

#matrix of significative OTUs
matrix.otu= t(sign.otu[1:76])

#matrix of significative VOCs
newdf= subset(data_rot, pvalue<=0.05)
matrix.metab= t(newdf[1:76])

# install psych package
install.package("psych")
library(psych)

#Find the correlation matrix and p-value adjusted
cor.1= corr.test(matrix.metab,matrix.otu, method="spearman", adjust="BH")


library(reshape2)
finite.matrix=merge(melt(cor.1$r, value.name="cor"), melt(cor.1$p.adj, value.name="q_value"), by=c("Var1", "Var2"))

#Subset the statistically significant correlations 
sign.matrix= subset(finite.matrix, q_value<=0.05)

write.csv(sign.matrix, "C../corr_sign.csv")


