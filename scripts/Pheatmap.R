
#Heatmap 
# install pheatmap package
install.package("pheatmap")

data_rot=as.data.frame(t(dat_expr[-1]))

samples=data0[2]
annotation<-(data.frame(SampleType=samples))
rownames(annotation) <- colnames(data_rot[c(1:76)])
names(annoCol) <- unique(samples$Group)
annoCol <- list(Group=c(ASD="royalblue1",CTRL="orange1"))

names(vect) = unique(samples$Group)

#Heatmap with all VOCs
pheatmap(data_rot[c(1:76)], scale = "row",
         border_color = NA,
         cluster_cols = T,
         cluster_rows = T,
         #cex=1,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         annotation_col = annotation,
         annotation_colors = annoCol,
         color = colorRampPalette(colors = c("blue", "blue3", "black", "yellow3", "yellow"))(100),
         show_rownames = F,
         show_colnames = F
         #cutree_cols = 2,
         #cutree_rows = 2,
)




#with anamnestic variable
metadata1=metadata
colnames(metadata1)[2]<-"Gender"

#Create double annotation
samples2=metadata1[c(4,2)]
samples2$Gender[samples2$Gender==0]="Male"
samples2$Gender[samples2$Gender==1] <- "Female"
annotation2<-(data.frame(SampleType=samples2))
rownames(annotation2) <- colnames(data_rot[c(1:76)])
colnames(annotation2)= c("Group","Gender")

annoCol2 <- list(
  Group=c(ASD="royalblue1",CTRL="orange1"),
  Gender=c("Male"="lightblue", "Female"="lightpink")
)
names(annoCol2) <- c("Group","Gender")



#Consider VOCs with p-value adjusted <=0.05
newdf2=subset(data_rot, pvalue_adj<=0.05)

#Heatmap of only statistically significant VOCs between ASD and CTRL with 2 annotation
pheatmap(newdf2[c(1:76)], scale = "row",
         border_color = NA,
         cluster_cols = T,
         cluster_rows = T,
         clustering_method = "ward.D2",
         annotation_col = annotation2,
         annotation_colors = annoCol2,
         color = colorRampPalette(colors = c("blue", "blue3", "black", "yellow3", "yellow"))(100),
         show_rownames = T,
         show_colnames = T,
         fontsize_row = 15,
         fontsize = 12,
         fontsize_col = 10
         #cutree_cols = 2,
         #cutree_rows = 2,
)


