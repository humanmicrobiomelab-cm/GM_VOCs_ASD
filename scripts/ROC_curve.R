

df=t(newdf2[-c(77,78,79)])
df2= as.data.frame(df2)
df2$Group= dat_expr$Group


#CREO TRAINING E TEST
install.packages(caTools)
library(caTools)
split <- sample.split(df2, SplitRatio = 0.7) 

train <- subset(df2, split == "TRUE") 
test <- subset(df2, split == "FALSE") 

model1 <- glm(Group ~ var_1, data=train, family="binomial")

pred1=predict(model1, newdata=test, type="response")
pr1=prediction(pred1, test$Group)
auc1 <- performance(pr1, measure = "auc")
auc1 <- auc1@y.values[[1]]

roc1 <- performance(pr1, measure = "tpr", x.measure = "fpr")
plot(roc1,
     lwd=2,
     col="black")
abline(coef = c(0, 1),
       col = "grey",
       lwd = 2,
       lty=2)
