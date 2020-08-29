#################################################################################################
#Step 1: Pre-processing our raw data: here we annotate RNAseq expression samples with the cancer  
#cohorts they were assigned to, as well as common gene IDs
#################################################################################################
#update R if necessary
#install.packages("installr")
#library(installr)
#updateR()

#installing packages and initializing them
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("cowplot")
#install.packages("grid")
#install.packages("e1071")
#install.packages("caret")
#install.packages("gapminder")
#install.packages("pROC")
#install.packages("data.table")
#install.packages("party")
#install.packages("gplots")
library("dplyr")
library("ggplot2")
library("cowplot")
library("grid")
library("e1071")
library("caret")
library("gapminder")
library("pROC")
library("data.table")
library("party")
library("gplots")

#import raw data and label samples with their cancer cohort, as well as common gene IDs
raw_data<-fread("data.csv",header=T)
sample_labels<-fread("labels.csv",header=T)
dt<-inner_join(sample_labels,raw_data,by=NULL,copy=FALSE)
gene_labels<-fread("gene_labels.tsv",header=FALSE)
tempmatrix<-rbind(c("sample_num",0),c("cancer_cohort",0))
new_gene_labels <- rbind(tempmatrix, gene_labels)
colnames(dt)<-new_gene_labels$V1

#saving our annotated data table as an .RData file
save(dt,file = "dt.RData")
#alternatively saving our annotated data table as a pair of two .RData files which can be merged later
dt_part1<-dt[1:round(nrow(dt)/2),]
dt_part2<-dt[(round(nrow(dt)/2)+1):nrow(dt),]
save(dt_part1,file = "dt_part1.RData")
save(dt_part2,file = "dt_part2.RData")




#################################################################################################
#Step 2: Computing RNAseq expression differences between patients from different cancer cohorts  
#for each gene (p-values and expression fold differences)
#################################################################################################

#loading our annotated data table (from line 40) if needed
#load("dt.RData")
#alternatively loading our annotated data from the pair of split tables 
#load("dt_part1.RData")
#load("dt_part2.RData")
#dt<-rbind(dt_part1,dt_part2)

#code to calculate raw p-values, adjusted p-values, and fold differences in expression between breast and prostate cancer patients
cancer_list<-c("BRCA","COAD","KIRC","LUAD","PRAD") 
if (exists("expr_significance")==TRUE) {
  rm(expr_significance)
}
expr_significance<-data.table(cancer1=as.character(), cancer2=as.character(), raw_pval=as.numeric(), expr_fold=as.numeric())

for(x in 1:(length(cancer_list)-1)) {
  for(y in (x+1):length(cancer_list)) {
    for (i in 32:ncol(dt)) {
      cancer1<-cancer_list[x]
      cancer2<-cancer_list[y]
      raw_pval<-t.test(dt[dt$cancer_cohort==cancer_list[x],..i],dt[dt$cancer_cohort==cancer_list[y],..i],var.equal = FALSE)$p.value
      expr_fold<-mean(as.numeric(unlist(dt[dt$cancer_cohort==cancer_list[x],..i])))-mean(as.numeric(unlist(dt[dt$cancer_cohort==cancer_list[y],..i])))
      tempDT<-data.table(cancer1,cancer2,raw_pval,expr_fold)
      expr_significance<-rbind(expr_significance,tempDT)
    }
  }
}
expr_significance<-replace(expr_significance,is.na(expr_significance),1)
expr_significance$adj_pval<-p.adjust(expr_significance$raw_pval,method = "bonferroni")

#saving our data table containing raw/adjusted p-values and expression fold difference as an RData file
save(expr_significance,file="expr_significance.RData")



#################################################################################################
#Step 3: Building a Naive Bayes Classifier based on genes that were differentially expressed
#################################################################################################

#loading our annotated data table (from line 75) if needed
#load("expr_significance.RData")

#filtering our dataset to only include genes that differentially expressed in the pair-wise
#comparisons via the following two criteria:
#1) the p-value adjusted for multiple test corrections must be below 0.05
#2) the fold difference in gene expression must be greater than 3
#we identified that 
filtered_expr_significance<-expr_significance %>% filter(abs(expr_fold)>log2(3),adj_pval<0.05)
cancer_pairs_diff_expr_counts<-as.data.table(filtered_expr_significance %>% group_by(cancer1,cancer2) %>% tally())
save(cancer_pairs_diff_expr_counts, file = "cancer_pairs_diff_expr_counts")
load("cancer_pairs_diff_expr_counts")

#heatmap of the pairwise comparisons of the number of differentially expressed genes between cancer cohorts
ggplot(cancer_pairs_diff_expr_counts, aes(cancer1, cancer2, fill= n)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red") + theme_minimal() +  
  labs(title="Differential Gene Expression Among Cancer Cohorts", fill="Number of Genes", x=element_blank() , y=element_blank()) + 
  theme(legend.position="right", plot.title = element_text(hjust = 0.5))

#identifying the index of differentially expressed genes between 2 cancers with the most similar gene expression profiles
brca_luad<-expr_significance %>% filter(cancer1=="BRCA",cancer2=="LUAD")
criteria_adj_pval<-which(brca_luad$adj_pval < 0.05)
criteria_expr_fold<-c(which(abs(brca_luad$expr_fold) > log2(3)))
target_gene_columns<-c(intersect(criteria_adj_pval,criteria_expr_fold)+31)

#assembling our response vector and feature matrix
target_matrix_columns<-c(2,target_gene_columns)
feature_matrix<-dt %>% select(target_matrix_columns) %>% filter(cancer_cohort == "BRCA" | cancer_cohort == "LUAD")
feature_matrix$cancer_cohort <- as.factor(feature_matrix$cancer_cohort)
save(feature_matrix,file="feature_matrix.RData")

#running a naive bayes classifier with 5-fold cross validation
n_folds <- 5
folds_i <- sample(rep(1:n_folds, length.out = nrow(feature_matrix)))
table(folds_i)

for(k in 1:n_folds){
  #select which group is our test data and define our training/testing data
  test_i <- which(folds_i == k)
  train_data <- feature_matrix[-test_i, ]
  test_data <- feature_matrix[test_i, ]
  
  #train classifier
  classifier_nb <- naiveBayes(train_data[,-1], train_data$cancer_cohort)
  
  #use classifier to predict classes of test data
  nb_pred_probabilities <- predict(classifier_nb, newdata = test_data, type = "raw")
  
  nb_pred_classes <- predict(classifier_nb, newdata = test_data, type = "class")
  
  #assess the accuracy of the classifier on the test data using a confusion matrix
  print(confusionMatrix(nb_pred_classes,test_data$cancer_cohort))
  
  #assess the accuracy of the classifier on the test data using a ROC (function defined above)
  ROC <- roc(predictor=nb_pred_probabilities[,1],
             response=test_data$cancer_cohort)
  print(ROC$auc)
  plot(ROC,main="ROC Curve for Naive Bayes Classifier")
}



#################################################################################################
#Step 4: Building a Decision Tree based on genes that were differentially expressed
#################################################################################################

#running a decision tree with 5-fold cross validation

for(k in 1:n_folds){
  #select which group is our test data and define our training/testing data
  test_i <- which(folds_i == k)
  train_data <- feature_matrix[-test_i, ]
  test_data <- feature_matrix[test_i, ]
  
  #train classifier
  classifier_ctree <- ctree(cancer_cohort ~ ., data=train_data)
  
  #plot our decision tree model
  print(classifier_ctree)
  plot(classifier_ctree) 
  plot(classifier_ctree,main="Classification tree for BRCA/LUAD based on gene expression")
  
  #run our decision tree on test data to get predicted classes/probabilities
  predicted_class<-predict(classifier_ctree, test_data, type="response")    
  predicted_prob<-sapply(predict(classifier_ctree, test_data,type="prob"),'[[',1) 
  
  #print confusion matrix
  print(confusionMatrix(predicted_class,test_data$cancer_cohort))
  
  #plot ROC curve
  ROC <- roc(predictor=predicted_prob,
             response=test_data$cancer_cohort)
  print(ROC$auc)
  plot(ROC,main="ROC Curve for Decision Tree Model")
}        




#################################################################################################
#Step 5: Determining how long each model takes to get trained and how long each model takes to 
#predict classes of test data 
#################################################################################################
#naive bayes times
start.time <- Sys.time()
classifier_nb <- naiveBayes(train_data[,-1], train_data$cancer_cohort)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste("Time to train Naive Bayes Classifier:",time.taken,"seconds"))
  
start.time <- Sys.time()
nb_pred_classes <- predict(classifier_nb, newdata = test_data, type = "class")
end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste("Time to predict classes of test data Using Naive Bayes Classifier:",time.taken,"seconds"))
  
#decision tree times
start.time <- Sys.time()
classifier_ctree <- ctree(cancer_cohort ~ ., data=train_data)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste("Time to train Decision Tree:",time.taken,"seconds"))
  
start.time <- Sys.time()
classifier_ctree <- ctree(cancer_cohort ~ ., data=train_data)
predicted_class<-predict(classifier_ctree, test_data, type="response")
end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste("Time to predict classes of test data Using Decision Tree Model:",time.taken,"seconds"))
