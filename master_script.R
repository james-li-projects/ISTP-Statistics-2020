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
#install.packages("tidyverse")
#install.packages("Hmisc")
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
library("tidyverse")
library("Hmisc")

#import raw data and label samples with their cancer cohort, as well as common gene IDs
raw_data<-fread("data.csv",header=T)
sample_labels<-fread("labels.csv",header=T)
dt<-inner_join(sample_labels,raw_data,by=NULL,copy=FALSE)
gene_labels<-fread("gene_labels.tsv",header=FALSE)
tempmatrix<-rbind(c("sample_num",0),c("cancer_cohort",0))
new_gene_labels <- rbind(tempmatrix, gene_labels)
colnames(dt)<-new_gene_labels$V1

#saving our annotated data table as an .RData file
#save(dt,file = "dt.RData")
#alternatively saving our annotated data table as a pair of two .RData files which can be merged later in order to meet GitHub max individual file size
dt_part1<-dt[1:round(nrow(dt)/2),]
dt_part2<-dt[(round(nrow(dt)/2)+1):nrow(dt),]
save(dt_part1,file = "dt_part1.RData")
save(dt_part2,file = "dt_part2.RData")




#################################################################################################
#Step 2: Computing RNAseq expression differences between patients from different cancer cohorts  
#for each gene (p-values and expression fold differences)
#################################################################################################

#loading our annotated data table (from line 48) if needed
#load("dt.RData")
#alternatively loading our annotated data from the pair of split tables (from lines 52 and 53) 
load("dt_part1.RData")
load("dt_part2.RData")
dt<-rbind(dt_part1,dt_part2)

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

#loading our annotated data table (from line 93) if needed
load("expr_significance.RData")

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

#plotting confusion matrix function: courtesy of Cybernetic
draw_confusion_matrix <- function(cm) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, 'Class1', cex=1.2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, 'Class2', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, 'Class1', cex=1.2, srt=90)
  text(140, 335, 'Class2', cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  


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
  nb_cm<-confusionMatrix(nb_pred_classes,test_data$cancer_cohort)
  draw_confusion_matrix(nb_cm)
  
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
  ctree_cm<-confusionMatrix(predicted_class,test_data$cancer_cohort)
  draw_confusion_matrix(ctree_cm)
  
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
#sampling some run time data
if (exists("run_time_data")==TRUE) {
  rm(run_time_data)
}
run_time_data<-data.table(Train_Naive_Bayes=as.numeric(), Test_Naive_Bayes=as.numeric(), Train_Decision_Tree=as.numeric(), Test_Decision_Tree=as.numeric())

for (i in 1:20) {
  #naive bayes times
  start.time <- Sys.time()
  classifier_nb <- naiveBayes(train_data[,-1], train_data$cancer_cohort)
  end.time <- Sys.time()
  time.taken_train_NB <- end.time - start.time
  Train_Naive_Bayes <- as.numeric(time.taken_train_NB)
  print(paste("Time to train Naive Bayes Classifier:",Train_Naive_Bayes,"seconds"))
  
  start.time <- Sys.time()
  nb_pred_classes <- predict(classifier_nb, newdata = test_data, type = "class")
  end.time <- Sys.time()
  time.taken_test_NB <- end.time - start.time
  Test_Naive_Bayes <- as.numeric(time.taken_test_NB)
  print(paste("Time to predict classes of test data Using Naive Bayes Classifier:",Test_Naive_Bayes,"seconds"))
  
  #decision tree times
  start.time <- Sys.time()
  classifier_ctree <- ctree(cancer_cohort ~ ., data=train_data)
  end.time <- Sys.time()
  time.taken_train_ctree <- end.time - start.time
  Train_Decision_Tree <- as.numeric(time.taken_train_ctree)
  print(paste("Time to train Decision Tree:",Train_Decision_Tree,"seconds"))
  
  start.time <- Sys.time()
  classifier_ctree <- ctree(cancer_cohort ~ ., data=train_data)
  predicted_class<-predict(classifier_ctree, test_data, type="response")
  end.time <- Sys.time()
  time.taken_test_ctree <- end.time - start.time
  Test_Decision_Tree <- as.numeric(time.taken_test_ctree)
  print(paste("Time to predict classes of test data Using Decision Tree Model:",Test_Decision_Tree,"seconds"))
  
  temp_run_time_data<-data.table(Train_Naive_Bayes,Test_Naive_Bayes,Train_Decision_Tree,Test_Decision_Tree)
  run_time_data<-rbind(run_time_data,temp_run_time_data)
}

#####plotting run times#####
train_times<-run_time_data %>% select(Train_Naive_Bayes,Train_Decision_Tree)
test_times<-run_time_data %>% select(Test_Naive_Bayes,Test_Decision_Tree)
require(reshape2) 

#plot model train times
train_times_melt<-melt(train_times)
colnames(train_times_melt)<-c("Stage_in_Model","Time_in_seconds")
ggplot(data = train_times_melt, aes(x=Stage_in_Model, y=Time_in_seconds)) + geom_boxplot(aes(fill=Stage_in_Model)) + theme(legend.position="none") + ggtitle("Boxplots of Time in Seconds it takes to train Naive Bayes vs Decision Tree Models (n=20)")

#plot model test times
test_times_melt<-melt(test_times)
colnames(test_times_melt)<-c("Stage_in_Model","Time_in_seconds")
ggplot(data = test_times_melt, aes(x=Stage_in_Model, y=Time_in_seconds)) + geom_boxplot(aes(fill=Stage_in_Model)) + theme(legend.position="none") + ggtitle("Boxplots of Time in Seconds it takes to test Naive Bayes vs Decision Tree Models (n=20)")




#################################################################################################
#Step 6: Visualizing trends in the gene expression data
#################################################################################################

# PLOTS OF ORIGINAL DATA (BRCA VS. LUAD)
# Plot expression fold difference against p value.
ggplot(brca_luad, aes(x=adj_pval, y=expr_fold)) + 
  geom_point(shape=1) +
  geom_hline(yintercept=log(3,base=2)) +
  geom_hline(yintercept=-log(3,base=2)) +
  geom_vline(xintercept=.05) + 
  scale_x_log10() +
  labs(x="Adjusted p-value", y="log2[Expression Fold Difference]") +
  theme()

# Prepare violin plots of expression of selected genes among patients.
brca_dt<-dt[dt$cancer_cohort=="BRCA",]
luad_dt<-dt[dt$cancer_cohort=="LUAD",]
combined_brcaluad_dt<-rbind(brca_dt,luad_dt)

# The first step is to convert data in combined_brcaluad_dt to a more useful format.

# Create vector (genes) of genes of interest.
# Create vector (subset_names) to store names of data tables for each gene.
# Create data.table (data_to_plot) which will contain all the data to plot.
# Individual data subsets will be sequentially added to data_to_plot.
genes <- c("KRAS","EGFR","MET","RET","ROS1")
subset_names <- rep("",length(genes))
data_to_plot = data.table(cancer_cohort = character(), 
                          expression_level = numeric(), gene = character())

# Summary: For each gene of interest in genes, 
# we add relevant data from combined_brcaluad_dt to data_to_plot.
# Details: 
#   indices is a vector containing the columns of interest,
#       specifically cancer_cohort and (expression level of) genes[i].
#   subset_names[i] is updated to contain the correct name of the 
#       current subset we are working with.
#   A data.table is created with name subset_names[i] and three columns:
#       cancer_cohort, (expression level of) genes[i], and a new column (gene).
#       The new column contains the name of the current gene (genes[i]).
#       cbind() combines the first two columns with the new third column.
#   The second column is renamed from genes[i] to "expression_level"
#       which is necessary for consistency.
#       eval(str2lang()) evaluates subset_names[i] as the name of a data.table.
#   The data.table subset_names[i] is added to data_to_plot with rbind().
#       Again, eval(str2lang()) is used.
for(i in 1:length(genes)){
  indices <- c(2, grep(paste("\\b", genes[i], "\\b", sep=""), 
                       colnames(combined_brcaluad_dt)))
  subset_names[i] <- paste("subset_", genes[i], sep="")
  assign(subset_names[i], 
         cbind(combined_brcaluad_dt[,..indices],
               data.table(gene = rep(genes[i],nrow(combined_brcaluad_dt)))))
  setnames(eval(str2lang((subset_names[i]))), genes[i], "expression_level")
  data_to_plot <- rbind(data_to_plot, eval(str2lang((subset_names[i]))))
}

# Now that the data in data_to_plot is in a useful format to make violin plots,
# i.e. three columns with gene name and expression level as separate variables,
# we can proceed with creating violin plots.
windows(11,5)
plot_violin <- ggplot(data_to_plot, aes(x=gene, y=expression_level, fill=cancer_cohort)) + 
  geom_violin(trim=TRUE, position=position_dodge(1)) + 
  labs(x="Gene of Interest", y="Normalized Expression (a.u.)",
       title="Distributions in Gene Expression among Patients",
       fill="Cancer Type") +
  geom_boxplot(width=0.1, color="black", position=position_dodge(1)) +
  theme_classic() +
  theme(legend.position="right", plot.title = element_text(hjust = 0.5))
plot_violin

# Repeat the same process as above for a second set of genes:
genes <- c("ATM","BARD1","BRIP1","CDH1","CHEK2","NBN","NF1","PALB2","PTEN", 
           "RAD51C","RAD51D","STK11","TP53")
subset_names <- rep("",length(genes))
data_to_plot = data.table(cancer_cohort = character(), 
                          expression_level = numeric(), gene = character())

for(i in 1:length(genes)){
  indices <- c(2, grep(paste("\\b", genes[i], "\\b", sep=""), 
                       colnames(combined_brcaluad_dt)))
  subset_names[i] <- paste("subset_", genes[i], sep="")
  assign(subset_names[i], 
         cbind(combined_brcaluad_dt[,..indices],
               data.table(gene = rep(genes[i],nrow(combined_brcaluad_dt)))))
  setnames(eval(str2lang((subset_names[i]))), genes[i], "expression_level")
  data_to_plot <- rbind(data_to_plot, eval(str2lang((subset_names[i]))))
}

windows(11,5)
plot_violin <- ggplot(data_to_plot, aes(x=gene, y=expression_level, fill=cancer_cohort)) + 
  geom_violin(trim=TRUE, position=position_dodge(1)) + 
  labs(x="Gene of Interest", y="Normalized Expression (a.u.)",
       title="Distributions in Gene Expression among Patients",
       fill="Cancer Type") +
  geom_boxplot(width=0.1, color="black", position=position_dodge(1)) +
  theme_classic() +
  theme(legend.position="right", plot.title = element_text(hjust = 0.5))
plot_violin
