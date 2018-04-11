#### STATISTICS WITH R - MSC BIOINFORMATICS AND SYSTEMS BIOLOGY #####
#### Alberto Gil - student no 2595259 ####
# Load the required libraries
library(broom)
library(knitr)
library(class) # for knn() function
library(randomForest) # for randomForest() package
library(e1071) # for naiveBayes() function
library(caret) # for accuracies estimation
library(ggplot2)
library(gplots) # for doing the heatmaps

######################### 1 ################################

# INPUT p: required number of partitions 
#       f: list of one or more factors to stratify over
# OUTPUT integer vectorcontaining partition labels (1 .. p).
partition <- function(p, f) {
  if (is.atomic(f)) {
    lenf <- length(f)
  } else {
    lenf < length(f[[1]])
  }
  part <- vector(mode="integer", length=lenf)
  tapply(X=1:length(part), INDEX=f, FUN=function(i) {
    part[i] <<- as.integer((sample(1:length(i)) + sample(1:p,1))%%p + 1)
  }
  )
  return(part)
}

# Crossvalidates a selected model with k-fold stratified cross
# validation, and reports the model performance
# INPUT: model: 'knn', 'RF' or'naive'
#        train: training set dataframe
#        labels: list of labels (target variable) of training set
#        npart: number of partitions of dataframe
#        modelsperformance: Dataframe for storing performance of previous models (values
#                         will be added in a new column)
#       (optional) k: k parameter for k-NN model
#       (optional) print_mean: Boolean; set to yes if mean and sd are desired to be in the final dataframe
# OUTPUT: modelsperformance (dataframe containing performance of the models )
crossvalidate_model <- function(model,train,trainlabels,npart,
                                models_performance,k=1,print_mean=FALSE){
  set.seed(2) # Set a seed for generating random numbers
  part <- partition(npart,trainlabels) # Create stratified partitions
  if (is.null(models_performance)){
    # Create a new dataframe for storing model performance
    models_performance <- data.frame(partition = min(part):max(part))
  }
  for (i in min(part):max(part)){
    if (model == 'knn'){
      set.seed(2)
      predictions <- knn(train[part != i,],train[part == i,],
                         trainlabels[part != i],k)
    }
    if (model == 'RF'){
      set.seed(2)
      rf <- randomForest(x=train[part != i,], y=trainlabels[part != i],
                         proximity=TRUE,importance=TRUE)
      predictions <- predict(rf, train[part == i,])
    }
    if (model == 'naive'){
      nb <- naiveBayes(x=train[part != i,], y=trainlabels[part != i])
      predictions <- predict(nb, train[part == i,])
    }
    # Calculate the error 
    error <- mean(predictions != trainlabels[part == i])
    models_performance[i,model] <- error
  }
  # If desired by the user, print the mean performance in the df 
  if(print_mean){
    g <- as.data.frame(lapply(models_performance,mean))
    g[1,1] <- 'mean'
    h <- as.data.frame(lapply(models_performance,sd))
    h[1,1] <- 'sd'
    models_performance[length(rownames(models_performance))+1,] <- g
    models_performance[length(rownames(models_performance))+1,] <- h
  }
  return(models_performance)
}

# Load the datasets
ecolidataset <- read.table('data/ecoli.data',header=TRUE,sep="",
                           comment.char = "#", row.names=1)
yeastdataset <- read.table('data/yeast.data',header=TRUE,sep="",
                           comment.char = "#")
# Retrieve the labels from the dataset
ecoli.labels <- as.list(subset(ecolidataset,select = c(location)))$location
ecoli.train <- subset(ecolidataset, select = -c(location))
# Retrive the performance from different models (E. coli dataset)
ecoli_performance <- crossvalidate_model('knn',ecoli.train,ecoli.labels,4,NULL,7)
ecoli_performance <- crossvalidate_model('RF',ecoli.train,ecoli.labels,4,
                                         ecoli_performance)
ecoli_performance <- crossvalidate_model('naive',ecoli.train,ecoli.labels,4,
                                         ecoli_performance,7,FALSE)
# Perform paired t-tests
t.test(as.list(ecoli_performance['knn'])$knn,
       as.list(ecoli_performance['RF'])$RF, paired=TRUE,conf.level=0.95)
t.test(as.list(ecoli_performance['knn'])$knn,
       as.list(ecoli_performance['naive'])$'naive', paired=TRUE,conf.level=0.95)
t.test(as.list(ecoli_performance['RF'])$RF,
       as.list(ecoli_performance['naive'])$'naive', paired=TRUE,conf.level=0.95)

# Explore the correlation between variables
ecoli.train <- data.frame(apply(ecoli.train, 2, 
                                function(x) as.numeric(as.character(x))))
pears_cor <- cor(ecoli.train,method="pearson")
hmcol <- colorRampPalette(c("red", "white","blue"))
heatmap.2(pears_cor,col = hmcol,main="Pearson correlation")
spear_cor <- cor(ecoli.train,method="spearman")
heatmap.2(spear_cor,col = hmcol,main="Spearman correlation")

# Create dataframes for the yeast dataset
yeastdataset <- yeastdataset[unique(yeastdataset[,1]),]
yeast.labels <- as.list(subset(yeastdataset,select = c(location)))$location
yeast.train <- subset(yeastdataset, select = -c(SpAccNr,location))

# Evaluate performance on models built from yeast dataset
yeast_performance <- crossvalidate_model('knn',yeast.train,yeast.labels,10,NULL,21)
yeast_performance <- crossvalidate_model('RF',yeast.train,yeast.labels,10,
                                         yeast_performance)
yeast_performance <- crossvalidate_model('naive',ecoli.train,yeast.labels,10,
                                         yeast_performance,7,TRUE)

# Perform a paired t-test of the performances of yeast models
t.test(as.list(yeast_performance['knn'])$knn,
       as.list(yeast_performance['RF'])$RF, paired=TRUE)
t.test(as.list(yeast_performance['knn'])$knn,
       as.list(yeast_performance['naive'])$'naive', paired=TRUE)
t.test(as.list(yeast_performance['RF'])$RF,
       as.list(yeast_performance['naive'])$'naive', paired=TRUE)

# Explore the correlation between variables
yeast.train <- data.frame(apply(yeast.train, 2, 
                                function(x) as.numeric(as.character(x))))
pears_cor <- cor(yeast.train,method="pearson")
hmcol <- colorRampPalette(c("red", "white","blue"))
heatmap.2(pears_cor,col = hmcol,main="Pearson correlation")
spear_cor <- cor(yeast.train,method="spearman")
heatmap.2(spear_cor,col = hmcol,main="Spearman correlation")


############################### 2. ###############################
# Read the datasets files
soildataset <- read.table('heavy_metal_data/Log2GR31_Q_Ave_DM.txt',header=TRUE,sep="",
                           comment.char = "#", row.names=1)
soildataset.samplesinfo <- read.table('heavy_metal_data/Samples.txt',header=TRUE,sep="",
                                      comment.char = "#",fill = TRUE)
# Retrieve only the information from the metal on the sample
soildataset.samplesinfo <- soildataset.samplesinfo$metal 
soildataset <- soildataset[-1,] # Delete first row

# Convert the dataframe to numeric (will avoid further problems in randomForest() func.)
df2 <- data.frame(apply(soildataset, 2, function(x) as.numeric(as.character(x))))
rownames(df2) <- rownames(soildataset)


# Creates a hierarchical tree from a distance matrix
# INPUT: dmatrix:    distance matrix
#       labels_func: list of labels
#       title:       title of the tree
createhierarchicaltree <- function(dmatrix,labels_func,title=NULL){
  clusters <- hclust(dmatrix,method="complete")
  dend <- as.dendrogram(clusters)
  mycol <- rainbow(length(unique(labels_func)))
  samplenum <- as.numeric(labels_func[order.dendrogram(dend)])
  leafcolor <- function(node) {
    if (is.leaf(node)) {
      # The <<- operator must be used instead of <- so that this function searches outside its own environment for the globally defined i.
      i <<- i + 1
      attr(node, "nodePar") <- list(pch=20, cex=1.5, col=mycol[samplenum[i]] )
    }
    return(node)
  }
  i <- 0
  dend <- dendrapply(dend, leafcolor)
  plot(dend, leaflab="none", axes=FALSE,main=title)
  legend("topright",levels(labels_func), pch=19, col=mycol,cex=0.5)
}

# Create hierarchical trees using Pearson correlation as distance matrix
createhierarchicaltree(as.dist(1-cor(df2,method="pearson")),
                       soildataset.samplesinfo,
                       'Distance function: Pearon correlation')
# Create hierarchical trees using Spearman correlation as distance matrix
createhierarchicaltree(as.dist(1-cor(df2,method="spearman")),
                       soildataset.samplesinfo,
                       'Distance function: Spearman correlation')

# B. Separate UNLABELLED DATA from labelled data
# Transpose the original df for having samples in rows, and genes (features) in columns
df2 <- t(df2)
# Retrieve rows with unlabelled data and create a new df with this data 
unlabeled_data_rows <- as.logical((soildataset.samplesinfo == 'Noy') +
                                    (soildataset.samplesinfo == 'PLO' ))
unlabeled_data <- df2[unlabeled_data_rows,]
# Delete rows with unlabelled data from original dataframe
soildataset <- df2[unlabeled_data_rows == FALSE,]

# C. Random forest. PARAMETER TUNING
# Train a random forest with maximumg 20.000 number of trees, and default mtry parameter
set.seed(2)
rf <- randomForest(x=soildataset, 
                   y=droplevels(as.factor(soildataset.samplesinfo[unlabeled_data_rows == FALSE])),
                   importance=TRUE,proximity=TRUE,ntree=15000)
save(rf,file = "rf_ooberror.RData")
load("rf_ooberror.RData")
# Retrieve the OOB error rate per number of trees used in the forest
ooberror <- rf$err.rate[,1]
# Plot the OOB error vs number of trees used
plot(ooberror,pch=19,cex=0.05,xlab="Number of Trees",ylab="OOB error",type="lines")
# Find whih number of trees scored the lower OOB error
minerror_ntrees <- which(ooberror == min(ooberror))
ntrees_new <- 4010 # Set optimum number of trees for future Random Forests trained
# Tune the mtry parameter
set.seed(2)
done = TRUE # avoid running the following line again
if (!done){
  tuned <- tuneRF(x=soildataset,
                y=droplevels(as.factor(soildataset.samplesinfo[unlabeled_data_rows == FALSE])),
                mtryStart = 1,ntreeTry = ntrees_new,stepFactor=20,plot=TRUE)
  save(tuned,file = "tuned.RData")
}
load("tuned.RData")

# D. Train the final model with the best parameters
done <- TRUE
if (!done){
  final_rf <- randomForest(x=soildataset, 
                   y=droplevels(as.factor(soildataset.samplesinfo[unlabeled_data_rows == FALSE])),
                   importance=TRUE,proximity=TRUE,ntree=ntrees_new,mtry=tuned[4,1])
  save(final_rf,file = "finalmodel_randomforest.RData")
}
load("finalmodel_randomforest.RData")

# retrieve labels from training set
labels <- droplevels(as.factor(soildataset.samplesinfo[unlabeled_data_rows == FALSE]))
# Create a hierarchical tree using the proximity matrix from the random forest 
createhierarchicaltree(as.dist(1-final_rf$proximity),labels,title=NULL)

# E. Generate the variable importance plot
varImpPlot(final_rf,cex=0.4,main=NULL)

# F. Retrieve the 20 most important genes according to MeanDecreaseAccuracy
genes_importance_sorted = sort(rf$importance[,8],decreasing=TRUE)
top20_genes = genes_importance_sorted[1:20]

# G. predict NOY and PLO samples
predictions_unlabeled_data <- predict(final_rf, unlabeled_data)

# H. Study the 10 most important gene activity variables
import10 <- rownames(as.data.frame(genes_importance_sorted[1:10]))
labels_toplot <- droplevels(as.factor(soildataset.samplesinfo[unlabeled_data_rows == FALSE]))

# Do a boxplot for every variable from the list of 10 most important variables
i <- 1
# Loop through the 10 most important variables
for (var in import10){
  if (i == 1){
    plot.new()
    par(mfrow=c(2,2))
  }
  if (i == 5){
    plot.new()
    par(mfrow=c(2,2))
  }
  if (i == 9){
    plot.new()
    par(mfrow=c(1,2))
  }
  j <- 1
  x <- list()
  # Study distribution of the variables among every metal
  for (metal in levels(labels_toplot)){
    t1 <- soildataset[labels_toplot==metal,var]
    x[[j]] <- t1
    j <- j + 1
  }
  boxplot(x,names=levels(labels_toplot),cex.axis=0.6, cex.names=0.6,
          col=rainbow(length(levels(labels_toplot))),main=var)
  i <- i + 1 
}
