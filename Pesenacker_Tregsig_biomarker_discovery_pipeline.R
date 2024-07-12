### Load libraries

library("stringr")                         
library("dplyr")
library("plyr")
library("readr")
library("edgeR")
library("ISLR")
library("elasticnet")
library("tidyverse")
library("caret")
library("glmnet")
library("MASS")
library(limma)

library(caret)
library(glmnet)
library(mlbench)
library(psych)
library(tidyverse)
library(tidyr)

dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

setwd(PrimaryDirectory)
setwd("../data/")
InputDirectory <- getwd()
setwd(PrimaryDirectory)

###################3
setwd(InputDirectory)
list.files(InputDirectory, ".csv")
log2_noflag_data <- read.csv(".csv")   # Enter normalised log2 data file name within ""
log2_noflag_data <- data.frame(log2_noflag_data)
str(log2_noflag_data)
Rdata <- read.csv(".csv")  # Enter file name within ""


Sdata <- read.csv("sample.details.csv")  # read metadata (only needed if not included currently in normalised data or groups need renaming)
Jdata <- join(Rdata,Sdata, by="file_name")  # maps metadata onto dataframe
Jdata <- Jdata[-c(1)]                       #removes extra dataframe column from Rdata- check it is not cutting off any GENES!
write.csv(Jdata, 'Jdata.csv') # check joined data is correct

#### read data #################
data <- read.csv(".csv") # Enter file name within if no joining of metadata required, otherwise read Jdata ""
data <- data.frame(data)
data
data <- data[-c(1)]  
str(data) #note observations and variables
#variable of interest= activity, prediction model to predict activity as function of other variables



############################################################# test/train split and model generation e.g. HC vs active JIA ###

set.seed(222) #random number generation for reproducible results
ind <- sample(2, nrow(data), replace = T, prob = c(0.5, 0.5)) #proportion of data used for training and test
train <- data[ind==1,] 
test <- data[ind==2,]

custom <-  trainControl(method = "LOOCV",
                        verboseIter = T,
                        savePredictions = TRUE) # verboseIter allows us to see whats going on

pre_proc_val <- preProcess(train, method = c("center", "scale"))

train = predict(pre_proc_val, train)
test = predict(pre_proc_val, test)

write.csv(train, "train.csv") # saves what is in your train and test sets
write.csv(test, "test.csv")

set.seed(1234)
elastic <- train(group ~.,
                 train,
                 method = 'glmnet', family= 'binomial',
                 metric = "Accuracy",
                 tuneGrid = expand.grid(alpha= seq(0,1, length=11),
                                        lambda = seq(0.1,1, length=11)),
                 trControl = custom)

plot(elastic)  #x= mixing percentage for alpha between 0-1, lines =regularization parameter for lambda. Higher y=high error 
#can then adjust lambda values to within range with less error
elastic
plot(elastic$finalModel, xvar= "lambda", label = T) #x= logLambda, y= coefficients/ sum squares. Top: variables. 
plot(elastic$finalModel, xvar = 'dev', label=T) #see any over-fitting with larger variance- spot parameters which are not fitting model well
importance_elastic <- plot(varImp(elastic, scale=T)) #importance of each variable (gene)
importance_elastic
varImp(elastic)
fm <- elastic$finalModel

score<- predict(elastic, test, type = 'prob') #predict scores of test set
newcoefs<-coef(elastic$finalModel, elastic$finalModel$lambdaOpt)
mycoefs <- as.data.frame(matrix(newcoefs))
mycoefs <- data.frame(mycoefs)
coefs<- data.frame(matrix(,nrow=38,ncol=0)) #number of genes
coefs<- cbind(coefs, mycoefs)

JIAscore <- cbind(score, test$group)

write.csv(JIAscore, file=".csv")
write.csv(coefs, file=".csv")

saveRDS(elastic, "final_model_.rds")
fm <- readRDS("final_model.rds")

#################################################### plugging in a new cohort into same model as new test set e.g. Inactive:
data <- read.csv(".csv") # Enter file name within ""
data <- data.frame(data)
data
data <- data[-c(1)]  
str(data) #note observations and variables
Rdata_test <- data
Sdata_test <- read.csv(".csv") # read metadata (only needed if not included currently in normalised data or groups need renaming)
Jdata_test <- join(Rdata_test,Sdata_test, by="file_name") # maps metadata onto dataframe
Jdata_test <- Jdata_test[-c(1)]                       #removes extra dataframe column from Rdata- check it is not cutting off Genes!
write.csv(Jdata_test, 'Jdata_test.csv')

test_inactive <- data
pre_proc_val <- preProcess(test_inactive, method = c("center", "scale"))  #scaling
test_inactive = predict(pre_proc_val, test_inactive)

score <- predict(elastic, test_inactive, type= 'prob')
JIAscores <- cbind(score, data$Activity)
write.csv(score, file=".csv")
