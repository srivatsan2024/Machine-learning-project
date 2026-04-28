#Final Project 1 and 2 done
#Objective of the project is to build predictive models that classify chemical 
#compounds based on their structural properties.

#Load the data set:Since the columns has no proper title, we create them first

colnames <- c("SpMax_L", "J_Dz_e", "nHM", "F01_NN", "F04_CN", 
              "NssssC", "nCb", "C_pct", "nCp", "nO", 
              "F03_CN", "SdssC", "HyWi_B", "LOC", "SM6_L", 
              "F03_CO", "Me", "Mi", "nN_N", "nArNO2", 
              "nCRX3", "SpPosA_B", "nCIR", "B01_CBr", "B03_CCl", 
              "N073", "SpMax_A", "Psi_i_1d", "B04_CBr", "SdO", 
              "TI2_L", "nCrt", "C026", "F02_CN", "nHDon", 
              "SpMax_Bm", "Psi_i_A", "nN", "SM6_Bm", "nArCOOR", 
              "nX", "class")
bioData <- read.csv("biodeg.csv",sep = ";", col.names = colnames, stringsAsFactors = FALSE, header = FALSE)
head(bioData)

#Check is any values are missing
sum(is.na(bioData))
str(bioData)
table(bioData$class)

#Descriptice statisitics
summary(bioData)
summary(bioData[,1:10])
#Variablity
sapply(bioData[,-42], sd)
bioData$class <- as.factor(bioData$class)#This our target variable
head(bioData)



#Now treat Constinous variables with Histogroams, boxplots and Sparse/binary with classs comparision and frequency table
#Remove variable who are near constant variables
library(caret)
removeVar <- nearZeroVar(bioData)
bioData_clean <- bioData[,-removeVar] #It removed 12 variables who were close to a zero-variance
bioData_clean


#Visualization
#Plot the newly cleaned data
library(ggplot2)

#This is a distribution plot of Class variable
ggplot(bioData_clean, aes(x=class, fill = class))+geom_bar() 
#Add text afterwards

#Histogram
par(mfrow=c(4,4))
ggplot(bioData_clean, aes(x=C_pct, fill = class))+ geom_histogram(alpha=0.5, bins = 30, position = "identity")

ggplot(bioData_clean, aes(x=nX, fill = class))+ geom_histogram(alpha=0.5, bins = 15, position = "identity")

ggplot(bioData_clean, aes(x=nO, fill = class))+ geom_histogram(alpha=0.5, bins = 25, position = "identity")


#BoxPlots
ggplot(bioData_clean, aes(x=class, y=C_pct, fill = class))+geom_boxplot()

ggplot(bioData_clean, aes(x=class, y=nX, fill = class))+geom_boxplot()

ggplot(bioData_clean, aes(x=class, y=nO, fill = class))+geom_boxplot()

ggplot(bioData_clean, aes(x=class, y=nHDon, fill = class))+geom_boxplot()

ggplot(bioData_clean, aes(x=class, y=nCIR, fill = class))+geom_boxplot()


#Correlation
install.packages('corrplot')
library(corrplot)

correlationMap <- bioData_clean[,c("C_pct","nO","nCIR","nN","nHDon")]
corrplot(cor(correlationMap), method = "number")

pca <- prcomp(bioData_clean[,-which(names(bioData_clean)=="class")],scale = TRUE)
par(mar=c(1,1,1,1))
plot(pca$x[,1:2], col=as.factor(bioData_clean$class),pch=19)

#If the PCA shows heave overlap, models like logistice regression, Linear SVM will struggle
# and models like random forest, gradient boosting, and nonlineat svm will be a
# better choice.

## Srivatsan's part
library(randomForest)
library(MASS)
library(glmnet)
trainindex <- createDataPartition(bioData_clean$class, p=0.8, list=FALSE)
train <- bioData_clean[trainindex, ]
test  <- bioData_clean[-trainindex, ]


modellog <- glm(class ~ ., data = bioData_clean, family = "binomial")
summary(modellog)
problog <- predict(modellog, newdata = test, type = "response")
predlog <- ifelse(problog > 0.5, "RB", "NRB")
predlog <- factor(predlog, levels = levels(test$class))
conflog <- confusionMatrix(predlog, test$class)
print(conflog)

lda1 <- lda(class ~ ., data = train)
print(lda1)
ldaraw <- predict(lda1, newdata = test)
ldapred <- ldaraw$class
conflda <- confusionMatrix(ldapred, test$class)
print(conflda)

modelrf <- randomForest(class ~ ., data = train, ntree = 500, mtry = 6, importance = TRUE)
print(modelrf)
importance(modelrf)
varImpPlot(modelrf, main = "Random Forest feature 20", n.var = 20)
predrf <- predict(modelrf, newdata = test)
confrf <- confusionMatrix(predrf, test$class)
print(confrf)
library(e1071)

modelsvm <- svm(class ~., data = train, kernel = "radial", cost = 1, gamma = 0.01, probability = TRUE)
print(modelsvm)
predsvm <- predict(modelsvm, newdata = test)
confsvm <- confusionMatrix(predsvm, test$class)
print(confsvm)

accuracy_summary <- data.frame(
  Model    = c("Logistic Regression", "LDA", "Random Forest", "SVM"),
  Accuracy = c(conflog$overall["Accuracy"],
               conflda$overall["Accuracy"],
               confrf$overall["Accuracy"],
               confsvm$overall["Accuracy"])
)
print(accuracy_summary)


ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
modellogtunned <- train(class ~. , data = train, method = "glmnet", family = "binomial", trControl = ctrl, metric = "ROC", tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 0.1, length=20)))
print(modellogtunned)
plot(modellogtunned, main = "Logistic Regression: Lambda Tuning")
lasso <- coef(modellogtunned$finalModel, modellogtunned$bestTune$lambda)
print(lasso)


nzv_cols <- nearZeroVar(train)
if (length(nzv_cols) > 0) {
  train_qda <- train[, -nzv_cols]
  test_qda  <- test[, -nzv_cols]
  cat("Removed near-zero variance columns:", length(nzv_cols), "\n")
} else {
  train_qda <- train
  test_qda  <- test
}
cor_matrix  <- cor(train_qda[, -ncol(train_qda)])
high_cor    <- findCorrelation(cor_matrix, cutoff = 0.8)
if (length(high_cor) > 0) {
  train_qda <- train_qda[, -high_cor]
  test_qda  <- test_qda[, -high_cor]
  cat("Removed highly correlated columns:", length(high_cor), "\n")
}
cat("Features remaining for QDA:", ncol(train_qda) - 1, "\n")
modelqda <- train(class ~. , data = train_qda, method = "qda", trControl = ctrl, metric = "ROC")
print(modelqda)
predqda <- predict(modelqda, newdata = test)
confqda <- confusionMatrix(predqda, test$class)
print(confqda)

modelrftuned <- train(class ~. , data = train, method = "rf", trControl = ctrl, metric = "ROC", tuneGrid = expand.grid(mtry = c(3, 5, 6, 8, 10, 15, 20))) 
print(modelrftuned)
plot(modelrftuned, main = "Random Forest: mtry Tuning")
bestmtry <- modelrftuned$bestTune$mtry
cat("Best mtry:", bestmtry, "\n")
ntreevals <- c(100, 200, 300, 500, 700, 1000)
ntreeacc <- numeric(length(ntreevals))

set.seed(123)
for (i in seq_along(ntreevals)) {
  rftemp      <- randomForest(class ~ ., data  = train, ntree = ntreevals[i], mtry  = bestmtry)
  predtemp    <- predict(rftemp, newdata = test)
  ntreeacc[i] <- mean(predtemp == test$class)
}
ntreeresults <- data.frame(ntree = ntreevals, accuracy = ntreeacc)
print(ntreeresults)
bestntree <- ntreevals[which.max(ntreeacc)]
cat("Best ntree:", bestntree, "\n")

set.seed(123)
modelrffinal <- randomForest(class ~ ., data = train,  ntree  = bestntree, mtry = bestmtry, importance = TRUE)
predrftuned <- predict(modelrffinal, newdata = test)
confrftuned <- confusionMatrix(predrftuned, test$class)

modelsvmtuned <- train(class ~ ., data  = train, method = "svmRadial", trControl = ctrl, metric    = "ROC", tuneGrid  = expand.grid(C = c(0.1, 0.5, 1, 2, 5, 10), sigma = c(0.001, 0.005, 0.01, 0.05, 0.1)))
## print(modelsvmtuned)
## plot(modelsvmtuned)
predsvmtuned <- predict(modelsvmtuned, newdata = test)
confsvmtuned <- confusionMatrix(predsvmtuned, test$class)
##print(confsvmtuned)

