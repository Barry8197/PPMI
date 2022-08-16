library("rpart")
library("rpart.plot")
library("caTools")
library(randomForest)
library(adabag)
library(glmnet)

setwd('~/PPMI_IP/R_Code/')

# Load Normalised Dataset -------------------------------------------------
load('./../Data/preprocessedData/preprocessed_data.RData')

de_genes <- rownames(genes_info %>% data.frame %>% filter((abs(log2FoldChange) > 0.1) & (padj < 0.05)))

data_input <- merge(t(datExpr) , (datMeta %>% select(PD_status)) , by = 0)
rownames(data_input) <- data_input$Row.names
data_input <- data_input[,2:length(colnames(data_input))]

set.seed(123)
split <- sample.split(data_input$PD_status, SplitRatio=0.8)
train <- subset(data_input, split==TRUE)
test <- subset(data_input, split==FALSE)

trainw <- train
testw <- test


# Classification Tree -----------------------------------------------------
CT <- rpart(PD_status~., data=trainw, method='class', 
            control=rpart.control(maxdepth=20))
CT
rpart.plot(CT, box.palette = "RdBu", digits = 3)

testw$pred <- predict(CT, testw, type="class") 
CM <- (table(testw$PD_status, testw$pred))
CM
41/58


# Random Forrest ----------------------------------------------------------
BaggingCT <- randomForest(PD_status~., data=train)
RF <- randomForest(x = train[, colnames(train) != 'PD_status'] , y = train$PD_status)
testw$baggingct <- predict(RF, test)
(CMBagCT <- table(testw$PD_status, testw$baggingct))
53/58

RFCT <- randomForest(PD_status~., data=train, ntree=1000 , mtry = 100)
testw$rfct <- predict(RFCT, test)
(CMRFCT <- table(testw$PD_status, testw$rfct))
54/58


# ADA Boost ---------------------------------------------------------------
adaboostCT <- boosting(PD_status~., data=train,method='class',
                       control=rpart.control(maxdepth=6), boos=TRUE) 
pred_adaCT <- predict(adaboostCT, testw)
table(pred_adaCT$class, testw$PD_status)
tree1 <- adaboostCT$trees[[1]]
rpart.plot(tree1, box.palette = "RdBu", digits = -3, roundint=FALSE)


# Glmnet ------------------------------------------------------------------
expr_fit <- glmnet(x = train[, colnames(train) != 'PD_status'] , y = train$PD_status , family = 'binomial')
plot(expr_fit , label = TRUE)

pred1 <- predict(expr_fit , newx = as.matrix(test[, colnames(test) != 'PD_status'])  , type = "class")
table(pred1[,100] , test$PD_status)

expr_fit <- cv.glmnet(x = as.matrix(train[, colnames(train) != 'PD_status']) , y = train$PD_status , family = 'binomial')
plot(expr_fit , label = TRUE)

pred1 <- predict(expr_fit , newx = as.matrix(test[, colnames(test) != 'PD_status']) , s = 'lambda.min'  , type = "response")
table(pred1 > 0.5 , test$PD_status)
boxplot(pred1 ~ test$PD_status)

expr_coef <- coef(expr_fit , s = 'lambda.min')
non_zero <- data.frame(expr_coef[(expr_coef[,1] != 0) , ])

genes_info[rownames(non_zero) , ]
