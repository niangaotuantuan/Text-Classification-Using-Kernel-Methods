library(kernlab)
## simple example using the spam data set
data(spam)
 
## create test and training set
index <- sample(1:dim(spam)[1])
spamtrain <- spam[index[1:floor(dim(spam)[1]/2)], ]
spamtest <- spam[index[((ceiling(dim(spam)[1]/2)) + 1):dim(spam)[1]], ]
 
## train a support vector machine
filter <- ksvm(type~.,data=spamtrain,kernel="rbfdot",
               kpar=list(sigma=0.05),C=5,cross=3)
filter
 
## predict mail type on the test set
mailtype <- predict(filter,spamtest[,-58]) ## col 58 is type, should be removed when predicting
 
## Check results
table(mailtype,spamtest[,58]) ## the confusion matrix
 
 
## Another example with the famous iris data
data(iris)
 
## Create a kernel function using the build in rbfdot function
rbf <- rbfdot(sigma=0.1)
rbf
 
## train a bound constraint support vector machine
irismodel <- ksvm(Species~.,data=iris,type="C-bsvc",
                  kernel=rbf,C=10,prob.model=TRUE)
 
irismodel
 
## get fitted values
fitted(irismodel)
 
## Test on the training set with probabilities as output
predict(irismodel, iris[,-5], type="probabilities")
 
