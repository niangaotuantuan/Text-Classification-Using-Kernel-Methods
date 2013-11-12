library(kernlab) 
### Use kernelMatrix
K <- as.kernelMatrix(crossprod(t(x)))
svp2 <- ksvm(K, y, type="C-svc") 
svp2
 
# test data
xtest <- rbind(matrix(rnorm(20),,2),matrix(rnorm(20,mean=3),,2))
# test kernel matrix i.e. inner/kernel product of test data with
# Support Vectors
Ktest <- as.kernelMatrix(crossprod(t(xtest),t(x[SVindex(svp2), ]))) 
predict(svp2, Ktest)
 
 
