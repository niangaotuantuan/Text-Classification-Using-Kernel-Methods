## Basic string kernels ##
##----------------------##

# We use the kernlab toolbox where several string kernels are implemented #
library(kernlab)

# Look at the descriptions and explanations

# Example: a simple spectrum kernel
sk <- stringdot(type="spectrum", length=2, normalized=FALSE)
sk('radar','abracadabra')

# Note the unexpected behavior
sk('aa','aa') # should be 1 but is 2
sk('aa','ba') # should be 0 but is 1
sk('aa','ab') # should be 0 and is 0
# This is because the strings contain one additional hidden character at the end!

# Note that all kernels created by stringdots only extract substrings without blocks
# To consider gappy substrings, we use another package
# library(stringkernels)

# Example: a simple gappy substring kernels
# sgk <- gapweightkernel(length=2,lambda=0.75,normalized=FALSE,use_characters=TRUE)

# Note that contrary to the spectrum kernel above, this one only computes a kernel on the characters you see
# sgk('aa','aa') # should be 1 ok
# sgk('aa','aaa') # should be 1*(1+1+0.75) ok



## Application to Reuter datasets ##
##--------------------------------##

# Load the data
data(reuters)
y <- rlabels
x <- reuters

# Visualize with kPCA for a particular kernel
sk <- stringdot(type="spectrum", length=2, normalized=FALSE)
kpc <- kpca(x,kernel=sk,scale=c())
plot(rotated(kpc),col=ifelse(y==levels(y)[1],1,2))

# Note that it is usually useful to normalize the kernel to remove the length effect
sk <- stringdot(type="spectrum", length=2, normalized=TRUE)
kpc <- kpca(x,kernel=sk,scale=c())
plot(rotated(kpc),col=ifelse(y==levels(y)[1],1,2))

# The gappy kernel...
#sgk <- gapweightkernel(length=2,lambda=0.75,normalized=TRUE,use_characters=TRUE)
#kpc <- kpca(x,kernel=sgk,scale=c())
#plot(rotated(kpc),col=ifelse(y==levels(y)[1],1,2))


# Train a SVM and estimate CV error with different kernels

# First we focus on spectrum and bounded kernels with different k
kmax <- 20
errspectrum <- numeric(kmax)
errboundrange <- numeric(kmax)
for (k in seq(kmax)) {
	cat('.')
	sk <- stringdot(type="spectrum", length=k, normalized=TRUE)
	svp <- ksvm(x,y,kernel=sk,scale=c(),cross=5)
	errspectrum[k] <- cross(svp)
	sk <- stringdot(type="boundrange", length=k, normalized=TRUE)
	svp <- ksvm(x,y,kernel=sk,scale=c(),cross=5)
	errboundrange[k] <- cross(svp)	
}

# Look at the results
plot(c(1,kmax),c(0,1),type='n',xlab="k",ylab='error')
lines(errspectrum,col=1)
lines(errboundrange,col=2)
grid()
legend("topleft",c('Spectrum','Bounded'),lty=1,col=seq(2))

# Test the exponential kernel with different lambdas
lambdaseq <- 1.1^seq(10)
nlambda <- length(lambdaseq)
errexponential <- numeric(nlambda)
for (i in seq(nlambda)) {
	cat('.')
	sk <- stringdot(type="exponential", lambda=lambdaseq[i], normalized=TRUE)
	svp <- ksvm(x,y,kernel=sk,scale=c(),cross=5)
	errexponential[i] = cross(svp)
}
plot(c(1,max(lambdaseq)),c(0,1),type='n',xlab="lambda",ylab='error')
lines(lambdaseq,errexponential)
grid()




## Protein localization ##
##----------------------##

# We need a library to parse sequence files in FASTA format
library(seqinr)

# Load localization data from http://www.psort.org/dataset/dataset1_0.txt
protdata <- read.fasta("dataset1_0.txt",seqtype="AA",as.string=TRUE)

# Extract protein localization information
annotation <- getName(protdata)

# We get the location information by parsing the annotation as follows
extractlocationfromannotation <- function(s){strsplit(s,'|',fixed=TRUE)[[1]][3]}
y <- factor(unlist(lapply(annotation,extractlocationfromannotation)))

# Extract the protein sequences
x <- unlist(getSequence(protdata,as.string=TRUE),recursive=FALSE)


## Discriminate membrane (inner or outer) vs non-membrane proteins ##

# We get the label as follows
ymembrane <- factor((y=="Inner") | (y=="Outer"))

# Take a random sample of size 100 for training, a random sample of size 100 for test
ntrain <- 100
ntest <- 100
itrain <- sample(length(x),ntrain)
itest <- sample(seq(length(x))[-itrain],ntest)

# Test spectrum kernel with various values for k and C
klist <- seq(5)
clist <- 3^seq(-1,6)
err <- matrix(0,length(klist),length(clist))

# Loop over k
for (ik in seq(length(klist))) {
	k <- klist[ik]
	cat('k=',k)
	
	# kernel definition
	sk <- stringdot(type="spectrum", length=k, normalized=TRUE)
	
	# Pre-compute the kernel matrices (this will be much faster after)
	cat('.')
	Ktrain <- kernelMatrix(sk,x[itrain])
	cat('.')
	Ktest <- kernelMatrix(sk,x[itest],x[itrain])
	ytrain <- ymembrane[itrain]
	
	# Loop over C
	for (ic in seq(length(clist))) {
		cat('.')
		
		# Train a SVM on training set
		svp <- ksvm(Ktrain,ytrain,type="C-svc",C=clist[ic])
		
		# Predict on test set
		ypred <- predict(svp,as.kernelMatrix(Ktest[,SVindex(svp)]))
		
		# Store accuracy
		err[ik,ic] <- sum(ypred != ymembrane[itest])
	}
	cat('\n')
}

err <- err/ntest
plot(c(min(clist),max(clist)),c(0,0.5),log="x",type='n', xlab="C",ylab="Error",main="Spectrum kernel")
nk <- length(klist)
for (ik in seq(nk)) {
	lines(clist,err[ik,],col=ik,lwd=2)
	}
legend("topright",paste("k =",seq(nk)),col=seq(nk),lwd=rep(2,nk))
grid(col='darkgray')


## Example of multiclass classification ##

sk <- stringdot(type="spectrum", length=2, normalized=TRUE)
Ktrain <- kernelMatrix(sk,x[itrain])
ytrain <- y[itrain]
m <- ksvm(Ktrain,ytrain,type="C-svc")
ypred <- predict(m,as.kernelMatrix(Ktest[,SVindex(m)]))

table(ypred,y[itest])
ypred <- predict(m,Ktest)


