##############
## CPH-TREE ##
##############
*Usage

cph.tree(X.train,Y.train,X.test,n.iter,sigma.prior.bart=10,ntree,jump.omega=1,burn.in=1,thinning=1)


*Arguments:

X.train 
Explanatory variables for training (in sample) data. Must be a matrix (typeof double) with (as usual) rows corresponding to 
observations and columns to variables. Note that for a categorical variable you need to use dummies and if there are more than two 
categories, you need to put all the dummies in (unlike linear regression). 
Function will generate draws of f(x) for each x which is a row of x.train.  

Y.train
Dependent variable for training (in sample) data.
Must be a vector (typeof double) with length equal to the number of observations (equal to the number of rows of X.train).

X.test 
Explanatory variables for test (out of sample) data.
Must be a (typeof double) matrix with the same number of columns as X.train.
Function will generate draws of f(x) for each x which is a row of X.test. 

n.iter
Number of MCMC iterations

sigma.prior.bart
The prior for the error variance (sigma^2) is inverted chi-squared (the standard conditionally conjugate prior). Most applications worked 
well with this argument set to 10 (default), however we recommend users to check the sensitivity for large data sets.

ntree
The number of trees in the sum.

jump.omega
"Jump" for the proposal distribution of omega used in the Metropolis-Hastings algorithm.

burn.in
Burning.in parameter for the posterior distribution. Default is set to one, i.e., no burning-in.

thinning
Thinning parameter for the posterior distribution. Default is set to one, i.e., no thinning.


*Value:

Output is divided in 3 slots: predictions, parameters, and variable-use.

$predictions:

	S.train
	Posterior survival distribution estimated (after burning-in and thinning) for training samples along with its 90% credibility interval.

$parameters:

	omega.train
	A matrix with (n.iter after burning-in and thinning) rows and nrow(X.train) columns corresponding to draws from the posteriors of
	the vector of omegas for training (in sample) data. 

	omega.jump.mat
	A matrix with (n.iter after burning-in and thinning) rows and nrow(X.train) columns corresponding to the probability of accepting the proposal 
	distribution in each draw. 

	omega.test
	A matrix with (n.iter after burning-in and thinning) rows and nrow(X.test) columns corresponding to draws from the posteriors of
	the vector of omegas for test (out of sample) data.

	sigma.error
	After burning-in and thinning draws of sigma.

$variable.use:
	a matrix with (n.iter after burning-in and thinning) rows and ncol(X.train) columns. Each row is the relative frequency of variable use 
	for a draw. For each variable (corresponding to the columns), the relative frequency is obtained by dividing the total count of the number 
	of times that variable is used in a tree decision rule (over all trees) by the total number of variables used in each draw.


*Examples:
#loading data
X.train<-read.table("X.train.txt")
X.test<-read.table("X.test.txt")
Y.train<-read.table("Y.train.txt")
Y.test<-read.table("Y.test.txt")

#modelling
source("cph.tree.R",echo=FALSE)
cph.tree.out<-cph.tree(X.train=X.train,Y.train=Y.train,X.test=X.test,n.iter=500,ntree=5,burn.in=250,thinning=5)

#FDR
source("FDR.R",echo=FALSE)
FDR.out<-FDR(tree.output=cph.tree.out,FDR=0.05,plot=TRUE)