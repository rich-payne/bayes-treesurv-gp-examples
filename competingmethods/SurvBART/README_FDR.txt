#########
## FDR ##
#########
*Usage

FDR(tree.output,FDR=0.05,plot=TRUE)

*Arguments:

tree.output
The output from any survival tree functions: aft.tree, wei.tree, or cph.tree.

FDR
The FDR desired. Ranges from 0 to 1.

plot
Plots a barplot with the posterior probabilities of each covariate along with the cutoff value. 


*Value:

post.prob
Posterior probabilities of each covariate.

cutoff
Cutoff value.

selected
Names of the selected covariates.


*Example
#loading data
X.train<-read.table("X.train.txt")
X.test<-read.table("X.test.txt")
Y.train<-read.table("Y.train.txt")
Y.test<-read.table("Y.test.txt")

#modelling
source("aft.tree.R",echo=FALSE)
aft.tree.out<-aft.tree(X.train=X.train,Y.train=Y.train,X.test=X.test,n.iter=5000,ntree=5,burn.in=2500,thinning=5)

#FDR
source("FDR.R",echo=FALSE)
FDR.out<-FDR(tree.output=aft.tree.out,FDR=0.05,plot=TRUE)


