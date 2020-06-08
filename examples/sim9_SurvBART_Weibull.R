setwd("~/Documents/mallick/prophaz2/examples/")
dat <- read.csv('../data/sim9.csv')
source('../competingmethods/SurvBART/wei.tree.R')

# Must expand categorical variables to dummy variables
X.train <- cbind(dat$x1,model.matrix(~-1 + dat$x2))
Y.train <- dat[,c('y','cens')]
X.test <- X.train[1:2,,drop=FALSE]

colnames(Y.train) <- c('Survival','Status')

res <- wei.tree.out<-wei.tree(X.train=X.train,Y.train=Y.train,
                              X.test=X.test,n.iter=500,
                              ntree=5,burn.in=250,thinning=5)

# I believe this is the baseline survival function
plot(res$predictions$S.train$t,res$predictions$S.train$`50%`,type='l')
                      

omegameans <- apply(res$parameters$omega.train,2,mean)

# How I tried to install the Design package (which was archived)
#  It fails due to the fact that it does not have a proper namespace file
setwd('../competingmethods/')
url <- "https://cran.r-project.org/src/contrib/Archive/Design/Design_2.3-0.tar.gz"
pkgFile <- "Design_2.3-0.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)




x <- read.table('../competingsoftware/SurvBART/X.train.txt')
y <- read.table('../competingsoftware/SurvBART/Y.train.txt')


setwd('../competingsoftware/SurvBART/')
X.train<-read.table("X.train.txt")
X.test<-read.table("X.test.txt")
Y.train<-read.table("Y.train.txt")
Y.test<-read.table("Y.test.txt")

#modelling
source("cph.tree.R",echo=FALSE)
cph.tree.out<-cph.tree(X.train=X.train,Y.train=Y.train,X.test=X.test,n.iter=500,ntree=5,burn.in=250,thinning=5)
