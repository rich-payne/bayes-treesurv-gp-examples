wei.tree<-function(X.train,Y.train,X.test,n.iter,sigma.prior.bart=10,ntree,tau.zero=1,k.zero=1,sigma.jump.omega=1,sigma.jump.theta=1,init.theta=NULL,burn.in=1,thinning=1){

library(BayesTree)
X.train<-as.matrix(X.train)
X.test<-as.matrix(X.test)
delta<-Y.train$Status
d<-sum(delta)
n.samples<-dim(X.train)[1]
n.samples.test<-dim(X.test)[1]
ndpost<-100
nskip<-50
time<-as.vector(Y.train$Survival)
sigma.error<-matrix(NA,n.iter,1)
g.of.x.train<-matrix(NA,n.iter,n.samples)
g.of.x.test<-matrix(NA,n.iter,n.samples.test)
rel.freq.wei<-matrix(NA,n.iter,dim(X.train)[2])
omega.mat<-matrix(NA,(n.iter+1),n.samples)
omega.mat[1,]<-rnorm(n.samples,0,1)
omega.mat.test<-matrix(NA,(n.iter+1),n.samples.test)
omega.jump.mat<-matrix(NA,(n.iter+1),n.samples)
omega.jump.mat.test<-matrix(NA,(n.iter+1),n.samples.test)
theta<-matrix(NA,(n.iter+1),1)
theta.jump.mat<-matrix(NA,n.iter)
if(is.null(init.theta)==TRUE){
	library(Design)
	init.theta<-exp(psm(Surv(Y.train$Survival,Y.train$Status)~as.matrix(X.train),dist="weibull",iter=100)$scale)}
theta[1]<-init.theta

log.post.omega<-function(theta,d,delta,omega,time,sigma.error,g.of.x){
		sum(theta*d + delta*omega + delta*(exp(theta)-1)*log(time)-exp(omega)*time^exp(theta))+
		sum(dnorm(omega,g.of.x,sigma.error,log=TRUE))
	}

omega.update<-function(sigma.jump.omega){
		omega.star<-rnorm(1,omega.mat[wei,ii],sigma.jump.omega)
		log.post.old<-log.post.omega(theta[wei],d,delta[ii],omega.mat[wei,ii],time[ii],sigma.error[wei],g.of.x.train[wei,ii])
		log.post.star<-log.post.omega(theta[wei],d,delta[ii],omega.star,time[ii],sigma.error[wei],g.of.x.train[wei,ii])
		r<-exp(log.post.star-log.post.old)
		omega<-ifelse(runif(1)<r,omega.star,omega.mat[wei,ii])
		if(is.na(omega)==TRUE) omega<-1
		p.jump<-min(r,1)

		list(omega=omega,p.jump=p.jump)
	}


log.post.theta<-function(theta,d,delta,omega,time,tau.zero,k.zero){
		summation1<-0
		summation2<-0
		for (ii in 1:length(time)){
			temp1<-((delta[ii]*omega[ii])+delta[ii]*(exp(theta)-1)*log(time[ii]))
			summation1<-summation1+temp1
			temp2<-(exp(omega[ii])*time[ii]^exp(theta))
			summation2<-summation2+temp2
		}
		sum(theta*d + summation1 - summation2 + theta*(tau.zero - 1) - (k.zero*exp(theta)) + theta)
	}


theta.update<-function(sigma.jump.theta){
		theta.star<-rnorm(1,theta[wei],sigma.jump.theta)
		log.post.old<-log.post.theta(theta[wei],d,delta,omega.mat[wei+1,],time,sigma.error[wei],g.of.x.train[wei,])
		log.post.star<-log.post.theta(theta.star,d,delta,omega.mat[wei+1,],time,sigma.error[wei],g.of.x.train[wei,])
		r<-exp(log.post.star-log.post.old)
		theta<-ifelse(runif(1)<r,theta.star,theta[wei])
		if(is.na(theta)==TRUE) theta<-theta[wei]
		p.jump.theta<-min(r,1)
	
		list(theta=theta,p.jump.theta=p.jump.theta)
	}



#MCMC
for (wei in 1:n.iter){
		bart1<-bart(X.train,omega.mat[wei,],x.test=X.test,sigest=sigma.prior.bart,ntree=ntree,ndpost=ndpost,nskip=nskip,verbose=FALSE)
		sigma.error[wei]<-sqrt(bart1$sigma[ndpost])
		g.of.x.train[wei,]<-bart1$yhat.train.mean
		omega.mat.test[wei,]<-bart1$yhat.test.mean

		for (ii in 1:n.samples){
			temp<-omega.update(sigma.jump.omega)
			omega.mat[wei+1,ii]<-temp$omega
			omega.jump.mat[wei+1,ii]<-temp$p.jump
		}


		temp2<-theta.update(sigma.jump.theta)
		theta[wei+1]<-temp2$theta
		theta.jump.mat[wei]<-temp2$p.jump.theta

		varcount<-apply(bart1$varcount,2,sum)
		rel.freq.wei[wei,]<-varcount/(sum(varcount))
		print(sprintf("DRAW_WEIBULL.TREE_%i",wei))
	}


#post-processing
draws<-seq(burn.in,n.iter,thinning)
omega.mat<-as.data.frame(omega.mat)
names(omega.mat)<-row.names(X.train)
omega.mat.test<-as.data.frame(omega.mat.test)
names(omega.mat.test)<-row.names(X.test)
rel.freq.wei<-as.data.frame(rel.freq.wei)
colnames(rel.freq.wei)<-colnames(X.train)

#predictions
omega.wei.quantiles<-apply(omega.mat[draws,],2,quantile,probs=c(0.05,0.5,0.95))
exp.theta<-exp(quantile(theta[draws],probs=c(0.05,0.5,0.95)))
S.train<-sort(exp(-exp(omega.wei.quantiles)*(time^exp.theta)),decreasing=TRUE)
S.train<-as.data.frame(cbind(sort(time),t(matrix(S.train,3,length(Y.train$Survival)))))
names(S.train)<-c("t","5%","50%","95%")

list(predictions=list(S.train=S.train),parameters=list(omega.train=omega.mat[draws,],omega.jump.mat=omega.jump.mat[draws,],omega.test=omega.mat.test[draws,],theta=theta[draws],theta.jump.mat=theta.jump.mat[draws],sigma.error=sigma.error[draws]),variable.use=data.frame(rel.freq.wei[draws,]))
}