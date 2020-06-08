cph.tree<-function(X.train,Y.train,X.test,n.iter,sigma.prior.bart=10,ntree,jump.omega=1,burn.in=1,thinning=1){

library(BayesTree)
library(survival)
X.train<-as.matrix(X.train)
X.test<-as.matrix(X.test)
delta<-Y.train$Status
d<-sum(delta)
n.samples<-dim(X.train)[1]
n.samples.test<-dim(X.test)[1]
ndpost<-100
nskip<-50
time<-as.vector(Y.train$Survival)
order<-order(time)
time<-time[order]
delta<-delta[order]
censored.id<-which(delta==0)
d<-sum(delta)
X.train<-X.train[order,,drop=FALSE]  # RDP added drop=FALSE
n.censored<-length(censored.id)
fit<-survreg(Surv(time,delta)~1,dist="weibull")
eta.zero<-fit$coeff
kk.zero<-fit$scale
capital.lambda.star<-eta.zero*(time)^kk.zero
lambda.star<-eta.zero*kk.zero*(time)^(kk.zero-1)
a<-10

#matrix declaration
sigma.error<-matrix(NA,n.iter,1)
B<-matrix(0,n.samples,1)
omega.cph<-matrix(NA,(n.iter+1),n.samples)
omega.jump.mat<-matrix(NA,(n.iter+1),n.samples)

rel.freq.cph<-matrix(NA,n.iter,dim(X.train)[2])
omega.cph[1,]<-rnorm(n.samples,0,1)
omega.cph.test<-matrix(NA,n.iter,n.samples.test)


log.post.omega<-function(omega.cph){
	A.i<-sum(exp(omega.cph[i:length(time)]))
	B[i]<-(-log(1-exp(omega.cph[i])/(a + A.i)))
	summation<-sum(a*B*capital.lambda.star)
	sum(-summation + delta[i]*log(a*lambda.star[i]*B[i])+dnorm(omega.cph[i],g.of.x[i],sqrt(sigma.error[cph]),log=TRUE))
	}

omega.update<-function(sigma.jump.omega){
	omega.star<-omega.cph[cph+1,]
	omega.star[i]<-rnorm(1,omega.cph[cph+1,i],sqrt(sigma.jump.omega))
	log.post.old<-log.post.omega(omega.cph[cph+1,])
	log.post.star<-log.post.omega(omega.star)
	r<-exp(log.post.star-log.post.old)
	omega.cph.hat<-ifelse(runif(1)<r,omega.star[i],omega.cph[cph+1,i])
	p.jump<-min(r,1)
	
	list(omega.cph.hat=omega.cph.hat,p.jump=p.jump)
	}



#MCMC
for (cph in 1:n.iter){

	bart1<-bart(X.train,omega.cph[cph,],x.test=X.test,sigest=sigma.prior.bart,ntree=ntree,ndpost=ndpost,nskip=nskip,verbose=FALSE)
	sigma.error[cph]<-bart1$sigma[ndpost]^2
	g.of.x<-bart1$yhat.train.mean
	omega.cph[cph+1,]<-bart1$yhat.train.mean
	omega.cph.test[cph,]<-bart1$yhat.test.mean
	
	for (i in 1:length(time)){
		temp2<-omega.update(jump.omega)
		omega.cph[cph+1,i]<-temp2$omega.cph.hat
		omega.jump.mat[cph,i]<-temp2$p.jump
	}
	
  varcount<-apply(bart1$varcount,2,sum)
	rel.freq.cph[cph,]<-varcount/(sum(varcount))
	print(sprintf("DRAW_CPH.TREE_%i",cph))
	}

#post-processing
draws<-seq(burn.in,n.iter,thinning)
omega.cph<-as.data.frame(omega.cph)
names(omega.cph)<-row.names(X.train)
omega.cph.test<-as.data.frame(omega.cph.test)
names(omega.cph.test)<-row.names(X.test)
rel.freq.cph<-as.data.frame(rel.freq.cph)
colnames(rel.freq.cph)<-colnames(X.train)

#predictions
omega.cph.quantiles<-apply(omega.cph[draws,],2,quantile,probs=c(0.05,0.5,0.95))
S.train<-sort((a/(a+exp(omega.cph.quantiles)))^(a*capital.lambda.star),decreasing=TRUE)
S.train<-as.data.frame(cbind(time,t(matrix(S.train,3,length(Y.train$Survival)))))
names(S.train)<-c("t","5%","50%","95%")

list(predictions=list(S.train=S.train),parameters=list(omega.train=omega.cph[draws,],omega.test=omega.cph.test[draws,],sigma.error=sigma.error[draws],omega.jump.mat=omega.jump.mat[draws,]),variable.use=data.frame(rel.freq.cph[draws,]))
}