FDR<-function(tree.output,FDR=0.05,plot=TRUE){

usage.mat<-as.matrix(tree.output$variable.use)
post.prob<-apply(usage.mat,2,mean)
sorted<-sort(post.prob,decreasing=TRUE,index.return=TRUE)
cumsum<-cumsum(sorted$x)
sel.features<-which(cumsum<FDR)
names.sel.features<-names(cumsum)[sel.features]
last<-sorted$ix[length(names.sel.features)]
cutoff<-post.prob[last]

if(plot==TRUE){
	names(post.prob)<-seq(1,length(post.prob),1)
	barplot(post.prob,col="gray",ylab="Posterior Probability",xlab="Covariates",las=2,cex=0.3,cex.axis=0.8)
	abline(h=cutoff, col="red", lwd=2)
	}

list(post.prob=as.matrix(post.prob),cutoff=cutoff,selected=names.sel.features)
}
