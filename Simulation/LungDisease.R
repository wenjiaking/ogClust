rm(list=ls()) 
library(sparcl) #to implement sparse K means
library(adaHuber)
library(parallel)
library(Brobdingnag)
library(glmnet)
library(hms)

source("FUNSpool.R")

load("Data/GSE47460_1_GSE47460_2_GSE37147_GSE42057.Rdata")
g.mean<-apply(X.GSE47460_1,1,mean)
cut.mean=quantile(g.mean,probs=0.5)
X.GSE47460_1=X.GSE47460_1[g.mean>cut.mean,] # remove 50% lowest mean expression genes

g.sd=apply(X.GSE47460_1,1, sd) #cv in original scale
cut.sd=quantile(g.sd,probs=0.5)
X.GSE47460_1=X.GSE47460_1[g.sd>=cut.sd,] # further remove 50% lowest variance genes
dim(X.GSE47460_1) #3240  331

X.GSE47460_1<-t(X.GSE47460_1)
X.GSE47460_1<-scale(X.GSE47460_1)

Z.GSE47460_1<-as.matrix(Z.GSE47460_1)
Z.GSE47460_1[,1]<-(Z.GSE47460_1[,1]-mean(Z.GSE47460_1[,1]))/sd(Z.GSE47460_1[,1])
Y1<-scale(Y.GSE47460_1)

X.GSE47460_2=t(X.GSE47460_2)
X.GSE47460_2=X.GSE47460_2[,match(colnames(X.GSE47460_1),colnames(X.GSE47460_2))] #dim 87 1620
X.GSE47460_2<-scale(X.GSE47460_2)
Z.GSE47460_2<-as.matrix(Z.GSE47460_2)
Z.GSE47460_2[,1]<-(Z.GSE47460_2[,1]-mean(Z.GSE47460_2[,1]))/sd(Z.GSE47460_2[,1])
Y2<-scale(Y.GSE47460_2)

######Tune parameters for WJL############

#-----------tune w ------
set.seed(12315)
mod.kmeans<-KMeansSparseCluster(X.GSE47460_1,K = 3,nstart = 50)
center<-matrix(NA,nrow=ncol(X.GSE47460_1),ncol=3)
center[,1]<-apply(X.GSE47460_1[which(mod.kmeans[[20]]$Cs==1),],2,mean)
center[,2]<-apply(X.GSE47460_1[which(mod.kmeans[[20]]$Cs==2),],2,mean)
center[,3]<-apply(X.GSE47460_1[which(mod.kmeans[[20]]$Cs==3),],2,mean)

s_G<-300
K=3

w_vector<-seq(0,1,0.1) 

wcs.lambda<-mclapply(1:length(w_vector),function(i){
  w<-w_vector[i]
  w<-(s_G*w)/(s_G*w+1-w)
  res.lambda<-region_lambda(lambda1 =40,lambda2=0,iteration = 20,Y=Y1,G=X.GSE47460_1,X=as.matrix(Z.GSE47460_1),center=center,w=w,K=K)
  return(res.lambda)
},mc.cores=51)

set.seed(12315)
index.sample<-sample(rep(1:10,length.out=nrow(X.GSE47460_1)))

res1<-mclapply(1:length(w_vector),function(ind_w){
  w<-w_vector[ind_w]
  w<-(s_G*w)/(s_G*w+1-w)
  lambda.vector<-wcs.lambda[[ind_w]]$lambda
  temp.R2<-c()
  for(ind_lambda in 1:length(lambda.vector)){
    lambda<-lambda.vector[ind_lambda]
    mod.all<- ogclust_con_Yujia(x=Z.GSE47460_1,G=t(X.GSE47460_1),y=Y1,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=3,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
    cluster1<-predict_test(mod.all,K=3,D.test=as.matrix(X.GSE47460_1),X1=as.matrix(Z.GSE47460_1),p=ncol(X.GSE47460_1))
    cluster1<-cluster1$clus
    centers<-cbind(apply(X.GSE47460_1[which(cluster1==1),],2,mean),
                   apply(X.GSE47460_1[which(cluster1==2),],2,mean),
                   apply(X.GSE47460_1[which(cluster1==3),],2,mean))
    cluster<-c()
    num1<-c()
    for(ind_cv in 1:10){
      ind_sample<-which(index.sample==ind_cv)
      G.test1<-X.GSE47460_1[ind_sample,]
      X.test1<-Z.GSE47460_1[ind_sample,]
      Y.test1<-Y1[ind_sample]
      G.train1<-X.GSE47460_1[(1:nrow(X.GSE47460_1))[-ind_sample],]
      X.train1<-Z.GSE47460_1[(1:nrow(X.GSE47460_1))[-ind_sample],]
      Y.train1<-Y1[(1:nrow(X.GSE47460_1))[-ind_sample]]
      mod<-ogclust_con_Yujia(x=X.train1,G=t(G.train1),y=Y.train1,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=3,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
      num1[ind_cv]<-sum(apply(mod$result_list$mu,1,function(x){length(unique(x))})!=1)
      res1<-predict_test(mod,K=3,D.test=as.matrix(G.test1),X1=as.matrix(X.test1),p=ncol(G.test1))
      distance.ave<-apply(centers,2,function(x){
        distance<-apply(G.test1,1,function(y){
          sqrt(sum((x-y)^2))
        })
        dis.ave<-rep(NA,K)
        for(i in 1:K){
          dis.ave[i]<-mean(distance[which(res1$clus==i)])
        }
        return(dis.ave)
      })
      rownames(distance.ave)<-paste("cluster",1:K)
      colnames(distance.ave)<-paste("center",1:K)
      
      cluster.index<-list()
      index1<-c()
      for(i in 1:K){
        cluster.index[[i]]<-which(res1$clus==i)
        index1[i]<-!is.na(distance.ave[i,])[1]
      }
      for(i in 1:K){
        if(index1[i]){
          res1$clus[cluster.index[[i]]]<-which.min(distance.ave[i,])
        }
      }
      cluster[ind_sample]<-res1$clus
      #cluster[ind_sample]<-apply(distance.ave,1,which.min)[res1$clus]
    }
    res<-Rsquare(cluster = cluster,Y = Y1,X = Z.GSE47460_1,G = t(X.GSE47460_1))
    temp.R2[ind_lambda]<-res$outcome
  }
  
  #return(temp.R2)
  return(res)
  
  #res$gene
},mc.cores = 51)

R2.outcome=sapply(res1,max)
w_vector[which.max(R2.outcome)]
#-----------------------------select lambda----------------------
w0=w_vector[which.max(R2.outcome)] 
w<-(s_G*w0)/(s_G*w0+1-w0)

lambda_vector1<-seq(15,50,length.out = 22)
set.seed(12315)
index.sample<-sample(rep(1:10,length.out=nrow(X.GSE47460_1)))

wcs<-mclapply(1:length(lambda_vector1),function(ind_lambda){
  lambda<-lambda_vector1[ind_lambda]
  mod.all<-ogclust_con_Yujia(x=Z.GSE47460_1,G=t(X.GSE47460_1),y=Y1,c_center=center,
                             lambda=lambda,v_int=NULL,pi_int=NULL,K=3,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
  select.feature<-as.numeric(apply(mod.all$result_list$mu,1,function(x){length(unique(x))})!=1)
  cluster1<-predict_test(mod.all,K=3,D.test=as.matrix(X.GSE47460_1),X1=as.matrix(Z.GSE47460_1),p=ncol(X.GSE47460_1))
  cluster1<-cluster1$clus
  numk<-length(table(cluster1))
  print(paste0("The number of predicted clusters for lambda=",lambda," is ",numk))
  if(length(table(cluster1))!=3){
    num<-sum(select.feature==1)
    res<-list(R2.outcome=NA,R2.gene=NA,num=num)
    return(res)
  }
  
  centers<-cbind(apply(X.GSE47460_1[which(cluster1==1),],2,mean),
                 apply(X.GSE47460_1[which(cluster1==2),],2,mean),
                 apply(X.GSE47460_1[which(cluster1==3),],2,mean))
  
  
  cluster<-rep(NA,nrow(X.GSE47460_1))
  for(ind_cv in 1:length(unique(index.sample))){
    ind_sample<-which(index.sample==ind_cv)
    G.test11<-X.GSE47460_1[ind_sample,]
    X.test1<-Z.GSE47460_1[ind_sample,]
    Y.test1<-Y1[ind_sample]
    G.train11<-X.GSE47460_1[(1:nrow(X.GSE47460_1))[-ind_sample],]
    X.train1<-Z.GSE47460_1[(1:nrow(X.GSE47460_1))[-ind_sample],]
    Y.train1<-Y1[(1:nrow(X.GSE47460_1))[-ind_sample]]
    
    # label.test1<-label.train[ind_sample]
    # label.train1<-label.train[(1:nrow(X.train))[-ind_sample]]
    #temp.num[ind_w,ind_cv]<-num.vector[index]
    mod<-ogclust_con_Yujia(x=X.train1,G=t(G.train11),y=Y.train1,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=3,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
    
    res1<-predict_test(mod,K=3,D.test=as.matrix(G.test11),X1=as.matrix(X.test1),p=ncol(G.test11))
    
    distance.ave<-apply(centers,2,function(x){
      distance<-apply(G.test11,1,function(y){
        sqrt(sum((x-y)^2))
      })
      dis.ave<-rep(NA,K)
      for(i in 1:K){
        dis.ave[i]<-mean(distance[which(res1$clus==i)])
      }
      return(dis.ave)
    })
    rownames(distance.ave)<-paste("cluster",1:K)
    colnames(distance.ave)<-paste("center",1:K)
    
    cluster.index<-list()
    index1<-c()
    for(i in 1:K){
      cluster.index[[i]]<-which(res1$clus==i)
      index1[i]<-!is.na(distance.ave[i,])[1]
    }
    for(i in 1:K){
      if(index1[i]){
        res1$clus[cluster.index[[i]]]<-which.min(distance.ave[i,])
      }
    }
    cluster[ind_sample]<-res1$clus
    
  }
  res1<-Rsquare(cluster = cluster,Y = Y1,X = Z.GSE47460_1,G = t(X.GSE47460_1))
  R2.outcome<-res1$outcome
  R2.gene<-mean(res1$gene[which(select.feature==1)])
  num<-sum(select.feature==1)
  res<-list(R2.outcome=R2.outcome,R2.gene=R2.gene,num=num,numk=numk)
  return(res)
},mc.cores = 22)

R2.outcome<-unlist(lapply(wcs,function(x){x$R2.outcome}))
R2.gene<-unlist(lapply(wcs,function(x){x$R2.gene}))
num<-unlist(lapply(wcs,function(x){x$num}))
lambda=lambda_vector1[which.max(sqrt(R2.gene*R2.outcome))]

mod = ogclust_con_Yujia(x=Z.GSE47460_1,G=t(X.GSE47460_1),y=Y1,c_center=center,
                        lambda=lambda,v_int=NULL,pi_int=NULL,K=3,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
pred.train<-predict.ogclust.test.select(mod,K=3,D.test=as.matrix(X.GSE47460_1),D.train=t(X.GSE47460_1),X1=Z.GSE47460_1,p=ncol(X.GSE47460_1),O.test = Y1,s_G = s_G,w = w)

res.test1<-predict_test(mod,K=3,D.test=as.matrix(X.GSE47460_2),X1=Z.GSE47460_2,p=ncol(X.GSE47460_1))
mod.WJL=list(mod=mod,res.test1=res.test1)
save(mod.WJL,file="mod.WJL.tune.RData")

#######Tune parameter for ogCLust ######

set.seed(12315)
K=3
lambda_vector1<-seq(0.01,0.1,length.out = 22)
mod.kmeans<-KMeansSparseCluster(X.GSE47460_1,K = K,nstart = 50,silent = T)

cluster<-mod.kmeans[[20]]$Cs
index1<-which(cluster==1)
data1<-data.frame(y=Y1[index1],x1=Z.GSE47460_1[index1,1],x2=Z.GSE47460_1[index1,2])
mod1<-lm(y~.,data=data1)
mod1<-summary(mod1) #train predict model for cluster 1

index2<-which(cluster==2)
data2<-data.frame(y=Y1[index2],x1=Z.GSE47460_1[index2,1],x2=Z.GSE47460_1[index2,2])
mod2<-lm(y~.,data=data2)
mod2<-summary(mod2) #train predict model for cluster 2

index3<-which(cluster==3)
data3<-data.frame(y=Y1[index3],x1=Z.GSE47460_1[index3,1],x2=Z.GSE47460_1[index3,2])
mod3<-lm(y~.,data=data3)
mod3<-summary(mod3) #train predict model for cluster 3

beta0_int = c(mod1$coefficients[1,1],mod2$coefficients[1,1],mod3$coefficients[1,1]) #clsuter specific interscept

data<-data.frame(y=Y1,Z.GSE47460_1)
mod<-lm(y~.,data=data)
mod<-summary(mod)
#mod$coefficients
beta_int<-mod$coefficients[2:3,1] #unified coefficients estimates

mod<-glmnet(X.GSE47460_1,as.factor(cluster),family = "multinomial",lambda=0)
mod1<-coef(mod)
gamma_int<-c(as.numeric(mod1$`1`),as.numeric(mod1$`2`))
sigma2_int<-1
theta_int = c(beta_int, gamma_int, beta0_int, sigma2_int) #initial 

set.seed(12315)
index.sample<-sample(rep(1:10,length.out=nrow(X.GSE47460_1))) #10-fold cv within the leave-one-out train data?

wcs.GM<-mclapply(1:length(lambda_vector1),function(ind_lambda){
  lambda<-lambda_vector1[ind_lambda]
  X.GSE47460_1<-as.data.frame(X.GSE47460_1)
  n1=nrow(X.GSE47460_1) 
  NG=ncol(X.GSE47460_1)
  np=ncol(Z.GSE47460_1) 
  mod.all<-fit.ogClust(n=n1, K=K, np=np, NG=NG, lambda=lambda,
                       alpha=0.5, G=X.GSE47460_1, Y=Y1, X=Z.GSE47460_1, theta_int=theta_int)
  
  beta=mod.all$par$gamma 
  beta<-as.data.frame(beta)
  rownames(beta)[1]<-"intercept"
  rownames(beta)[2:nrow(beta)]<-colnames(X.GSE47460_1)
  select.feature<-as.numeric(apply(beta,1,function(x){sum(x==0)})!=K)[-1] #genes with non-zero weight of at least one clsuter
  #num.vector[ind_lambda]<-sum(select.feature==1)
  cluster1<-mod.all$grp_assign 
  centers<-cbind(apply(X.GSE47460_1[which(cluster1==1),,drop=F],2,mean),
                 apply(X.GSE47460_1[which(cluster1==2),,drop=F],2,mean),
                 apply(X.GSE47460_1[which(cluster1==3),,drop=F],2,mean))
  #print(paste0("Begin the 10-fold CV with temporary lambda=",lambda))
  cluster<-rep(NA,nrow(X.GSE47460_1))
  for(ind_cv in 1:length(unique(index.sample))){
    ind_sample<-which(index.sample==ind_cv)
    G.test1<-X.GSE47460_1[ind_sample,]
    X.test1<-Z.GSE47460_1[ind_sample,]
    Y.test1<-Y1[ind_sample]
    G.train1<-X.GSE47460_1[(1:nrow(X.GSE47460_1))[-ind_sample],]
    X.train1<-Z.GSE47460_1[(1:nrow(X.GSE47460_1))[-ind_sample],]
    Y.train1<-Y1[(1:nrow(X.GSE47460_1))[-ind_sample]]
    n1sub=nrow(G.train1) 
    NGsub=ncol(G.train1) 
    np=ncol(X.train1)
    mod<-fit.ogClust(n=n1sub, K=K, np=np, NG=NGsub, lambda=lambda,
                     alpha=0.5, G=G.train1, Y=Y.train1, X=X.train1, theta_int=theta_int)
    
    
    cluster1<-mod$grp_assign 
    gamma_est_matrix=mod$par$gamma
    
    G.test1.add = as.matrix(cbind(1, G.test1))
    pai_est.num = sapply(1:K, function(k) exp(G.test1.add %*% gamma_est_matrix[,  k, drop = F])/rowSums(exp(G.test1.add %*% gamma_est_matrix)))
    #pai_est.num = sapply(1:K, function(k) exp(G.test %*% gamma_est_matrix[,  k, drop = F])/rowSums(exp(G.test %*% gamma_est_matrix)))
    cluster.test<-apply(pai_est.num,1,which.max) #soft cluster assignment
    
    distance.ave<-apply(centers,2,function(x){
      distance<-apply(G.test1,1,function(y){
        sqrt(sum((x-y)^2))
      })
      dis.ave<-rep(NA,K)
      for(i in 1:K){
        dis.ave[i]<-mean(distance[which(cluster.test==i)])
      }
      return(dis.ave)
    })
    rownames(distance.ave)<-paste("cluster",1:K)
    colnames(distance.ave)<-paste("center",1:K)
    
    cluster.index<-list()
    index1<-c()
    for(i in 1:K){
      cluster.index[[i]]<-which(cluster.test==i)
      index1[i]<-!is.na(distance.ave[i,])[1]
    }
    for(i in 1:K){
      if(index1[i]){
        cluster.test[cluster.index[[i]]]<-which.min(distance.ave[i,])
      }
    }
    cluster[ind_sample]<-cluster.test
    
  }
  if(length(table(cluster))!=3){
    R2.outcome<-NA
    R2.gene<-NA
  }else{
    res1<-Rsquare(cluster = cluster,Y = Y1,X = Z.GSE47460_1,G = t(X.GSE47460_1))
    R2.outcome<-res1$outcome
    R2.gene<-mean(res1$gene[which(select.feature==1)])
  }
  num<-sum(select.feature==1)
  res<-list(R2.outcome=R2.outcome,R2.gene=R2.gene,num=num)
  return(res)
},mc.cores=22)

R2.outcome=sapply(wcs.GM,"[[",1)
R2.gene=sapply(wcs.GM,"[[",2)
lambda<-lambda_vector1[which.max(sqrt(R2.outcome*R2.gene))]
fit.res.num<-fit.ogClust(n=n1, K=K, np=np, NG=NG, lambda=lambda,
                         alpha=0.5, G=X.GSE47460_1, Y=Y1, X=Z.GSE47460_1, theta_int=theta_int)

########### Tune wbounds for Sparse K means###########
K=3
set.seed(12315)
optimal_wbounds<-KMeansSparseCluster.permute(X.GSE47460_1,K=K,wbounds=2:50,silent = T)$bestw
mod.kmeans<-KMeansSparseCluster(X.GSE47460_1,K=K,wbounds=optimal_wbounds,silent = T)

########### Tune lambda PMBC ########
K=3
set.seed(12315)
data.kmeans<-kmeans(as.matrix(X.GSE47460_1),K)
c_size<-data.kmeans$size
n=nrow(X.GSE47460_1)
pi_int<-c_size/n
miu_int<-sapply(1:K,function(x) apply(as.matrix(X.GSE47460_1[data.kmeans$cluster==x,]),2,mean))
sigma_c<-sapply(1:K,function(x) apply(as.matrix(X.GSE47460_1[data.kmeans$cluster==x,]),2,var))
sum_squares<-sapply(1:K,function(x) c_size[x]*sigma_c[,x])
sigma_int<-apply(sum_squares,1,sum)/n

lambda.all=seq(30,60,length.out = 22) 

wcs.MBC=mclapply(1:length(lambda.all),function(ind_lambda) {
  lambda=lambda.all[ind_lambda]
  BIC<-tryCatch({
    BIC<-em_mbc(t(X.GSE47460_1), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                lambda=lambda, K=K, no_init=1, max_iter = 200)$BIC 
  }, error=function(e){
    BIC=Inf
    return(BIC)
  })
  return(BIC)
},mc.cores=22)


lambda=lambda.all[which.min(unlist(wcs.MBC))]

cltrain<-em_mbc(t(X.GSE47460_1), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                lambda=lambda, K=K, no_init=1, max_iter = 200)
modfit=list(mod.GM=fit.res.num,mod.SKM=mod.kmeans,mod.PMBC=cltrain)
save(modfit,file="mod.all.tune.RData")

######Fixed number of genes for WJL ########

wcs.fixgenes.fit=function(n_Genes=400,X=X.GSE47460_1,Y=Y1,Z=Z.GSE47460_1,w_vector,s_G=300,K=3) {
  set.seed(12315)
  mod.kmeans<-KMeansSparseCluster(X,K = K,nstart = 50,silent = T)
  center<-matrix(NA,nrow=ncol(X),ncol=K) #used for initial center
  center[,1]<-apply(X[which(mod.kmeans[[20]]$Cs==1),],2,mean)
  center[,2]<-apply(X[which(mod.kmeans[[20]]$Cs==2),],2,mean)
  center[,3]<-apply(X[which(mod.kmeans[[20]]$Cs==3),],2,mean)
  
  wcs.lambda<-mclapply(1:length(w_vector),function(i){
    w<-w_vector[i]
    w<-(s_G*w)/(s_G*w+1-w)
    res.lambda<-region_lambda(lambda1 =50,lambda2=0,iteration = 20,Y=Y,G=X,X=as.matrix(Z),center=center,w=w,K=K) #what is s_G?
    
    return(res.lambda)
  },mc.cores = length(w_vector))

  lambda.input<-sapply(wcs.lambda,function(x){
    lambda<-x$lambda[min(which(x$num>n_Genes))]
  })
  print(paste(c("The number of selected genes for different weights is ",sapply(wcs.lambda,function(x){lambda<-x$num[min(which(x$num>n_Genes))]})),collapse = " "))

  start_time <- Sys.time()
  set.seed(12315)
  index.sample<-sample(rep(1:10,length.out=nrow(X)))
  R2.outcome.vector<-c()
  for(ind_w in 1:length(w_vector)){
    w<-w_vector[ind_w]
    w<-(s_G*w)/(s_G*w+1-w)
    lambda<-lambda.input[ind_w]
    mod.all<-ogclust_con_Yujia(x=Z,G=t(X),y=Y,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=3,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
    cluster1<-apply(mod.all$result_list$z,1,which.max)
    if(length(table(cluster1))!=3){
      R2.outcome<-NA
      R2.outcome.vector[ind_w]<-R2.outcome
      next
    }
    centers<-cbind(apply(X[which(cluster1==1),],2,mean),
                   apply(X[which(cluster1==2),],2,mean),
                   apply(X[which(cluster1==3),],2,mean))
    cluster<-c()
    for(ind_cv in 1:10){
      #print(paste0(ind_cv,"th-fold CV within the ",ind_w,"th weight and lambda=",lambda))
      ind_sample<-which(index.sample==ind_cv)
      G.test1<-X[ind_sample,]
      X.test1<-Z[ind_sample,]
      Y.test1<-Y[ind_sample]
      G.train1<-X[(1:nrow(X))[-ind_sample],]
      X.train1<-Z[(1:nrow(X))[-ind_sample],]
      Y.train1<-Y[(1:nrow(X))[-ind_sample]]
      mod<-ogclust_con_Yujia(x=X.train1,G=t(G.train1),y=Y.train1,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
      res1<-predict_test(mod,K=K,D.test=as.matrix(G.test1),X1=X.test1,p=ncol(G.test1))
      distance.ave<-apply(centers,2,function(x){
        distance<-apply(G.test1,1,function(y){
          sqrt(sum((x-y)^2)) #calculate the distance to the center trained by all train data (leave-one-out). Does not make sense? should use the center trained by 9-fold train data
        })
        dis.ave<-rep(NA,K)
        for(i in 1:K){
          dis.ave[i]<-mean(distance[which(res1$clus==i)]) #calculate the average distance within each center trained by 9-fold train data?
        }
        return(dis.ave)
      })
      rownames(distance.ave)<-paste("cluster",1:K)
      colnames(distance.ave)<-paste("center",1:K) #distance.ave is a 3x3 matrix

      cluster.index<-list()
      index1<-c()
      for(i in 1:K){
        cluster.index[[i]]<-which(res1$clus==i)
        index1[i]<-!is.na(distance.ave[i,])[1]
      }
      for(i in 1:K){
        if(index1[i]){
          res1$clus[cluster.index[[i]]]<-which.min(distance.ave[i,])
        }
      }
      cluster[ind_sample]<-res1$clus #the test cluster for each fold
      #cluster[ind_sample]<-apply(distance.ave,1,which.min)[res1$clus]
    }
    if(length(table(cluster))!=3){

      R2.outcome<-NA
      R2.outcome.vector[ind_w]<-R2.outcome
      next
    }
    res<-Rsquare(cluster = cluster,Y = Y,X = Z,G = t(X)) #based on the R2 of all trained samples to tune w? wired!
    R2.outcome<-res$outcome
    #just use the R square of outcome fitness conditional on the covariates to select the value of w

    R2.outcome.vector[ind_w]<-R2.outcome

  }
  index<-which.max(R2.outcome.vector)
  lambda<-lambda.input[index]
  w0<-w_vector[index]
  print(paste0("By 10-fold CV, the final selected value of weight is ",w0," and lambda=",lambda))


  w<-(s_G*w0)/(s_G*w0+1-w0)
  lambda_vector1<-seq(15,50,length.out = 22)
  num.all=mclapply(lambda_vector1,function(lambda) {
    mod.temp<-ogclust_con_Yujia(x=as.matrix(Z),G=t(X),y=Y,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
    select.feature<-as.numeric(apply(mod.temp$result_list$mu,1,function(x){length(unique(x))})>1)
    num<-sum(select.feature==1)
    return(num)
  },mc.cores=length(lambda_vector1))
  
  num.all=unlist(num.all)
  print(paste0(c("With w0=",w0,"the number of genes given different lambda is",num.all),collapse = " "))
  lambda<-lambda_vector1[max(which(num.all>n_Genes))]
  print(paste0("With w0= ",w0," select lambda=",lambda," to generate ",n_Genes," genes."))
  
  start_time <- Sys.time()
  
  mod<-ogclust_con_Yujia(x=as.matrix(Z),G=t(X),y=Y,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
  
  end_time <- Sys.time()
  train_time=as_hms(end_time-start_time)
  print(paste0("Time consumed for tuning w and model fitting is ",train_time))
  res<-list(mod=mod,lambda=lambda,w0=w0,train_time=train_time)
  return(res)
}

wcs.fixgenes.pred=function(modfit=modfit.WJL,X=X.GSE47460_1,Y=Y1,Z=Z.GSE47460_1,X.test1=X.GSE47460_2,Z.test1=Z.GSE47460_2,s_G=300) {
  mod=modfit$mod
  w0=modfit$w0
  w<-(s_G*w0)/(s_G*w0+1-w0)
  pred.train<-predict.ogclust.test.select(mod,K=K,D.test=as.matrix(X),D.train=t(X),X1=Z,p=ncol(X),O.test = Y,s_G = s_G,w = w)
  res.test1<-predict_test(mod,K=K,D.test=as.matrix(X.test1),X1=as.matrix(Z.test1),p=ncol(X))
  modfit$res.test1=res.test1
  modfit$pred.train=pred.train
  return(modfit)
}

set.seed(12315)
start_time <- Sys.time()
modfit.WJL=wcs.fixgenes.fit(n_Genes=300,X=X.GSE47460_1,Y=Y1,Z=Z.GSE47460_1,w_vector=seq(0,1,0.1),s_G=300,K=3)
end_time <- Sys.time()
mod.WJL<-wcs.fixgenes.pred(modfit=modfit.WJL)
save(mod.WJL,file="mod.WJL.fix300.Rdata")

start_time <- Sys.time()
modfit.WJL=wcs.fixgenes.fit(n_Genes=400,X=X.GSE47460_1,Y=Y1,Z=Z.GSE47460_1,w_vector=seq(0,1,0.1),s_G=300,K=3)
end_time <- Sys.time()
mod.WJL<-wcs.fixgenes.pred(modfit=modfit.WJL)
save(mod.WJL,file="mod.WJL.fix400.Rdata")

######Fixed number of genes for ogCLust, SKM, and PMBC ########
fixgenes.fit=function(X=X.GSE47460_1,Y=Y1,Z=Z.GSE47460_1,K=3,lambda.GM,wbound,lambda.PMBC) {
  start_time=Sys.time()
  mod.kmeans<-KMeansSparseCluster(X,K = K,nstart = 50,silent = T)
  cluster<-mod.kmeans[[20]]$Cs
  
  index1<-which(cluster==1)
  data1<-data.frame(y=Y[index1],x1=Z[index1,1],x2=Z[index1,2])
  mod1<-lm(y~.,data=data1)
  mod1<-summary(mod1) #train predict model for cluster 1
  
  index2<-which(cluster==2)
  data2<-data.frame(y=Y[index2],x1=Z[index2,1],x2=Z[index2,2])
  mod2<-lm(y~.,data=data2)
  mod2<-summary(mod2)
  
  index3<-which(cluster==3)
  data3<-data.frame(y=Y[index3],x1=Z[index3,1],x2=Z[index3,2])
  mod3<-lm(y~.,data=data3)
  mod3<-summary(mod3) #train predict model for cluster 3
  
  beta0_int = c(mod1$coefficients[1,1],mod2$coefficients[1,1],mod3$coefficients[1,1]) #clsuter specific interscept
  
  data<-data.frame(y=Y,Z)
  mod<-lm(y~.,data=data)
  mod<-summary(mod)
  #mod$coefficients
  beta_int<-mod$coefficients[2:3,1] #unified coefficients estimates
  
  mod<-glmnet(X,as.factor(cluster),family = "multinomial",lambda=0)
  mod1<-coef(mod)
  gamma_int<-c(as.numeric(mod1$`1`),as.numeric(mod1$`2`))
  
  sigma2_int<-1
  theta_int = c(beta_int, gamma_int, beta0_int, sigma2_int) #initial 
  #X<-as.data.frame(X)
  
  n1=nrow(X) # number of samples
  NG=ncol(X) # number of genes
  np=ncol(Z) # number of covariates
  
  mod.ogClust<-fit.ogClust(n=n1, K=K, np=np, NG=NG, lambda=lambda.GM,
                           alpha=0.5, G=as.data.frame(X), Y=Y, X=Z, theta_int=theta_int)
  end_time <- Sys.time()
  train_time.GM=as_hms(end_time-start_time)
  
  start_time=Sys.time()
  mod.SKM<-KMeansSparseCluster(X,K = K,nstart = 50,wbounds=wbound,silent = T)
  end_time <- Sys.time()
  train_time.SKM=as_hms(end_time-start_time)
  
  start_time=Sys.time()
  data.kmeans<-kmeans(as.matrix(X),K)
  c_size<-data.kmeans$size
  n=nrow(X)
  pi_int<-c_size/n
  miu_int<-sapply(1:K,function(x) apply(as.matrix(X[data.kmeans$cluster==x,]),2,mean))
  
  sigma_c<-sapply(1:K,function(x) apply(as.matrix(X[data.kmeans$cluster==x,]),2,var))
  sum_squares<-sapply(1:K,function(x) c_size[x]*sigma_c[,x])
  sigma_int<-apply(sum_squares,1,sum)/n
  
  cltrain<-em_mbc(t(X), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                  lambda=lambda.PMBC, K=K, no_init=1, max_iter = 200)
  end_time <- Sys.time()
  train_time.PMBC=as_hms(end_time-start_time)
  
  res<-list(mod.GM=mod.ogClust,mod.SKM=mod.SKM,mod.PMBC=cltrain,
            train_time.GM=train_time.GM,train_time.SKM=train_time.SKM,train_time.PMBC=train_time.PMBC)
  return(res) 
}

# For GM model:
# lambda<-0.014  select about 300 genes 
# lambda<- 0.004  select about 400 genes

# For SKM model:
# wbound=14.2 select about 300 genes 
# wbound=16.5 select about 400 genes 

# For PMBC model:
# lambda=101 select about 300 genes 
# lambda=94.28571 select about 400 genes 

set.seed(12315)
start_time <- Sys.time()
modfit=fixgenes.fit(lambda.GM=0.004,wbound=16.5,lambda.PMBC=94.28571)
end_time <- Sys.time()
save(modfit,file="mod.all.fix400.Rdata") 
set.seed(12315)
start_time <- Sys.time()
modfit=fixgenes.fit(lambda.GM=0.014,wbound=14.2,lambda.PMBC=101)
end_time <- Sys.time()
save(modfit,file="mod.all.fix300.Rdata") 

########Generate Figures###############
readIN.EXT=function(ngenes=100) {
  load("Data/GSE47460_1_GSE47460_2_GSE37147_GSE42057.Rdata")
  g.mean<-apply(X.GSE47460_1,1,mean)
  cut.mean=quantile(g.mean,probs=0.5)
  X.GSE47460_1=X.GSE47460_1[g.mean>cut.mean,] # remove 50% lowest mean expression genes
  
  g.sd=apply(X.GSE47460_1,1, sd) #cv in original scale
  cut.sd=quantile(g.sd,probs=0.5)
  X.GSE47460_1=X.GSE47460_1[g.sd>=cut.sd,] # further remove 50% lowest variance genes
  dim(X.GSE47460_1) #3240  331
  
  X.GSE47460_1<-t(X.GSE47460_1)
  X.GSE47460_1<-scale(X.GSE47460_1)
  
  Z.GSE47460_1<-as.matrix(Z.GSE47460_1)
  Z.GSE47460_1[,1]<-(Z.GSE47460_1[,1]-mean(Z.GSE47460_1[,1]))/sd(Z.GSE47460_1[,1])
  Y1<-scale(Y.GSE47460_1)

  load(paste0("mod.WJL.fix",ngenes,".RData"))
  load(paste0("mod.all.fix",ngenes,".RData"))
  
  X.GSE47460_2=t(X.GSE47460_2)
  X.GSE47460_2=X.GSE47460_2[,match(colnames(X.GSE47460_1),colnames(X.GSE47460_2))] #dim 87 1620
  X.GSE47460_2<-scale(X.GSE47460_2)
  Z.GSE47460_2<-as.matrix(Z.GSE47460_2)
  Z.GSE47460_2[,1]<-(Z.GSE47460_2[,1]-mean(Z.GSE47460_2[,1]))/sd(Z.GSE47460_2[,1])
  Y2<-scale(Y.GSE47460_2)
 
  plotlist=EXTplot.gen(X=X.GSE47460_1,Y=Y1,Z=Z.GSE47460_1,mod.WJL=mod.WJL,modfit=modfit,ngenes=ngenes,K=3,
                       X.test1=X.GSE47460_2,Z.test1=Z.GSE47460_2,Y.test1=Y2)
  return(plotlist)
}

library(pheatmap)
library(gplots)
library(tibble)
library(ggplotify)
library(ggpubr)
library(gridExtra)

EXTplot.gen=function(X=X.GSE47460_1,Y=Y1,Z=Z.GSE47460_1,mod.WJL=mod.WJL,modfit=modfit,ngenes=100,K=3,
                     X.test1=X.GSE47460_2,Z.test1=Z.GSE47460_2,Y.test1=Y2) {
  
  mod=mod.WJL$mod
  res.test1=mod.WJL$res.test1
  pred.train=mod.WJL$pred.train
  select.feature<-as.numeric(apply(mod$result_list$mu,1,function(x){length(unique(x))})>1)
  
  num.yujia.num<-sum(select.feature==1) #300 genes for tuned parameters without screening, while 330 with screening
  
  cluster.train = factor(pred.train$clus)
  
  cluster.train1<-cluster.train
  index<-order(c(median(Y[cluster.train1==1]),median(Y[cluster.train1==2]),
                 median(Y[cluster.train1==3])),decreasing = F)
  cluster.train[which(cluster.train1==index[1])]<-1
  cluster.train[which(cluster.train1==index[2])]<-2
  cluster.train[which(cluster.train1==index[3])]<-3
  
  
  cluster<-cluster.train
  cluster<-as.numeric(as.character(cluster))
  cluster[which(cluster.train==1)]<-"low"
  cluster[which(cluster.train==2)]<-"medium"
  cluster[which(cluster.train==3)]<-"high"
  data=data.frame(Y=Y,cluster=cluster)
  data$cluster<-factor(data$cluster,levels = c("low","medium","high"))
  
  
  p1<-ggplot(data)+aes(x=cluster,y=Y,fill=cluster)+geom_boxplot()+theme_bw()+
    scale_fill_manual(values =  c("#61D04F", "#2297E6", "#28E2E5"))+ylab("Fev1")+
    #ggtitle("ogClustWJL without screening")+
    #ggtitle(title)+
    theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"),
          strip.text.x = element_text(size = 40),legend.text = element_text(size=25),
          legend.title = element_text(size=25),
          legend.key.size = unit(1.5, 'cm'),legend.position="none",panel.border = element_rect(colour = "black",size=2))
  
  
  mod1<-aov(Y~cluster,data=data)
  mod1<-summary(mod1)
  pvalue.train<-mod1[[1]]$`Pr(>F)`[1]
  RMSE.wcs=sqrt(sum((pred.train$Y.hard-Y)^2)/length(Y)) #0.5750711 without sc while  0.5751258/0.5752444 with sc
  
  data1<-data.frame(Y=Y,cluster=as.factor(cluster),x1=Z[,1],x2=Z[,2])
  res.lm <- lm(Y~ cluster+x1+x2, data = data1)
  RSS<-sum(res.lm$residuals^2)
  res.lm <- lm(Y ~ 1+x1+x2, data = data1)
  TSS<-sum(res.lm$residuals^2)
  R2.outcome.wcs<-1-RSS/TSS 
  
  R2.gene.wcs<-c()
  for(i in 1:ncol(X)){
    data1<-data.frame(gene=X[,i],cluster=as.factor(cluster))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.wcs[i]<-1-RSS/TSS 
  }
  R2.gene.wcs<-R2.gene.wcs[which(select.feature==1)]
  
  #For test data 1
  cluster.test<-res.test1$clus
  
  cluster.test1<-cluster.test
  index<-order(c(median(Y[cluster.train1==1]),median(Y[cluster.train1==2]),
                 median(Y[cluster.train1==3])),decreasing = F)
  cluster.test[which(cluster.test1==index[1])]<-1
  cluster.test[which(cluster.test1==index[2])]<-2
  cluster.test[which(cluster.test1==index[3])]<-3
  
  
  cluster<-cluster.test
  cluster<-as.numeric(as.character(cluster))
  cluster[which(cluster.test==1)]<-"low"
  cluster[which(cluster.test==2)]<-"medium"
  cluster[which(cluster.test==3)]<-"high"
  data<-data.frame(Y=Y.test1,cluster=cluster)
  data$cluster<-factor(data$cluster,levels = c("low","medium","high"))
  
  p2<-ggplot(data)+aes(x=cluster,y=Y,fill=cluster)+geom_boxplot()+theme_bw()+
    scale_fill_manual(values =  c("#61D04F", "#2297E6", "#28E2E5"))+ylab("Fev1")+
    theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"),
          strip.text.x = element_text(size = 40),legend.text = element_text(size=25),
          legend.title = element_text(size=25),
          legend.key.size = unit(1.5, 'cm'),legend.position="none",panel.border = element_rect(colour = "black",size=2))
  
  mod1<-aov(Y~cluster,data=data)
  mod1<-summary(mod1)
  pvalue.test<-mod1[[1]]$`Pr(>F)`[1]
  RMSE.test1.wcs=sqrt(sum((res.test1$Y-Y.test1)^2)/length(Y.test1)) 
  
  data1<-data.frame(Y=Y.test1,cluster=as.factor(cluster),x1=Z.test1[,1],x2=Z.test1[,2])
  res.lm <- lm(Y~ cluster+x1+x2, data = data1)
  RSS<-sum(res.lm$residuals^2)
  res.lm <- lm(Y ~ 1+x1+x2, data = data1)
  TSS<-sum(res.lm$residuals^2)
  R2.outcome.test1.wcs<-1-RSS/TSS 
  
  R2.gene.test1.wcs<-c()
  for(i in 1:ncol(X.test1)){
    data1<-data.frame(gene=X.test1[,i],cluster=as.factor(cluster))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.test1.wcs[i]<-1-RSS/TSS 
  }
  R2.gene.test1.wcs<-R2.gene.test1.wcs[which(select.feature==1)]
  
  
  
  p1 <- annotate_figure(p1,top = text_grob("(II.B1) outcome separation (discovery)",face = "bold", size = 20 ))
  p2 <- annotate_figure(p2,top = text_grob("(II.B2) outcome separation (validation)",face = "bold", size = 20 ))
  
  fit.res.num=modfit$mod.GM
  beta=fit.res.num$par$gamma 
  beta<-as.data.frame(beta)
  rownames(beta)[1]<-"intercept"
  rownames(beta)[2:nrow(beta)]<-colnames(X)
  select.feature<-as.numeric(apply(beta,1,function(x){sum(x==0)})!=K)[-1] #genes with non-zero weight of at least one clsuter
  num.Peng.num<-sum(select.feature==1) 
  
  pred.train<-fit.res.num$Y_prd
  cluster.train<-fit.res.num$grp_assign 
  cluster.train1<-cluster.train
  index<-order(c(median(Y[cluster.train1==1]),median(Y[cluster.train1==2]),
                 median(Y[cluster.train1==3])),decreasing = F)
  cluster.train[which(cluster.train1==index[1])]<-1
  cluster.train[which(cluster.train1==index[2])]<-2
  cluster.train[which(cluster.train1==index[3])]<-3
  
  cluster<-cluster.train
  cluster<-as.numeric(as.character(cluster))
  cluster[which(cluster.train==1)]<-"low"
  cluster[which(cluster.train==2)]<-"medium"
  cluster[which(cluster.train==3)]<-"high"
  data<-data.frame(Y=Y,cluster=cluster)
  data$cluster<-factor(data$cluster,levels = c("low","medium","high"))
  data$diagnosis<-rep(c("COPD","ILD"),length.out=nrow(data))
  p4<-ggplot(data)+aes(x=cluster,y=Y,fill=cluster)+geom_boxplot()+theme_bw()+
    scale_fill_manual(values =  c("#61D04F", "#2297E6", "#28E2E5"))+ylab("Fev1")+
    theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"),
          strip.text.x = element_text(size = 40),legend.text = element_text(size=25),
          legend.title = element_text(size=25),
          legend.key.size = unit(1.5, 'cm'),legend.position="none",panel.border = element_rect(colour = "black",size=2))
  
  mod1<-aov(Y~cluster,data=data)
  mod1<-summary(mod1)
  pvalue.train.GM<-mod1[[1]]$`Pr(>F)`[1]
  RMSE.GM=sqrt(sum((pred.train-Y)^2)/length(Y)) 
  
  data1<-data.frame(Y=Y,cluster=as.factor(cluster),x1=Z[,1],x2=Z[,2])
  res.lm <- lm(Y~ cluster+x1+x2, data = data1)
  RSS<-sum(res.lm$residuals^2)
  res.lm <- lm(Y ~ 1+x1+x2, data = data1)
  TSS<-sum(res.lm$residuals^2)
  R2.outcome.GM<-1-RSS/TSS
  
  R2.gene.GM<-c()
  for(i in 1:ncol(X)){
    data1<-data.frame(gene=X[,i],cluster=as.factor(cluster))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.GM[i]<-1-RSS/TSS 
  }
  R2.gene.GM<-R2.gene.GM[which(select.feature==1)]
  
  beta_est=fit.res.num$par$beta
  gamma_est_matrix=fit.res.num$par$gamma
  beta0_est=fit.res.num$par$beta0
  sigma2_est=fit.res.num$par$sigma2
  
  X.test1.add = cbind(1, X.test1)
  pai_est.num = sapply(1:K, function(k) exp(X.test1.add %*% gamma_est_matrix[,  k, drop = F])/rowSums(exp(X.test1.add %*% gamma_est_matrix)))
  
  cluster.test<-apply(pai_est.num,1,which.max)
  
  intercept<-rep(NA,nrow(Z.test1))
  for(i in 1:K){
    intercept[which(cluster.test==i)]<-beta0_est[i]
  }
  pred.test<-as.matrix(scale(Z.test1)) %*% beta_est+intercept
  
  cluster.test1<-cluster.test
  index<-order(c(median(Y[cluster.train1==1]),median(Y[cluster.train1==2]),
                 median(Y[cluster.train1==3])),decreasing = F)
  cluster.test[which(cluster.test1==index[1])]<-1
  cluster.test[which(cluster.test1==index[2])]<-2
  cluster.test[which(cluster.test1==index[3])]<-3
  
  cluster<-cluster.test
  cluster<-as.numeric(as.character(cluster))
  cluster[which(cluster.test==1)]<-"low"
  cluster[which(cluster.test==2)]<-"medium"
  cluster[which(cluster.test==3)]<-"high"
  data<-data.frame(Y=Y.test1,cluster=cluster)
  data$cluster<-factor(data$cluster,levels = c("low","medium","high"))
  
  p5<-ggplot(data)+aes(x=cluster,y=Y,fill=cluster)+geom_boxplot()+theme_bw()+
    scale_fill_manual(values =  c("#61D04F", "#2297E6", "#28E2E5"))+ylab("Fev1")+
    theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"),
          strip.text.x = element_text(size = 40),legend.text = element_text(size=25),
          legend.title = element_text(size=25),
          legend.key.size = unit(1.5, 'cm'),legend.position="none",panel.border = element_rect(colour = "black",size=2))
  
  mod1<-aov(Y~cluster,data=data)
  mod1<-summary(mod1)
  pvalue.test.GM<-mod1[[1]]$`Pr(>F)`[1]
  RMSE.test1.GM=sqrt(sum((pred.test-Y.test1)^2)/length(Y.test1)) #0.8014442 for sc, while 0.8031067 without sc
  
  data1<-data.frame(Y=Y.test1,cluster=as.factor(cluster),x1=Z.test1[,1],x2=Z.test1[,2])
  res.lm <- lm(Y~ cluster+x1+x2, data = data1)
  RSS<-sum(res.lm$residuals^2)
  res.lm <- lm(Y ~ 1+x1+x2, data = data1)
  TSS<-sum(res.lm$residuals^2)
  R2.outcome.test1.GM<-1-RSS/TSS 
  
  R2.gene.test1.GM<-c()
  for(i in 1:ncol(X.test1)){
    data1<-data.frame(gene=X.test1[,i],cluster=as.factor(cluster))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.test1.GM[i]<-1-RSS/TSS 
  }
  R2.gene.test1.GM<-R2.gene.test1.GM[which(select.feature==1)]
  
  print("finish GM plot")
  
  p4 <- annotate_figure(p4,top = text_grob("(II.C1) outcome separation (discovery)",face = "bold", size = 20 ))
  p5 <- annotate_figure(p5,top = text_grob("(II.C2) outcome separation (validation)",face = "bold", size = 20 ))
  
  mod.kmeans=modfit$mod.SKM
  select.feature<-as.numeric(mod.kmeans[[1]]$ws!=0)
  num.SKM.num<-sum(select.feature==1)
  
  weight<-mod.kmeans[[1]]$ws
  cluster.train<-mod.kmeans[[1]]$Cs
  centers<-matrix(NA,ncol(X),3)
  for(i1 in 1:3){
    if(sum(mod.kmeans[[1]]$Cs==i1)>1){
      centers[,i1]<-apply((X)[which(mod.kmeans[[1]]$Cs==i1),],2,mean)
    }else{
      centers[,i1]<-X[which(mod.kmeans[[1]]$Cs==i1),]
    }
    
  }
  
  cluster.train1<-cluster.train
  index<-order(c(median(Y[cluster.train1==1]),median(Y[cluster.train1==2]),
                 median(Y[cluster.train1==3])),decreasing = F)
  cluster.train[which(cluster.train1==index[1])]<-1
  cluster.train[which(cluster.train1==index[2])]<-2
  cluster.train[which(cluster.train1==index[3])]<-3
  
  
  cluster<-cluster.train
  cluster<-as.numeric(as.character(cluster))
  cluster[which(cluster.train==1)]<-"low"
  cluster[which(cluster.train==2)]<-"medium"
  cluster[which(cluster.train==3)]<-"high"
  data<-data.frame(Y=Y,cluster=cluster)
  data$cluster<-factor(data$cluster,levels = c("low","medium","high"))
  
  p7<-ggplot(data)+aes(x=cluster,y=Y,fill=cluster)+geom_boxplot()+theme_bw()+
    scale_fill_manual(values =  c("#61D04F", "#2297E6", "#28E2E5"))+ylab("Fev1")+
    theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"),
          strip.text.x = element_text(size = 40),legend.text = element_text(size=25),
          legend.title = element_text(size=25),
          legend.key.size = unit(1.5, 'cm'),legend.position="none",panel.border = element_rect(colour = "black",size=2))
  
  
  mod1<-aov(Y~cluster,data=data)
  mod1<-summary(mod1)
  pvalue.train.SKM<-mod1[[1]]$`Pr(>F)`[1]
  pred.train.kmeans<-rep(NA,length(Y))
  
  for(i in 1:3){
    data.temp<-data.frame(Y=Y[which(mod.kmeans[[1]]$Cs==i)],Z[which(mod.kmeans[[1]]$Cs==i),])
    mod.lm<-lm(Y~.,data=data.temp)
    pred.train.kmeans[which(mod.kmeans[[1]]$Cs==i)]<-mod.lm$fitted.values
  }
  RMSE.SKM<-sqrt(sum((pred.train.kmeans-Y)^2)/length(Y)) 
  data1<-data.frame(Y=Y,cluster=as.factor(cluster),x1=Z[,1],x2=Z[,2])
  res.lm <- lm(Y~ cluster+x1+x2, data = data1)
  RSS<-sum(res.lm$residuals^2)
  res.lm <- lm(Y ~ 1+x1+x2, data = data1)
  TSS<-sum(res.lm$residuals^2)
  R2.outcome.SKM<-1-RSS/TSS
  
  R2.gene.SKM<-c()
  for(i in 1:ncol(X)){
    data1<-data.frame(gene=X[,i],cluster=as.factor(cluster))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.SKM[i]<-1-RSS/TSS 
  }
  R2.gene.SKM<-R2.gene.SKM[which(select.feature==1)]
  
  
  distance<-apply(X.test1,1,function(x){
    distance1<-apply(centers,2,function(y){
      sum(weight*(x-y)^2)
    })
  })
  
  #Z.test<-as.matrix(Z.test)
  cluster.test<-apply(distance,2,which.min)
  
  cluster.test1<-cluster.test
  index<-order(c(median(Y[cluster.train1==1]),median(Y[cluster.train1==2]),
                 median(Y[cluster.train1==3])),decreasing = F)
  cluster.test[which(cluster.test1==index[1])]<-1
  cluster.test[which(cluster.test1==index[2])]<-2
  cluster.test[which(cluster.test1==index[3])]<-3
  
  
  cluster<-cluster.test
  cluster<-as.numeric(as.character(cluster))
  cluster[which(cluster.test==1)]<-"low"
  cluster[which(cluster.test==2)]<-"medium"
  cluster[which(cluster.test==3)]<-"high"
  data<-data.frame(Y=Y.test1,cluster=cluster)
  data$cluster<-factor(data$cluster,levels = c("low","medium","high"))
  
  p8<-ggplot(data)+aes(x=cluster,y=Y,fill=cluster)+geom_boxplot()+theme_bw()+
    scale_fill_manual(values =  c("#61D04F", "#2297E6", "#28E2E5"))+ylab("Fev1")+
    theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"),
          strip.text.x = element_text(size = 40),legend.text = element_text(size=25),
          legend.title = element_text(size=25),
          legend.key.size = unit(1.5, 'cm'),legend.position="none",panel.border = element_rect(colour = "black",size=2))
  
  mod1<-aov(Y~cluster,data=data)
  mod1<-summary(mod1)
  pvalue.test.SKM<-mod1[[1]]$`Pr(>F)`[1]
  
  pred.test.kmeans<-rep(NA,length(Y.test1))
  for(i in 1:3){
    data.temp<-data.frame(Y=Y[which(mod.kmeans[[1]]$Cs==i)],Z[which(mod.kmeans[[1]]$Cs==i),])
    mod.lm<-lm(Y~.,data=data.temp)
    #pred.train.kmeans[which(mod[[index]]$Cs==i)]<-mod.lm$fitted.values
    pred.test.kmeans[which(cluster.test==i)]<-predict(mod.lm,newdata = as.data.frame(as_tibble(Z.test1)[which(cluster.test==i),]))
  }
  
  RMSE.test1.SKM<-sqrt(sum((pred.test.kmeans-Y.test1)^2)/length(Y.test1)) 
  
  data1<-data.frame(Y=Y.test1,cluster=as.factor(cluster),x1=Z.test1[,1],x2=Z.test1[,2])
  res.lm <- lm(Y~ cluster+x1+x2, data = data1)
  RSS<-sum(res.lm$residuals^2)
  res.lm <- lm(Y ~ 1+x1+x2, data = data1)
  TSS<-sum(res.lm$residuals^2)
  R2.outcome.test1.SKM<-1-RSS/TSS
  
  R2.gene.test1.SKM<-c()
  for(i in 1:ncol(X.test1)){
    data1<-data.frame(gene=X.test1[,i],cluster=as.factor(cluster))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.test1.SKM[i]<-1-RSS/TSS 
  }
  R2.gene.test1.SKM<-R2.gene.test1.SKM[which(select.feature==1)]
  
  
  p7 <- annotate_figure(p7,top = text_grob("(II.D1) outcome separation (discovery)",face = "bold", size = 20 ))
  p8 <- annotate_figure(p8,top = text_grob("(II.D2) outcome separation (validation)",face = "bold", size = 20 ))
  print("finish SKM plot")
  
  cltrain=modfit$mod.PMBC
  select.feature<-cltrain$selected.index 
  num.PMBC.num=length(cltrain$selected.index) 
  
  pi_est_train<-cltrain$optimal_result$pi
  mu_est_train<-cltrain$optimal_result$mu
  sigma_est_train<-cltrain$optimal_result$sigma
  grp_assign_train<-apply(cltrain$optimal_result$z,1,which.max)
  cluster.train=grp_assign_train
  cluster.train1<-cluster.train
  index<-order(c(median(Y[cluster.train1==1]),median(Y[cluster.train1==2]),
                 median(Y[cluster.train1==3])),decreasing = F)
  cluster.train[which(cluster.train1==index[1])]<-1
  cluster.train[which(cluster.train1==index[2])]<-2
  cluster.train[which(cluster.train1==index[3])]<-3
  
  
  cluster<-cluster.train
  cluster<-as.numeric(as.character(cluster))
  cluster[which(cluster.train==1)]<-"low"
  cluster[which(cluster.train==2)]<-"medium"
  cluster[which(cluster.train==3)]<-"high"
  data<-data.frame(Y=Y,cluster=cluster)
  data$cluster<-factor(data$cluster,levels = c("low","medium","high"))
  
  p10<-ggplot(data)+aes(x=cluster,y=Y,fill=cluster)+geom_boxplot()+theme_bw()+
    scale_fill_manual(values =  c("#61D04F", "#2297E6", "#28E2E5"))+ylab("Fev1")+
    theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"),
          strip.text.x = element_text(size = 40),legend.text = element_text(size=25),
          legend.title = element_text(size=25),
          legend.key.size = unit(1.5, 'cm'),legend.position="none",panel.border = element_rect(colour = "black",size=2))
  
  mod1<-aov(Y~cluster,data=data)
  mod1<-summary(mod1)
  pvalue.train.mbc<-mod1[[1]]$`Pr(>F)`[1]
  pred.train.mbc<-rep(NA,length(Y))
  
  for(i in 1:3){
    data.temp<-data.frame(Y=Y[which(grp_assign_train==i)],Z[which(grp_assign_train==i),])
    mod.lm<-lm(Y~.,data=data.temp)
    pred.train.mbc[which(grp_assign_train==i)]<-mod.lm$fitted.values
  }
  RMSE.PMBC<-sqrt(sum((pred.train.mbc-Y)^2)/length(Y)) 
  
  data1<-data.frame(Y=Y,cluster=as.factor(cluster),x1=Z[,1],x2=Z[,2])
  res.lm <- lm(Y~ cluster+x1+x2, data = data1)
  RSS<-sum(res.lm$residuals^2)
  res.lm <- lm(Y ~ 1+x1+x2, data = data1)
  TSS<-sum(res.lm$residuals^2)
  R2.outcome.PMBC<-1-RSS/TSS 
  
  R2.gene.PMBC<-c()
  for(i in 1:ncol(X)){
    data1<-data.frame(gene=X[,i],cluster=as.factor(cluster))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.PMBC[i]<-1-RSS/TSS 
  }
  R2.gene.PMBC<-R2.gene.PMBC[select.feature]
  
  mult_pdf=t(apply(X.test1,1,function(x) {
    mult_pdf_v<-rep(NA, K)
    for(i in 1:K){
      mult_pdf_v[i]<-as.numeric(mult_density(x,mu=mu_est_train[,i],sigma=sigma_est_train))
    }
    return(mult_pdf_v)
  }))
  
  d<-apply(mult_pdf,1,function(x) brob(x)*pi_est_train)
  postprob.matrix<-sapply(d,function(x) as.numeric(x/sum(x)))
  
  grp_assign_test=apply(postprob.matrix,2,which.max)
  cluster.test<-grp_assign_test
  
  cluster.test1<-cluster.test
  index<-order(c(median(Y[cluster.train1==1]),median(Y[cluster.train1==2]),
                 median(Y[cluster.train1==3])),decreasing = F)
  cluster.test[which(cluster.test1==index[1])]<-1
  cluster.test[which(cluster.test1==index[2])]<-2
  cluster.test[which(cluster.test1==index[3])]<-3
  
  
  cluster<-cluster.test
  cluster<-as.numeric(as.character(cluster))
  cluster[which(cluster.test==1)]<-"low"
  cluster[which(cluster.test==2)]<-"medium"
  cluster[which(cluster.test==3)]<-"high"
  data<-data.frame(Y=Y.test1,cluster=cluster)
  data$cluster<-factor(data$cluster,levels = c("low","medium","high"))
  
  p11<-ggplot(data)+aes(x=cluster,y=Y,fill=cluster)+geom_boxplot()+theme_bw()+
    scale_fill_manual(values =  c("#61D04F", "#2297E6", "#28E2E5"))+ylab("Fev1")+
    theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"),
          strip.text.x = element_text(size = 40),legend.text = element_text(size=25),
          legend.title = element_text(size=25),
          legend.key.size = unit(1.5, 'cm'),legend.position="none",panel.border = element_rect(colour = "black",size=2))
  
  mod1<-aov(Y~cluster,data=data)
  mod1<-summary(mod1)
  pvalue.test.mbc<-mod1[[1]]$`Pr(>F)`[1]
  
  pred.test.mbc<-rep(NA,length(Y.test1))
  for(i in 1:K){
    data.temp<-data.frame(Y=Y[which(grp_assign_train==i)],Z[which(grp_assign_train==i),])
    mod.lm<-lm(Y~.,data=data.temp)
    new.data<-as.data.frame(as_tibble(Z.test1)[which(grp_assign_test==i),])
    pred.test.mbc[which(grp_assign_test==i)]<-predict(mod.lm,newdata = new.data)
    
  }
  
  RMSE.test1.PMBC<-sqrt(sum((pred.test.mbc-Y.test1)^2)/length(Y.test1)) 
  
  data1<-data.frame(Y=Y.test1,cluster=as.factor(cluster),x1=Z.test1[,1],x2=Z.test1[,2])
  res.lm <- lm(Y~ cluster+x1+x2, data = data1)
  RSS<-sum(res.lm$residuals^2)
  res.lm <- lm(Y ~ 1+x1+x2, data = data1)
  TSS<-sum(res.lm$residuals^2)
  R2.outcome.test1.PMBC<-1-RSS/TSS
  
  R2.gene.test1.PMBC<-c()
  for(i in 1:ncol(X.test1)){
    data1<-data.frame(gene=X.test1[,i],cluster=as.factor(cluster))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.test1.PMBC[i]<-1-RSS/TSS 
  }
  R2.gene.test1.PMBC<-R2.gene.test1.PMBC[select.feature]
  
  p10 <- annotate_figure(p10,top = text_grob("(II.E1) outcome separation (discovery)",face = "bold", size = 20 ))
  p11 <- annotate_figure(p11,top = text_grob("(II.E2) outcome separation (validation)",face = "bold", size = 20 ))
  print("finish PMBC plot")
  
  
  
  res.data=data.frame(R2.outcome=c(R2.outcome.wcs,R2.outcome.GM,R2.outcome.SKM,R2.outcome.PMBC),
                      RMSE=c(RMSE.wcs,RMSE.GM,RMSE.SKM,RMSE.PMBC),
                      R2.gene=c(mean(R2.gene.wcs),mean(R2.gene.GM),mean(R2.gene.SKM),mean(R2.gene.PMBC)),
                      R2.outcome.test1=c(R2.outcome.test1.wcs,R2.outcome.test1.GM,R2.outcome.test1.SKM,R2.outcome.test1.PMBC),
                      RMSE.test1=c(RMSE.test1.wcs,RMSE.test1.GM,RMSE.test1.SKM,RMSE.test1.PMBC),
                      R2.gene.test1=c(mean(R2.gene.test1.wcs),mean(R2.gene.test1.GM),mean(R2.gene.test1.SKM),mean(R2.gene.test1.PMBC)),
                      ngenes=c(num.yujia.num,num.Peng.num,num.SKM.num,num.PMBC.num))
  rownames(res.data)=c("ogClustWJL","ogClust","SKM","PMBC")
  return(list(p1,p2,p4,p5,p7,p8,p10,p11,res.data))
}

plotlist=readIN.EXT(ngenes=300)
plotlist[[9]]
plot.WJL <- ggarrange(plotlist[[1]], plotlist[[2]], nrow=2)
plot.GM <- ggarrange(plotlist[[3]],plotlist[[4]],nrow=2)
plot.SKM<-ggarrange(plotlist[[5]], plotlist[[6]], nrow=2)
plot.PMBC<-ggarrange(plotlist[[7]], plotlist[[8]],nrow=2)

plot_combined<-grid.arrange(arrangeGrob(plot.WJL,
                                        plot.GM,
                                        plot.SKM,
                                        plot.PMBC,
                                        ncol=4),heights=c(10, 1))

file.name<-"FigureIII.png"
ggsave(plot_combined,filename = file.name, width = 70, height = 40, units = "cm")

plotlist=readIN.EXT(ngenes=400)
plotlist[[9]]
plot.WJL <- ggarrange(plotlist[[1]], plotlist[[2]], nrow=2)
plot.GM <- ggarrange(plotlist[[3]],plotlist[[4]],nrow=2)
plot.SKM<-ggarrange(plotlist[[5]], plotlist[[6]], nrow=2)
plot.PMBC<-ggarrange(plotlist[[7]], plotlist[[8]],nrow=2)

plot_combined<-grid.arrange(arrangeGrob(plot.WJL,
                                        plot.GM,
                                        plot.SKM,
                                        plot.PMBC,
                                        ncol=4),heights=c(10, 1))

file.name<-"SupplementFigureVI_II.png"
ggsave(plot_combined,filename = file.name, width = 70, height = 40, units = "cm")

readIN.EXTtune=function() {
  load("Data/GSE47460_1_GSE47460_2_GSE37147_GSE42057.Rdata")
  g.mean<-apply(X.GSE47460_1,1,mean)
  cut.mean=quantile(g.mean,probs=0.5)
  X.GSE47460_1=X.GSE47460_1[g.mean>cut.mean,] # remove 50% lowest mean expression genes
  
  g.sd=apply(X.GSE47460_1,1, sd) #cv in original scale
  cut.sd=quantile(g.sd,probs=0.5)
  X.GSE47460_1=X.GSE47460_1[g.sd>=cut.sd,] # further remove 50% lowest variance genes
  dim(X.GSE47460_1) #3240  331
  
  X.GSE47460_1<-t(X.GSE47460_1)
  X.GSE47460_1<-scale(X.GSE47460_1)
  
  Z.GSE47460_1<-as.matrix(Z.GSE47460_1)
  Z.GSE47460_1[,1]<-(Z.GSE47460_1[,1]-mean(Z.GSE47460_1[,1]))/sd(Z.GSE47460_1[,1])
  Y1<-scale(Y.GSE47460_1)
  
  load("mod.WJL.tune.RData")
  load("mod.all.tune.RData")
  
  X.GSE47460_2=t(X.GSE47460_2)
  X.GSE47460_2=X.GSE47460_2[,match(colnames(X.GSE47460_1),colnames(X.GSE47460_2))] #dim 87 1620
  X.GSE47460_2<-scale(X.GSE47460_2)
  Z.GSE47460_2<-as.matrix(Z.GSE47460_2)
  Z.GSE47460_2[,1]<-(Z.GSE47460_2[,1]-mean(Z.GSE47460_2[,1]))/sd(Z.GSE47460_2[,1])
  Y2<-scale(Y.GSE47460_2)
  
  plotlist=EXTplot.gen(X=X.GSE47460_1,Y=Y1,Z=Z.GSE47460_1,mod.WJL=mod.WJL,modfit=modfit,ngenes=ngenes,K=3,
                       X.test1=X.GSE47460_2,Z.test1=Z.GSE47460_2,Y.test1=Y2)
  return(plotlist)
}

plotlist=readIN.EXTtune()
plotlist[[9]]
plot.WJL <- ggarrange(plotlist[[1]], plotlist[[2]], nrow=2)
plot.GM <- ggarrange(plotlist[[3]],plotlist[[4]],nrow=2)
plot.SKM<-ggarrange(plotlist[[5]], plotlist[[6]], nrow=2)
plot.PMBC<-ggarrange(plotlist[[7]], plotlist[[8]],nrow=2)
plot_combined<-grid.arrange(arrangeGrob(plot.WJL,
                                        plot.GM,
                                        plot.SKM,
                                        plot.PMBC,
                                        ncol=4),heights=c(10, 1))

file.name<-"SupplementFigureVI_I.png"
ggsave(plot_combined,filename = file.name, width = 70, height = 40, units = "cm")

