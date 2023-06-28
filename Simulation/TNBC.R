rm(list=ls())
library(sparcl)
library(adaHuber)
library(parallel)
library(Brobdingnag)
library(glmnet)
library(survival)
library(dplyr)

source("Functions.R")

load("Data/Data.Metabric.after.Combat_TNBC.Rdata")
X.train<-Data.Metabric.after.Combat$Expression #dim 275 18964
Y.train<-Data.Metabric.after.Combat$OS
delta.train<-Data.Metabric.after.Combat$OS.event
Z.train<-Data.Metabric.after.Combat$covariate
Z.train<-as.matrix(Z.train[,1])

Index.Mean<-order(apply(X.train,2,mean),decreasing = T)
X.train<-X.train[,Index.Mean[1:(ncol(X.train)/2)]] #dim 275 9482
Index.Sd<-order(apply(X.train,2,sd),decreasing = T)
X.train<-X.train[,Index.Sd[1:(ncol(X.train)/2)]]
X.train<-scale(X.train) #dim 275 4741

load("Data/Data.ScanB.after.Combat_TNBC.Rdata")
X.test<-Data.ScanB.after.Combat$Expression
Y.test<-Data.ScanB.after.Combat$OS
delta.test<-Data.ScanB.after.Combat$OS.event
Z.test<-Data.ScanB.after.Combat$covariate
Z.test<-as.matrix(Z.test$Age)
X.test<-scale(X.test)

index<-match(colnames(X.train),colnames(X.test))
X.test<-X.test[,index]
index<-which(apply(X.test,2,function(x){sum(is.na(x))})!=0)
X.test[,index]<-0


set.seed(12315)
mod.kmeans<-KMeansSparseCluster(X.train,K = 2,nstart = 50)
center<-matrix(NA,nrow=ncol(X.train),ncol=2)
center[,1]<-apply(X.train[which(mod.kmeans[[20]]$Cs==1),],2,mean)
center[,2]<-apply(X.train[which(mod.kmeans[[20]]$Cs==2),],2,mean)
c_center<-center
s_G<-300
#s_G<-400
K<-2


set.seed(12315)
index.sample<-sample(rep(1:10,length.out=nrow(X.train)))
w_vector<-seq(0,1,0.1)
wcs.lambda.surv<-mclapply(1:length(w_vector),function(i){
  print(i)
  w<-w_vector[i]
  w<-(s_G*w)/(s_G*w+1-w)
  res.lambda<-region.lambda.surv(lambda1 =45,lambda2=0,iteration = 20,Y=Y.train,G=X.train,X=Z.train,center=center,w=w,s_G=s_G,K=K,delta.train = delta.train )
  return(res.lambda)
},mc.cores=11)


res1<-mclapply(1:length(w_vector),function(ind_w){
  w<-w_vector[ind_w]
  w<-(s_G*w)/(s_G*w+1-w)
  lambda.vector<-wcs.lambda.surv[[ind_w]]$lambda
  temp.F<-c()
  for(ind_lambda in 1:length(lambda.vector)){
    lambda<-lambda.vector[ind_lambda]
    mod.all<-ogClust_Surv(x=Z.train,G=t(X.train),y=Y.train,y.ind=delta.train,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
    #sum(apply(mod.all$result_list$mu,1,function(x){length(unique(x))})!=1)
    cluster1<-predict_test(mod.all,K=K,D.test=as.matrix(X.train),X1=as.matrix(Z.train),p=ncol(X.train))
    cluster1<-cluster1$clus
    centers<-cbind(apply(X.train[which(cluster1==1),],2,mean),
                   apply(X.train[which(cluster1==2),],2,mean)
    )
    cluster<-c()
    num1<-c()
    for(ind_cv in 1:10){
      ind_sample<-which(index.sample==ind_cv)
      G.test1<-X.train[ind_sample,]
      X.test1<-as.matrix(Z.train[ind_sample,])
      Y.test1<-Y.train[ind_sample]
      delta.test1<-delta.train[ind_sample]
      G.train1<-X.train[(1:nrow(X.train))[-ind_sample],]
      X.train1<-as.matrix(Z.train[(1:nrow(X.train))[-ind_sample],])
      Y.train1<-Y.train[(1:nrow(X.train))[-ind_sample]]
      delta.train1<-delta.train[(1:nrow(X.train))[-ind_sample]]



      mod<-ogClust_Surv(x=X.train1,G=t(G.train1),y=Y.train1,y.ind=delta.train1,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
      num1[ind_cv]<-sum(apply(mod$result_list$mu,1,function(x){length(unique(x))})!=1)
      res1<-predict_test(mod,K=K,D.test=as.matrix(G.test1),X1=as.matrix(X.test1),p=ncol(G.test1))
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
    if(length(table(cluster))>1){
      data1<-data.frame(cluster=cluster,Y=Y.train,Y.ind=delta.train,X=Z.train)
      mod.cox<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+X, data = data1)
      mod.cox1<-summary(mod.cox)
      Fstat<-abs(mod.cox1$waldtest[1])
      temp.F[ind_lambda]<-Fstat
    }else{
      temp.F[ind_lambda]<-NA
    }

  }

  return(temp.F)
},mc.cores = 11)

max.F<-sapply(res1,function(x){max(x,na.rm=T)})
w0=w_vector[which.max(max.F)]
w<-(s_G*w0)/(s_G*w0+1-w0)
lambda_vector1<-seq(5,50,length.out = 22)

set.seed(12315)
res1=mclapply(1:length(lambda_vector1), function(ind_lambda) {
  lambda<-lambda_vector1[ind_lambda]
  mod.all<-ogClust_Surv(x=Z.train,G=t(X.train),y=Y.train,y.ind=delta.train,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
  select.feature<-as.numeric(apply(mod.all$result_list$mu,1,function(x){length(unique(x))})!=1)
  num.all<-sum(select.feature==1)
  if (num.all<100) {
    print(paste0(ind_lambda,"th lambda selects a small number of genes is ",num.all))
    final.res<-list(Fstat.outcome=NA,Fstat.gene=NA,num=0)
    return(final.res)
  }
  cluster1<-predict_test(mod.all,K=K,D.test=as.matrix(X.train),X1=as.matrix(Z.train),p=ncol(X.train))
  cluster1<-cluster1$clus
  centers<-cbind(apply(X.train[which(cluster1==1),],2,mean),
                 apply(X.train[which(cluster1==2),],2,mean)
  )
  cluster<-c()
  num1<-c()
  for(ind_cv in 1:10){
    ind_sample<-which(index.sample==ind_cv)
    G.test1<-X.train[ind_sample,]
    X.test1<-as.matrix(Z.train[ind_sample,])
    Y.test1<-Y.train[ind_sample]
    delta.test1<-delta.train[ind_sample]
    G.train1<-X.train[(1:nrow(X.train))[-ind_sample],]
    X.train1<-as.matrix(Z.train[(1:nrow(X.train))[-ind_sample],])
    Y.train1<-Y.train[(1:nrow(X.train))[-ind_sample]]
    delta.train1<-delta.train[(1:nrow(X.train))[-ind_sample]]


    mod<-ogClust_Surv(x=X.train1,G=t(G.train1),y=Y.train1,y.ind=delta.train1,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
    num1[ind_cv]<-sum(apply(mod$result_list$mu,1,function(x){length(unique(x))})!=1)
    res1<-predict_test(mod,K=K,D.test=as.matrix(G.test1),X1=as.matrix(X.test1),p=ncol(G.test1))
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
  #print("10 fold CV finish!")
  data1<-data.frame(cluster=cluster,Y=Y.train,Y.ind=delta.train,X=Z.train)
  mod.cox<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+X, data = data1)
  mod.cox1<-summary(mod.cox)
  Fstat.outcome<-abs(mod.cox1$waldtest[1])
  #res<-Rsquare(cluster = cluster,Y = Y.train,X = Z.train,G = t(X.train))
  Fstat.gene<-c()
  for(i in 1:ncol(X.train)){
    data1<-data.frame(gene=X.train[,i],cluster=as.factor(cluster))
    mod<-lm(gene~cluster,data=data1)
    mod<-summary(mod)
    Fstat.gene[i]<-abs(mod$coefficients[2,3])
  }
  Fstat.gene<-mean(Fstat.gene[which(select.feature==1)])
  final.res<-list(Fstat.outcome=Fstat.outcome,Fstat.gene=Fstat.gene,num=num.all)
  return(final.res)
},mc.cores=22)

index<-which(unlist(lapply(res1,is.list)))
res1<-res1[index]
Fstat.outcome<-unlist(lapply(res1,function(x){x$Fstat.outcome}))
Fstat.gene<-unlist(lapply(res1,function(x){x$Fstat.gene}))
lambda=lambda_vector1[which.max(sqrt(Fstat.gene*Fstat.outcome))]

mod<-ogClust_Surv(x=Z.train,G=t(X.train),y=Y.train,y.ind=delta.train,c_center=c_center,
                      lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=100,w_outcome=w,w_G=1-w)

res.test<-predict_test(mod,K=2,D.test=as.matrix(X.test),X1=Z.test,p=ncol(X.train))
mod.WJL=list(mod=mod,res.test=res.test)
save(mod.WJL,file="survmod.WJL.tune.RData")

##########GM######################

n1=nrow(X.train) # number of samples
NG=ncol(X.train) # number of genes
np=ncol(Z.train) # number of covariates

set.seed(12315)
index.sample<-sample(rep(1:10,length.out=nrow(X.train)))
mod.kmeans<-KMeansSparseCluster(X.train,K = 2,nstart = 10)
cluster<-mod.kmeans[[20]]$Cs

index1<-which(cluster==1)
data1<-data.frame(y=Y.train[index1],x1=Z.train[index1,],delta=delta.train[index1])
mod1<-survreg(Surv(y, delta) ~ ., data1,robust=T,
              dist="loglogistic")

index2<-which(cluster==2)
data2<-data.frame(y=Y.train[index2],x1=Z.train[index2,],delta=delta.train[index2])
mod2<-survreg(Surv(y, delta) ~ ., data2,robust=T,
              dist="loglogistic")


beta0_int<-c(mod1$coefficients[1],mod2$coefficients[1])
data<-data.frame(y=Y.train,x=Z.train,delta=delta.train)
mod<-survreg(Surv(y, delta) ~ ., data,robust=T,
             dist="loglogistic")
beta_int<-mod$coefficients[2]
#beta_int = runif(np, 0, 3)
K<-2

mod<-glmnet(X.train,as.factor(cluster),family = "binomial",lambda=0)
mod1<-coef(mod)
gamma_int<-c(as.numeric(mod1))


sigma2_int<-1
theta_int = c(beta_int, gamma_int, beta0_int, sigma2_int)
#lambda_vector_Peng<-c(0.001,0.005,0.008,seq(0.01,0.5,0.02))
lambda_vector_Peng<-seq(0.0001,0.015,0.0005)
wcs.GM<-mclapply(1:length(lambda_vector_Peng),function(ind_lambda){
  X.train<-scale(X.train)
  lambda<-lambda_vector_Peng[ind_lambda]
  X.train<-as.data.frame(X.train)
  mod.all<-fit.ogClust.surv(n=n1, K=K, np=np, NG=NG, lambda=lambda,delta = delta.train,
                            alpha=1, G=X.train, Y=Y.train, X=as.data.frame(Z.train), theta_int=theta_int)

  beta=mod.all$par$gamma
  beta<-as.data.frame(beta)
  rownames(beta)[1]<-"intercept"
  rownames(beta)[2:nrow(beta)]<-colnames(X.train)
  select.feature<-as.numeric(apply(beta,1,function(x){sum(x==0)})!=K)[-1] #genes with non-zero weight of at least one clsuter
  #num.vector[ind_lambda]<-sum(select.feature==1)
  num.all<-sum(select.feature==1)
  if (num.all<50) {
    print(paste0(ind_lambda,"th lambda selects a small number of genes is ",num.all))
    final.res<-list(Fstat.outcome=NA,Fstat.gene=NA,num=0)
    return(final.res)
  }
  cluster1<-mod.all$grp_assign
  if(length(table(cluster1))!=2){
    print(paste0(ind_lambda,"th lambda has less than 2 clusters!"))
    final.res<-list(Fstat.outcome=NA,Fstat.gene=NA,num=num.all)
    return(final.res)
  }
  centers<-cbind(apply(X.train[which(cluster1==1),],2,mean),
                 apply(X.train[which(cluster1==2),],2,mean))
  cluster<-c()
  #num1<-c()
  for(ind_cv in 1:10){
    ind_sample<-which(index.sample==ind_cv)
    G.test1<-X.train[ind_sample,]
    X.test1<-as.matrix(Z.train[ind_sample,])
    Y.test1<-Y.train[ind_sample]
    delta.test1<-delta.train[ind_sample]
    G.train1<-X.train[(1:nrow(X.train))[-ind_sample],]
    X.train1<-as.matrix(Z.train[(1:nrow(X.train))[-ind_sample],])
    Y.train1<-Y.train[(1:nrow(X.train))[-ind_sample]]
    delta.train1<-delta.train[(1:nrow(X.train))[-ind_sample]]

    n1sub=nrow(G.train1)
    NGsub=ncol(G.train1)
    np=ncol(X.train1)
    mod<-fit.ogClust.surv(n=n1sub, K=K, np=np, NG=NGsub, lambda=lambda,delta = delta.train1,
                          alpha=1, G=G.train1, Y=Y.train1, X=X.train1, theta_int=theta_int)


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
  #print("10 fold CV finish!")
  data1<-data.frame(cluster=cluster,Y=Y.train,Y.ind=delta.train,X=Z.train)
  mod.cox<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+X, data = data1)
  mod.cox1<-summary(mod.cox)
  Fstat.outcome<-abs(mod.cox1$waldtest[1])
  #res<-Rsquare(cluster = cluster,Y = Y.train,X = Z.train,G = t(X.train))
  Fstat.gene<-c()
  for(i in 1:ncol(X.train)){
    data1<-data.frame(gene=X.train[,i],cluster=as.factor(cluster))
    mod<-lm(gene~cluster,data=data1)
    mod<-summary(mod)
    Fstat.gene[i]<-abs(mod$coefficients[2,3])
  }
  Fstat.gene<-mean(Fstat.gene[which(select.feature==1)])
  final.res<-list(Fstat.outcome=Fstat.outcome,Fstat.gene=Fstat.gene,num=num.all)
  return(final.res)
},mc.cores=28)

Fstat.outcome<-unlist(lapply(wcs.GM,function(x){x$Fstat.outcome}))
Fstat.gene<-unlist(lapply(wcs.GM,function(x){x$Fstat.gene}))
num<-unlist(lapply(wcs.GM,function(x){x$num}))
num[which.max(sqrt(Fstat.gene*Fstat.outcome))]
lambda=lambda_vector_Peng[which.max(sqrt(Fstat.gene*Fstat.outcome))]

mod.GM<-fit.ogClust.surv(n=n1, K=K, np=np, NG=NG, lambda=lambda,delta = delta.train,
                                 alpha=1, G=as.data.frame(X.train), Y=Y.train, X=as.data.frame(Z.train), theta_int=theta_int)
########## SKM ######################
K=2
set.seed(12315)
optimal_wbounds<-KMeansSparseCluster.permute(X.train,K=K,wbounds=30:60,silent = T)$bestw
mod.SKM<-KMeansSparseCluster(X.train,K=K,wbounds=optimal_wbounds,silent = T)

########PMBC###########
K=2
set.seed(12315)
data.kmeans<-kmeans(as.matrix(X.train),K)
c_size<-data.kmeans$size
n=nrow(X.train)
pi_int<-c_size/n
miu_int<-sapply(1:K,function(x) apply(as.matrix(X.train[data.kmeans$cluster==x,]),2,mean))
sigma_c<-sapply(1:K,function(x) apply(as.matrix(X.train[data.kmeans$cluster==x,]),2,var))
sum_squares<-sapply(1:K,function(x) c_size[x]*sigma_c[,x])
sigma_int<-apply(sum_squares,1,sum)/n

lambda.all=seq(30,60,length.out = 22)

wcs.MBC=mclapply(1:length(lambda.all),function(ind_lambda) {
  lambda=lambda.all[ind_lambda]
  BIC<-tryCatch({
    BIC<-em_mbc(t(X.train), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                lambda=lambda, K=K, no_init=1, max_iter = 200)$BIC
  }, error=function(e){
    BIC=Inf
    return(BIC)
  })
  return(BIC)
},mc.cores=22)

lambda=lambda.all[which.min(unlist(wcs.MBC))]

mod.PMBC<-em_mbc(t(X.train), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                lambda=lambda, K=K, no_init=1, max_iter = 200)

modfit=list(mod.GM=mod.GM, mod.SKM=mod.SKM,mod.PMBC=mod.PMBC)
save(modfit,file="survmod.all.tune.RData")

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

######Fixed number of genes for GM,PMBC, SKM ########
fixgenes.survfit=function(X=X.train,Y=Y.train,Z=Z.train,delta.train=delta.train,K=2,lambda.GM,wbound,lambda.PMBC) {
  set.seed(12315)
  start_time=Sys.time()
  mod.kmeans<-KMeansSparseCluster(X,K = K,nstart = 50)
  cluster<-mod.kmeans[[20]]$Cs

  index1<-which(cluster==1)
  data1<-data.frame(y=Y[index1],x1=Z[index1,],delta=delta.train[index1])
  mod1<-survreg(Surv(y, delta) ~ ., data1,robust=T,
                dist="loglogistic")

  index2<-which(cluster==2)
  data2<-data.frame(y=Y[index2],x1=Z[index2,],delta=delta.train[index2])
  mod2<-survreg(Surv(y, delta) ~ ., data2,robust=T,
                dist="loglogistic")


  beta0_int<-c(mod1$coefficients[1],mod2$coefficients[1])
  data<-data.frame(y=Y,x=Z,delta=delta.train)
  mod<-survreg(Surv(y, delta) ~ ., data,robust=T,
               dist="loglogistic")
  beta_int<-mod$coefficients[2]

  mod<-glmnet(X,as.factor(cluster),family = "binomial",lambda=0)
  mod1<-coef(mod)
  gamma_int<-c(as.numeric(mod1))
  sigma2_int<-1
  theta_int = c(beta_int, gamma_int, beta0_int, sigma2_int)

  n1=nrow(X)
  NG=ncol(X)
  np=ncol(Z)
  mod.ogClust<-fit.ogClust.surv(n=n1, K=K, np=np, NG=NG, lambda=lambda.GM,delta = delta.train,
                                alpha=1, G=as.data.frame(X), Y=Y, X=as.data.frame(Z), theta_int=theta_int)

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
# lambda<-0.007  select about 300 genes
# lambda<- 0.0015  select about 400 genes

# For SKM model:
# wbound=11.9 select about 300 genes
# wbound=14 select about 400 genes

# For PMBC model:
# lambda=65.8 select about 300 genes
# lambda=62.7 select about 400 genes

start_time <- Sys.time()
modfit=fixgenes.survfit(X=X.train,Y=Y.train,Z=Z.train,delta.train=delta.train,K=2,lambda.GM=0.007,wbound=11.9,lambda.PMBC=65.8)
end_time <- Sys.time()
save(modfit,file="survmod.all.fix400.Rdata")
save(modfit,file="survmod.all.fix300.Rdata")
#########Generate figures##########
readIN.survEXT=function(ngenes=300) {
  load("~/OneDrive - University of Pittsburgh/AOASyujia_peng/ogclust_Wei_Wenjia/TNBC/Data.Metabric.after.Combat_TNBC.Rdata")
  X.train<-Data.Metabric.after.Combat$Expression #dim 275 18964
  Y.train<-Data.Metabric.after.Combat$OS
  delta.train<-Data.Metabric.after.Combat$OS.event
  Z.train<-Data.Metabric.after.Combat$covariate
  Z.train<-as.matrix(Z.train[,1])

  Index.Mean<-order(apply(X.train,2,mean),decreasing = T)
  X.train<-X.train[,Index.Mean[1:(ncol(X.train)/2)]] #dim 275 9482
  Index.Sd<-order(apply(X.train,2,sd),decreasing = T)
  X.train<-X.train[,Index.Sd[1:(ncol(X.train)/2)]]
  X.train<-scale(X.train) #dim 275 4741
  load(paste0("survmod.WJL.fix",ngenes,".RData"))
  load(paste0("survmod.all.fix",ngenes,".RData"))

  load("~/OneDrive - University of Pittsburgh/AOASyujia_peng/ogclust_Wei_Wenjia/TNBC/Data.ScanB.after.Combat_TNBC.Rdata")
  X.test<-Data.ScanB.after.Combat$Expression
  Y.test<-Data.ScanB.after.Combat$OS
  delta.test<-Data.ScanB.after.Combat$OS.event
  Z.test<-Data.ScanB.after.Combat$covariate
  Z.test<-as.matrix(Z.test$Age)
  X.test<-scale(X.test)

  index<-match(colnames(X.train),colnames(X.test))
  X.test<-X.test[,index]
  index<-which(apply(X.test,2,function(x){sum(is.na(x))})!=0)
  X.test[,index]<-0

  plotlist=survEXTplot.gen(X=X.train,Y=Y.train,Z=Z.train,delta.train=delta.train,mod.WJL=mod.WJL,modfit=modfit,ngenes=ngenes,K=2,
                           X.test=X.test,Z.test=Z.test,Y.test=Y.test,delta.test=delta.test)
  return(plotlist)
}
library(pheatmap)
library(gplots)
library(tibble)
library(ggplotify)
library(ggpubr)
library(mclust)
library(survminer)
library(survcomp)
library(survival)

survEXTplot.gen=function(X=X.train,Y=Y.train,Z=Z.train,delta.train=delta.train,mod.WJL=mod.WJL,modfit=modfit,ngenes=400,K=2,
                         X.test=X.test,Z.test=Z.test,Y.test=Y.test,delta.test=delta.test) {

  mod=mod.WJL$mod
  res.test=mod.WJL$res.test
  select.feature<-as.numeric(apply(mod$result_list$mu,1,function(x){length(unique(x))})>1)

  num.yujia.num<-sum(select.feature==1)
  cluster.train = factor(apply(mod$result_list$z,1,which.max))

  data1<-data.frame(cluster=cluster.train,Y=Y,Y.ind=delta.train,X=Z)
  mod.cox<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+X, data = data1)
  mod.cox1<-summary(mod.cox)
  coef<-mod.cox1$coefficients[1,1]

  index1<-which(cluster.train==1)
  index2<-which(cluster.train==2)
  if(coef>0){
    cluster.train[index2]<-1
    cluster.train[index1]<-2
  }


  cluster<-cluster.train
  cluster<-as.numeric(as.character(cluster))

  data<-data.frame(Y=Y,cluster=cluster,Y.ind=delta.train)
  surv_object <- Surv(time = data$Y,
                      event = data$Y.ind)

  fit.os<-surv_fit(surv_object ~ data$cluster,
                   data = data)

  p1<-ggsurvplot(fit.os, data = data,legend.title="Clusters",size=2,censor=FALSE
                 ,palette=1:2,title = "",risk.table=F,
                 pval=F,conf.int = F,xlab = "Time (Years)",legend="none",ylab="P(Surv)")%++%theme_bw()
  p1<-ggpar(p1,
            font.main = c(40, "bold"),
            font.x = c(40, "bold"),
            font.y = c(40, "bold"),
            font.caption = c(40, "bold"),
            font.legend = c(40, "bold"),
            font.tickslab = c(40, "bold"))
  p1<-p1%++%theme(legend.position="none",panel.border = element_rect(colour = "black",size=2))
  p1<-as.ggplot(p1$plot)


  data1<-data.frame(group=cluster.train,Y.ind=delta.train,Y=Y,x=Z[,1])
  mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(group)+x, data = data1)
  mod.cox2<-coxph(Surv(Y, Y.ind) ~ x, data = data1)
  mod.anova<-anova(mod.cox1,mod.cox2)
  F.WJL<-mod.anova$Chisq[2]
  p.WJL=mod.anova$`Pr(>|Chi|)`[2]
  intercept<-mod$result_list$int_coef[apply(mod$result_list$z,1,which.max)]
  pred<-mod$result_list$int_coef[3]*Z+intercept

  cindex.WJL<-concordance.index(x=1/as.numeric(pred), surv.time=Y, surv.event=delta.train, method="noether")$c.index
  cindex.adjust<-c()
  for(i in 1:500){
    cindex.adjust[i]<-concordance.index(x=1/as.numeric(sample(pred)), surv.time=Y, surv.event=delta.train, method="noether")$c.index
  }
  cindex.WJL<-(cindex.WJL-mean(cindex.adjust))/(1-mean(cindex.adjust))

  R2.gene.wcs<-c()
  for(i in 1:ncol(X)){
    data1<-data.frame(gene=X[,i],cluster=as.factor(cluster.train))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.wcs[i]<-1-RSS/TSS
  }
  R2.gene.wcs<-R2.gene.wcs[which(select.feature==1)]

  cluster.test<-res.test$clus
  index1<-which(cluster.test==1)
  index2<-which(cluster.test==2)
  if(coef>0){
    cluster.test[index2]<-1
    cluster.test[index1]<-2
  }

  cluster<-cluster.test
  cluster<-as.numeric(as.character(cluster))

  data<-data.frame(Y=Y.test,cluster=cluster,Y.ind=delta.test)
  #data$cluster<-factor(data$cluster,levels = c("low","high"))
  surv_object <- Surv(time = data$Y,
                      event = data$Y.ind)
  fit.os<-surv_fit(surv_object ~ data$cluster,
                   data = data)

  p2<-ggsurvplot(fit.os, data = data,legend.title="Clusters",size=2,censor=FALSE
                 ,palette=1:2,title = "",risk.table=F,
                 pval=F,conf.int = F,xlab = "Time (Years)",legend="none",ylab="P(Surv)")%++%theme_bw()
  p2<-ggpar(p2,
            font.main = c(40, "bold"),
            font.x = c(40, "bold"),
            font.y = c(40, "bold"),
            font.caption = c(40, "bold"),
            font.legend = c(40, "bold"),
            font.tickslab = c(40, "bold"))

  p2<-p2%++%theme(legend.position="none",panel.border = element_rect(colour = "black",size=2))
  p2<-as.ggplot(p2$plot)

  data1<-data.frame(group=cluster.test,Y.ind=delta.test,Y=Y.test,x=Z.test)
  mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(group)+x, data = data1)
  mod.cox2<-coxph(Surv(Y, Y.ind) ~ x, data = data1)
  mod.anova<-anova(mod.cox1,mod.cox2)
  F.test.WJL<-mod.anova$Chisq[2]
  p.test.WJL<-mod.anova$`Pr(>|Chi|)`[2]
  pred<-res.test$Y.hard

  cindex.WJL.test<-concordance.index(x=1/as.numeric(pred), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index #0.7437702 without SC
  cindex.adjust<-c()
  for(i in 1:500){
    cindex.adjust[i]<-concordance.index(x=1/as.numeric(sample(pred)), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index
  }
  cindex.WJL.test<-(cindex.WJL.test-mean(cindex.adjust))/(1-mean(cindex.adjust)) #0.4840802 for SC, 0.4851878 without SC

  R2.gene.test.wcs<-c()
  for(i in 1:ncol(X.test)){
    data1<-data.frame(gene=X.test[,i],cluster=as.factor(cluster.test))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.test.wcs[i]<-1-RSS/TSS
  }
  R2.gene.test.wcs<-R2.gene.test.wcs[which(select.feature==1)]
  p1 <- annotate_figure(p1,top = text_grob("(I.B1) outcome separation (discovery)",face = "bold", size = 20 ))
  p2 <- annotate_figure(p2,top = text_grob("(I.B2) outcome separation (validation)",face = "bold", size = 20 ))
  print("finish WJL plot")

  fit.res.num=modfit$mod.GM
  beta=fit.res.num$par$gamma
  beta<-as.data.frame(beta)
  rownames(beta)[1]<-"intercept"
  rownames(beta)[2:nrow(beta)]<-colnames(X)
  select.feature<-as.numeric(apply(beta,1,function(x){sum(x==0)})!=K)[-1]
  num.Peng.num<-sum(select.feature==1)

  cluster.train<-fit.res.num$grp_assign
  data1<-data.frame(cluster=cluster.train,Y=Y,Y.ind=delta.train,X=Z)
  mod.cox<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+X, data = data1)
  mod.cox1<-summary(mod.cox)
  coef<-mod.cox1$coefficients[1,1]

  index1<-which(cluster.train==1)
  index2<-which(cluster.train==2)
  if(coef>0){
    cluster.train[index2]<-1
    cluster.train[index1]<-2
  }

  cluster<-cluster.train
  cluster<-as.numeric(as.character(cluster))

  data<-data.frame(Y=Y,cluster=cluster,Y.ind=delta.train)
  surv_object <- Surv(time = data$Y,
                      event = data$Y.ind)
  fit.os<-surv_fit(surv_object ~ data$cluster,
                   data = data)

  p3<-ggsurvplot(fit.os, data = data,legend.title="Clusters",size=2,censor=FALSE
                 ,palette=1:2,title = "",risk.table=F,
                 pval=F,conf.int = F,xlab = "Time (Years)",legend="none",ylab="P(Surv)")%++%theme_bw()
  p3<-ggpar(p3,
            font.main = c(40, "bold"),
            font.x = c(40, "bold"),
            font.y = c(40, "bold"),
            font.caption = c(40, "bold"),
            font.legend = c(40, "bold"),
            font.tickslab = c(40, "bold"))
  p3<-p3%++%theme(legend.position="none",panel.border = element_rect(colour = "black",size=2))
  p3<-as.ggplot(p3$plot)

  data1<-data.frame(group=cluster.train,Y.ind=delta.train,Y=Y,x=Z[,1])
  mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x, data = data1)
  mod.cox2<-coxph(Surv(Y, Y.ind) ~ x, data = data1)
  mod.anova<-anova(mod.cox1,mod.cox2)
  F.GM<-mod.anova$Chisq[2]
  p.GM<-mod.anova$`Pr(>|Chi|)`[2]

  beta_est = fit.res.num$par$beta
  beta0_est = fit.res.num$par$beta0

  if(coef>0) beta0_est<-beta0_est[c(2,1)]

  intercept<-beta0_est[cluster.train]
  pred<-beta_est*Z+intercept

  cindex.GM<-concordance.index(x=1/as.numeric(pred), surv.time=Y, surv.event=delta.train, method="noether")$c.index
  cindex.adjust<-c()
  for(i in 1:500){
    cindex.adjust[i]<-concordance.index(x=1/as.numeric(sample(pred)), surv.time=Y, surv.event=delta.train, method="noether")$c.index
  }
  cindex.GM<-(cindex.GM-mean(cindex.adjust))/(1-mean(cindex.adjust))

  R2.gene.GM<-c()
  for(i in 1:ncol(X)){
    data1<-data.frame(gene=X[,i],cluster=as.factor(cluster.train))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.GM[i]<-1-RSS/TSS
  }
  R2.gene.GM<-R2.gene.GM[which(select.feature==1)]

  sigma2_est = fit.res.num$par$sigma2
  gamma_est_matrix =fit.res.num$par$gamma

  X.test.add = cbind(1, X.test)
  pai_est.num = sapply(1:K, function(k) exp(X.test.add %*% gamma_est_matrix[,  k, drop = F])/rowSums(exp(X.test.add %*% gamma_est_matrix)))
  assignment.test<-apply(pai_est.num,1,which.max)
  cluster.test<-apply(pai_est.num,1,which.max)

  index1<-which(cluster.test==1)
  index2<-which(cluster.test==2)
  if(coef>0){
    cluster.test[index2]<-1
    cluster.test[index1]<-2
  }

  data<-data.frame(Y=Y.test,cluster=cluster.test,Y.ind=delta.test)
  surv_object <- Surv(time = data$Y,
                      event = data$Y.ind)
  fit.os<-surv_fit(surv_object ~ data$cluster,
                   data = data)

  p4<-ggsurvplot(fit.os, data = data,legend.title="Clusters",size=2,censor=FALSE
                 ,palette=1:2,title = "",risk.table=F,
                 pval=F,conf.int = F,xlab = "Time (Years)",legend="none",ylab="P(Surv)")%++%theme_bw()

  p4<-ggpar(p4,
            font.main = c(40, "bold"),
            font.x = c(40, "bold"),
            font.y = c(40, "bold"),
            font.caption = c(40, "bold"),
            font.legend = c(40, "bold"),
            font.tickslab = c(40, "bold"))

  p4<-p4%++%theme(legend.position="none",panel.border = element_rect(colour = "black",size=2))
  p4<-as.ggplot(p4$plot)

  data1<-data.frame(group=cluster.test,Y.ind=delta.test,Y=Y.test,x=Z.test)
  mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(group)+x, data = data1)
  mod.cox2<-coxph(Surv(Y, Y.ind) ~ x, data = data1)
  mod.anova<-anova(mod.cox1,mod.cox2)
  F.test.GM<-mod.anova$Chisq[2]
  p.test.GM<-mod.anova$`Pr(>|Chi|)`[2]

  intercept<-beta0_est[cluster.test]
  pred<-beta_est*Z.test+intercept
  cindex.GM.test<-concordance.index(x=1/as.numeric(pred), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index
  cindex.adjust<-c()
  for(i in 1:500){
    cindex.adjust[i]<-concordance.index(x=1/as.numeric(sample(pred)), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index
  }
  cindex.GM.test<-(cindex.GM.test-mean(cindex.adjust))/(1-mean(cindex.adjust))

  R2.gene.test.GM<-c()
  for(i in 1:ncol(X.test)){
    data1<-data.frame(gene=X.test[,i],cluster=as.factor(cluster.test))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.test.GM[i]<-1-RSS/TSS
  }
  R2.gene.test.GM<-R2.gene.test.GM[which(select.feature==1)]

  p3 <- annotate_figure(p3,top = text_grob("(I.C1) outcome separation (discovery)",face = "bold", size = 20 ))
  p4 <- annotate_figure(p4,top = text_grob("(I.C2) outcome separation (validation)",face = "bold", size = 20 ))
  print("finish GM plot")


  mod.kmeans=modfit$mod.SKM
  select.feature<-as.numeric(mod.kmeans[[1]]$ws!=0)
  num.SKM.num<-sum(select.feature==1)

  weight<-mod.kmeans[[1]]$ws
  cluster.train<-mod.kmeans[[1]]$Cs
  centers<-matrix(NA,ncol(X),K)
  for(i1 in 1:K){
    if(sum(mod.kmeans[[1]]$Cs==i1)>1){
      centers[,i1]<-apply((X)[which(mod.kmeans[[1]]$Cs==i1),],2,mean)
    }else{
      centers[,i1]<-X[which(mod.kmeans[[1]]$Cs==i1),]
    }

  }

  data1<-data.frame(cluster=cluster.train,Y=Y,Y.ind=delta.train,X=Z)
  mod.cox<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+X, data = data1)
  mod.cox1<-summary(mod.cox)
  coef<-mod.cox1$coefficients[1,1]

  index1<-which(cluster.train==1)
  index2<-which(cluster.train==2)
  if(coef>0){
    cluster.train[index2]<-1
    cluster.train[index1]<-2
  }


  cluster<-cluster.train
  cluster<-as.numeric(as.character(cluster))

  data<-data.frame(Y=Y,cluster=cluster,Y.ind=delta.train)
  surv_object <- Surv(time = data$Y,
                      event = data$Y.ind)
  fit.os<-surv_fit(surv_object ~ data$cluster,
                   data = data)

  p5<-ggsurvplot(fit.os, data = data,legend.title="Clusters",size=2,censor=FALSE
                 ,palette=1:2,title = "",risk.table=F,
                 pval=F,conf.int = F,xlab = "Time (Years)",legend="none",ylab="P(Surv)")%++%theme_bw()
  p5<-ggpar(p5,
            font.main = c(40, "bold"),
            font.x = c(40, "bold"),
            font.y = c(40, "bold"),
            font.caption = c(40, "bold"),
            font.legend = c(40, "bold"),
            font.tickslab = c(40, "bold"))
  p5<-p5%++%theme(legend.position="none",panel.border = element_rect(colour = "black",size=2))
  p5<-as.ggplot(p5$plot)

  data1<-data.frame(group=cluster.train,Y.ind=delta.train,Y=Y,x=Z[,1])
  mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x, data = data1)
  mod.cox2<-coxph(Surv(Y, Y.ind) ~ x, data = data1)
  mod.anova<-anova(mod.cox1,mod.cox2)
  F.SKM<-mod.anova$Chisq[2] #6.247252 for SC, 3.001317 without SC
  p.SKM<-mod.anova$`Pr(>|Chi|)`[2] #0.01243861 for SC, 0.08319688 without SC
  pred<-rep(NA,length(cluster.train))
  for(i in 1:K){
    data.temp<-data.frame(Y=Y[which(cluster.train==i)],x=Z[which(cluster.train==i)],
                          delta=delta.train[which(cluster.train==i)])
    mod.lm<-survreg(Surv(Y,delta)~x,dist="loglogistic",data=data.temp)
    pred[which(cluster.train==i)]<-predict(mod.lm,newdata = data.frame(x=Z[which(cluster.train==i)]))
  }


  cindex.SKM<-concordance.index(x=1/as.numeric(pred), surv.time=Y, surv.event=delta.train, method="noether")$c.index
  cindex.adjust<-c()
  for(i in 1:500){
    cindex.adjust[i]<-concordance.index(x=1/as.numeric(sample(pred)), surv.time=Y, surv.event=delta.train, method="noether")$c.index
  }
  cindex.SKM<-(cindex.SKM-mean(cindex.adjust))/(1-mean(cindex.adjust))

  R2.gene.SKM<-c()
  for(i in 1:ncol(X)){
    data1<-data.frame(gene=X[,i],cluster=as.factor(cluster.train))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.SKM[i]<-1-RSS/TSS
  }
  R2.gene.SKM<-R2.gene.SKM[which(select.feature==1)]

  distance<-apply(X.test,1,function(x){
    distance1<-apply(centers,2,function(y){
      sum(weight*(x-y)^2)
    })
  })

  cluster.test<-apply(distance,2,which.min)
  index1<-which(cluster.test==1)
  index2<-which(cluster.test==2)
  if(coef>0){
    cluster.test[index2]<-1
    cluster.test[index1]<-2
  }


  data<-data.frame(Y=Y.test,cluster=cluster.test,Y.ind=delta.test)
  surv_object <- Surv(time = data$Y,
                      event = data$Y.ind)
  fit.os<-surv_fit(surv_object ~ data$cluster,
                   data = data)

  p6<-ggsurvplot(fit.os, data = data,legend.title="Clusters",size=2,censor=FALSE
                 ,palette=1:2,title = "",risk.table=F,
                 pval=F,conf.int = F,xlab = "Time (Years)",legend="none",ylab="P(Surv)")%++%theme_bw()
  p6<-ggpar(p6,
            font.main = c(40, "bold"),
            font.x = c(40, "bold"),
            font.y = c(40, "bold"),
            font.caption = c(40, "bold"),
            font.legend = c(40, "bold"),
            font.tickslab = c(40, "bold"))
  p6<-p6%++%theme(legend.position="none",panel.border = element_rect(colour = "black",size=2))
  p6<-as.ggplot(p6$plot)

  data1<-data.frame(group=cluster.test,Y.ind=delta.test,Y=Y.test,x=Z.test)
  mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(group)+x, data = data1)
  mod.cox2<-coxph(Surv(Y, Y.ind) ~ x, data = data1)
  mod.anova<-anova(mod.cox1,mod.cox2)
  F.test.SKM<-mod.anova$Chisq[2]
  p.test.SKM<-mod.anova$`Pr(>|Chi|)`[2]
  pred<-rep(NA,length(cluster.test))
  for(i in 1:K){
    data.temp<-data.frame(Y=Y[which(cluster.train==i)],x=Z[which(cluster.train==i)],
                          delta=delta.train[which(cluster.train==i)])
    mod.lm<-survreg(Surv(Y,delta)~x,dist="loglogistic",data=data.temp)
    pred[which(cluster.test==i)]<-predict(mod.lm,newdata = data.frame(x=Z.test[which(cluster.test==i)]))
  }
  cindex.SKM.test<-concordance.index(x=1/as.numeric(pred), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index
  cindex.adjust<-c()
  for(i in 1:500){
    cindex.adjust[i]<-concordance.index(x=1/as.numeric(sample(pred)), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index
  }
  cindex.SKM.test<-(cindex.SKM.test-mean(cindex.adjust))/(1-mean(cindex.adjust))

  R2.gene.test.SKM<-c()
  for(i in 1:ncol(X.test)){
    data1<-data.frame(gene=X.test[,i],cluster=as.factor(cluster.test))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.test.SKM[i]<-1-RSS/TSS
  }
  R2.gene.test.SKM<-R2.gene.test.SKM[which(select.feature==1)]

  p5 <- annotate_figure(p5,top = text_grob("(I.D1) outcome separation (discovery)",face = "bold", size = 20 ))
  p6 <- annotate_figure(p6,top = text_grob("(I.D2) outcome separation (validation)",face = "bold", size = 20 ))
  print("finish SKM plot")


  cltrain=modfit$mod.PMBC
  select.feature<-cltrain$selected.index
  num.PMBC.num=length(cltrain$selected.index)

  pi_est_train<-cltrain$optimal_result$pi
  mu_est_train<-cltrain$optimal_result$mu
  sigma_est_train<-cltrain$optimal_result$sigma
  grp_assign_train<-apply(cltrain$optimal_result$z,1,which.max)
  cluster.train=grp_assign_train

  data1<-data.frame(cluster=cluster.train,Y=Y,Y.ind=delta.train,X=Z)
  mod.cox<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+X, data = data1)
  mod.cox1<-summary(mod.cox)
  coef<-mod.cox1$coefficients[1,1]

  index1<-which(cluster.train==1)
  index2<-which(cluster.train==2)
  if(coef>0){
    cluster.train[index2]<-1
    cluster.train[index1]<-2
  }

  cluster<-cluster.train
  cluster<-as.numeric(as.character(cluster))

  data<-data.frame(Y=Y,cluster=cluster,Y.ind=delta.train)
  surv_object <- Surv(time = data$Y,
                      event = data$Y.ind)
  fit.os<-surv_fit(surv_object ~ data$cluster,
                   data = data)

  p7<-ggsurvplot(fit.os, data = data,legend.title="Clusters",size=2,censor=FALSE
                 ,palette=1:2,title = "",risk.table=F,
                 pval=F,conf.int = F,xlab = "Time (Years)",legend="none",ylab="P(Surv)")%++%theme_bw()
  p7<-ggpar(p7,
            font.main = c(40, "bold"),
            font.x = c(40, "bold"),
            font.y = c(40, "bold"),
            font.caption = c(40, "bold"),
            font.legend = c(40, "bold"),
            font.tickslab = c(40, "bold"))
  p7<-p7%++%theme(legend.position="none",panel.border = element_rect(colour = "black",size=2))
  p7<-as.ggplot(p7$plot)

  data1<-data.frame(group=cluster.train,Y.ind=delta.train,Y=Y,x=Z[,1])
  mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x, data = data1)
  mod.cox2<-coxph(Surv(Y, Y.ind) ~ x, data = data1)
  mod.anova<-anova(mod.cox1,mod.cox2)
  F.PMBC<-mod.anova$Chisq[2]
  p.PMBC<-mod.anova$`Pr(>|Chi|)`[2]
  pred<-rep(NA,length(cluster.train))
  for(i in 1:K){
    data.temp<-data.frame(Y=Y[which(cluster.train==i)],x=Z[which(cluster.train==i)],
                          delta=delta.train[which(cluster.train==i)])
    mod.lm<-survreg(Surv(Y,delta)~x,dist="loglogistic",data=data.temp)
    pred[which(cluster.train==i)]<-predict(mod.lm,newdata = data.frame(x=Z[which(cluster.train==i)]))
  }

  cindex.PMBC<-concordance.index(x=1/as.numeric(pred), surv.time=Y, surv.event=delta.train, method="noether")$c.index
  cindex.adjust<-c()
  for(i in 1:500){
    cindex.adjust[i]<-concordance.index(x=1/as.numeric(sample(pred)), surv.time=Y, surv.event=delta.train, method="noether")$c.index
  }
  cindex.PMBC<-(cindex.PMBC-mean(cindex.adjust))/(1-mean(cindex.adjust))

  R2.gene.PMBC<-c()
  for(i in 1:ncol(X)){
    data1<-data.frame(gene=X[,i],cluster=as.factor(cluster.train))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.PMBC[i]<-1-RSS/TSS
  }
  R2.gene.PMBC<-R2.gene.PMBC[select.feature]

  mult_pdf=t(apply(X.test,1,function(x) {
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
  index1<-which(cluster.test==1)
  index2<-which(cluster.test==2)
  if(coef>0){
    cluster.test[index2]<-1
    cluster.test[index1]<-2
  }

  data<-data.frame(Y=Y.test,cluster=cluster.test,Y.ind=delta.test)
  surv_object <- Surv(time = data$Y,
                      event = data$Y.ind)
  fit.os<-surv_fit(surv_object ~ data$cluster,
                   data = data)

  p8<-ggsurvplot(fit.os, data = data,legend.title="Clusters",size=2,censor=FALSE
                 ,palette=1:2,title = "",risk.table=F,
                 pval=F,conf.int = F,xlab = "Time (Years)",legend="none",ylab="P(Surv)")%++%theme_bw()
  p8<-ggpar(p8,
            font.main = c(40, "bold"),
            font.x = c(40, "bold"),
            font.y = c(40, "bold"),
            font.caption = c(40, "bold"),
            font.legend = c(40, "bold"),
            font.tickslab = c(40, "bold"))
  p8<-p8%++%theme(legend.position="none",panel.border = element_rect(colour = "black",size=2))
  p8<-as.ggplot(p8$plot)

  data1<-data.frame(group=cluster.test,Y.ind=delta.test,Y=Y.test,x=Z.test)
  mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(group)+x, data = data1)
  mod.cox2<-coxph(Surv(Y, Y.ind) ~ x, data = data1)
  mod.anova<-anova(mod.cox1,mod.cox2)
  F.test.PMBC<-mod.anova$Chisq[2]
  p.test.PMBC<-mod.anova$`Pr(>|Chi|)`[2]
  pred<-rep(NA,length(cluster.test))
  for(i in 1:K){
    data.temp<-data.frame(Y=Y[which(cluster.train==i)],x=Z[which(cluster.train==i)],
                          delta=delta.train[which(cluster.train==i)])
    mod.lm<-survreg(Surv(Y,delta)~x,dist="loglogistic",data=data.temp)
    pred[which(cluster.test==i)]<-predict(mod.lm,newdata = data.frame(x=Z.test[which(cluster.test==i)]))
  }
  cindex.PMBC.test<-concordance.index(x=1/as.numeric(pred), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index
  cindex.adjust<-c()
  for(i in 1:500){
    cindex.adjust[i]<-concordance.index(x=1/as.numeric(sample(pred)), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index
  }
  cindex.PMBC.test<-(cindex.PMBC.test-mean(cindex.adjust))/(1-mean(cindex.adjust))

  R2.gene.test.PMBC<-c()
  for(i in 1:ncol(X.test)){
    data1<-data.frame(gene=X.test[,i],cluster=as.factor(cluster.test))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene.test.PMBC[i]<-1-RSS/TSS
  }
  R2.gene.test.PMBC<-R2.gene.test.PMBC[select.feature]

  p7 <- annotate_figure(p7,top = text_grob("(I.E1) outcome separation (discovery)",face = "bold", size = 20 ))
  p8 <- annotate_figure(p8,top = text_grob("(I.E2) outcome separation (validation)",face = "bold", size = 20 ))
  print("finish PMBC plot")

  group<-c(rep("ogClust_WJL",length(R2.gene.wcs)),
           rep("ogClust_GM",length(R2.gene.GM)),
           rep("skmeans",length(R2.gene.SKM)),
           rep("PMBC",length(R2.gene.PMBC)))
  group<-factor(group,levels=c("ogClust_WJL","ogClust_GM","skmeans","PMBC"))
  df<-data.frame(stat=c(R2.gene.wcs,R2.gene.GM, R2.gene.SKM,R2.gene.PMBC),group=group)
  p9 <- ggplot(df, aes(x=group,y=stat)) + xlab("")+ylab("R2(genes)")+theme_bw()+coord_cartesian(ylim = c(0,1))+
    geom_boxplot(fill="gray",lwd=1)+theme(legend.position = "right",
                                          plot.title = element_text(size = rel(2), hjust = 0.5),
                                          axis.title = element_text(size = 40,face="bold"),axis.text = element_text(size=25,face = "bold"),
                                          legend.text = element_text(size=25),panel.border = element_rect(colour = "black",size=2),
                                          axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))+
    scale_x_discrete(labels=c("ogClust_WJL" = expression(ogClust[WJL]), "ogClust_GM" = expression(ogClust[GM]),
                              "skmeans" = expression(paste("sparse ",italic(K),"-means")),
                              "PMBC"="PMBC"))

  group<-c(rep("ogClust_WJL",length(R2.gene.test.wcs)),
           rep("ogClust_GM",length(R2.gene.test.GM)),
           rep("skmeans",length(R2.gene.test.SKM)),
           rep("PMBC",length(R2.gene.test.PMBC)))
  group<-factor(group,levels=c("ogClust_WJL","ogClust_GM","skmeans","PMBC"))
  df<-data.frame(stat=c(R2.gene.test.wcs, R2.gene.test.GM, R2.gene.test.SKM,R2.gene.test.PMBC),group=group)
  p10 <- ggplot(df, aes(x=group,y=stat)) + xlab("")+ylab("R2(genes)")+theme_bw()+coord_cartesian(ylim = c(0,1))+
    geom_boxplot(fill="gray",lwd=1)+theme(legend.position = "right",
                                          plot.title = element_text(size = rel(2), hjust = 0.5),
                                          axis.title = element_text(size = 40,face="bold"),axis.text = element_text(size=25,face = "bold"),
                                          legend.text = element_text(size=25),panel.border = element_rect(colour = "black",size=2),
                                          axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))+
    scale_x_discrete(labels=c("ogClust_WJL" = expression(ogClust[WJL]), "ogClust_GM" = expression(ogClust[GM]),
                              "skmeans" = expression(paste("sparse ",italic(K),"-means")),
                              "PMBC"="PMBC"))

  p9 <- annotate_figure(p9,top = text_grob("(I) outcome separation (discovery)",face = "bold", size = 20 ))
  p10 <- annotate_figure(p10,top = text_grob("(J) outcome separation (validation)",face = "bold", size = 20 ))
  print("finish R2_genes plot")


  res.data=data.frame(F.stat=c(F.WJL,F.GM,F.SKM,F.PMBC),
                      pval=c(p.WJL,p.GM,p.SKM,p.PMBC),
                      cindex=c(cindex.WJL,cindex.GM,cindex.SKM,cindex.PMBC),
                      R2.gene=c(mean(R2.gene.wcs),mean(R2.gene.GM),mean(R2.gene.SKM),mean(R2.gene.PMBC)),
                      F.stat.test=c(F.test.WJL,F.test.GM,F.test.SKM,F.test.PMBC),
                      pval.test=c(p.test.WJL,p.test.GM,p.test.SKM,p.test.PMBC),
                      cindex.test=c(cindex.WJL.test,cindex.GM.test,cindex.SKM.test,cindex.PMBC.test),
                      R2.gene.test=c(mean(R2.gene.test.wcs,na.rm=T),mean(R2.gene.test.GM,na.rm=T),mean(R2.gene.test.SKM,na.rm=T),mean(R2.gene.test.PMBC,na.rm=T)),
                      ngenes=c(num.yujia.num,num.Peng.num,num.SKM.num,num.PMBC.num))
  rownames(res.data)=c("ogClustWJL","ogClust","SKM","PMBC")
  return(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,res.data))
}

plotlist=readIN.EXT(ngenes=400)
plotlist[[11]]
plot.WJL <- ggarrange(plotlist[[1]], plotlist[[2]], nrow=2)
plot.GM <- ggarrange(plotlist[[3]],plotlist[[4]],nrow=2)
plot.SKM<-ggarrange(plotlist[[5]],plotlist[[6]],nrow=2)
plot.PMBC<-ggarrange(plotlist[[7]], plotlist[[8]],nrow=2)

plot_combined<-grid.arrange(arrangeGrob(plot.GM,
                                        plot.WJL,
                                        plot.SKM,
                                        plot.PMBC,
                                        ncol=4),heights=c(10, 1))
file.name<-paste("FigureIV.png")
ggsave(plot_combined,filename = file.name, width = 70, height = 40, units = "cm")


plotlist=readIN.survEXT(ngenes=300)
plot.WJL <- ggarrange(plotlist[[1]], plotlist[[2]], nrow=2)
plot.GM <- ggarrange(plotlist[[3]],plotlist[[4]],nrow=2)
plot.SKM<-ggarrange(plotlist[[5]],plotlist[[6]],nrow=2)
plot.PMBC<-ggarrange(plotlist[[7]], plotlist[[8]],nrow=2)
plotlist[[11]]
plot_combined<-grid.arrange(arrangeGrob(plot.GM,
                                        plot.WJL,
                                        plot.SKM,
                                        plot.PMBC,
                                        ncol=4),heights=c(10, 1))

file.name<-paste("SupplementFigureVII_II.png")
ggsave(plot_combined,filename = file.name, width = 70, height = 40, units = "cm")

readIN.survEXTtune=function() {
  load("~/OneDrive - University of Pittsburgh/AOASyujia_peng/ogclust_Wei_Wenjia/TNBC/Data.Metabric.after.Combat_TNBC.Rdata")
  X.train<-Data.Metabric.after.Combat$Expression #dim 275 18964
  Y.train<-Data.Metabric.after.Combat$OS
  delta.train<-Data.Metabric.after.Combat$OS.event
  Z.train<-Data.Metabric.after.Combat$covariate
  Z.train<-as.matrix(Z.train[,1])

  Index.Mean<-order(apply(X.train,2,mean),decreasing = T)
  X.train<-X.train[,Index.Mean[1:(ncol(X.train)/2)]] #dim 275 9482
  Index.Sd<-order(apply(X.train,2,sd),decreasing = T)
  X.train<-X.train[,Index.Sd[1:(ncol(X.train)/2)]]
  X.train<-scale(X.train) #dim 275 4741

  load("survmod1.WJL.tune.RData")
  load("survmod.all.tune.RData")

  load("~/OneDrive - University of Pittsburgh/AOASyujia_peng/ogclust_Wei_Wenjia/TNBC/Data.ScanB.after.Combat_TNBC.Rdata")
  X.test<-Data.ScanB.after.Combat$Expression
  Y.test<-Data.ScanB.after.Combat$OS
  delta.test<-Data.ScanB.after.Combat$OS.event
  Z.test<-Data.ScanB.after.Combat$covariate
  Z.test<-as.matrix(Z.test$Age)
  X.test<-scale(X.test)

  index<-match(colnames(X.train),colnames(X.test))
  X.test<-X.test[,index]
  index<-which(apply(X.test,2,function(x){sum(is.na(x))})!=0)
  X.test[,index]<-0

  plotlist=survEXTplot.gen(X=X.train,Y=Y.train,Z=Z.train,delta.train=delta.train,mod.WJL=mod.WJL,modfit=modfit,ngenes=ngenes,K=2,
                           X.test=X.test,Z.test=Z.test,Y.test=Y.test,delta.test=delta.test)
  return(plotlist)
}

plotlist=readIN.survEXTtune()
plot.WJL <- ggarrange(plotlist[[1]], plotlist[[2]], nrow=2)
plot.GM <- ggarrange(plotlist[[3]],plotlist[[4]],nrow=2)
plot.SKM<-ggarrange(plotlist[[5]],plotlist[[6]],nrow=2)
plot.PMBC<-ggarrange(plotlist[[7]], plotlist[[8]],nrow=2)
plotlist[[11]]
plot_combined<-grid.arrange(arrangeGrob(plot.GM,
                                        plot.WJL,
                                        plot.SKM,
                                        plot.PMBC,
                                        ncol=4),heights=c(10, 1))

file.name<-paste("SupplementFigureVII_I.png")
ggsave(plot_combined,filename = file.name, width = 70, height = 40, units = "cm")




