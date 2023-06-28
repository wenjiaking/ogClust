rm(list=ls())
library(sparcl) #to implement sparse K means
library(adaHuber)
library(parallel)
library(Brobdingnag)
library(glmnet)
library(mclust)
library(psych)
library(hms)
library(survcomp)
#library(survminer)
library(survival)
library(tibble)
library(dplyr)

source("Functions.R")

n<-99
q<-50
q1<-50
q2<-50
q3<-1850
K<-3
num.sim<-50
mu_vector<-c(0.9,1.2,1.5,1.8)
mu1_vector<-mu_vector+0.2
beta_vector<-0.5
sigma_y<-0.5
c1=2
var_g=2

for (ind_mu in 1:length(mu_vector)) {
  mu<-mu_vector[ind_mu]
  mu1<-mu1_vector[ind_mu]
  file.name<-paste("SimII_Data_c1=",c1,"_mu=",mu,".rds",sep="")
  data.temp=mclapply(1:num.sim,function(ind.data) {
    set.seed(ind.data)
    data<-Sim_surv(n=n,beta1=beta_vector,beta2=beta_vector,q=q,q1=q1,q2=q2,q3=q3,c1=c1,var_g=var_g,mu=mu,mu1=mu1,sigma_y=sigma_y,censor=200)
    X<-data$x
    Y<-data$y
    G<-data$G
    delta<-data$y.ind
    true.label<-c(rep(1,33),rep(2,33),rep(3,33))
    true.feature<-c(rep(1,q),rep(0,2000-q))
    return(list(X=X,Y=Y,G=G,delta=delta,true.label=true.label,true.feature=true.feature))
  },mc.cores=num.sim)
  saveRDS(data.temp,file.name)
}


resSurv.gen=function(K=3,SC=TRUE,num.sim=50,mu_vector=c(0.9,1.2,1.5,1.8),c1_vector_all=c(2)) {

  for(ind_c1 in 1:length(c1_vector_all)) {
    c1<-c1_vector_all[ind_c1]
    for (ind_mu in 1:length(mu_vector)) {
      mu=mu_vector[ind_mu]
      file.name<-paste("SimII_Data_c1=",c1,"_mu=",mu,".rds",sep="")
      print(paste0("Begin ",file.name))
      data.temp=readRDS(file.name)
      RES=mclapply(1:num.sim,function(ind.data) {
        data=data.temp[[ind.data]]
        X<-data$X
        Y<-data$Y
        G<-data$G
        delta=data$delta
        true.label<-data$true.label
        true.feature<-data$true.feature
        set.seed(ind.data)
        index.test<-sample(1:nrow(X),round(nrow(X)/3))
        index.train<-(1:nrow(X))[-index.test]
        X.train<-X[index.train,]
        G.train<-G[index.train,]
        G.train<-scale(G.train)
        Y.train<-Y[index.train]
        delta.train<-delta[index.train]
        X.test<-X[index.test,]
        G.test<-G[index.test,]
        G.test<-scale(G.test)
        Y.test<-Y[index.test]
        delta.test<-delta[index.test]
        label.train<-true.label[index.train]
        label.test<-true.label[index.test] #1/3 samples for testing and 2/3 for training

        if (SC) {
          t.stat<-rep(NA,ncol(G.train))
          for(i in 1:ncol(G.train)){
            data<-data.frame(y=Y.train,x=G.train[,i])
            mod<-lm(y~x,data=data)
            mod<-summary(mod)
            t.stat[i]<-mod$coefficients[2,3]
          }
          index<-order(abs(t.stat),decreasing = T)[1:500]
          true.feature<-true.feature[index]
          G.train<-G.train[,index]
          G.test<-G.test[,index]
        }

        set.seed(12315)
        index.sample<-sample(rep(1:10,length.out=nrow(X.train)))
        #--------------------------------------------------
        #WJL's method
        #--------------------------------------------------
        mod.kmeans<-kmeans(G.train,center=K,nstart=50)
        center<-t(mod.kmeans$centers)

        w_vector=seq(0,1,0.1)
        lambda_vector=seq(2,15,1)
        s_G<-50
        lambda_vector1<-seq(2,18,0.5)

        start_time=Sys.time()
        temp.R2<-rep(NA,length(w_vector))
        for(ind_w in 1:length(w_vector)){
          # print(ind_w)
          w<-w_vector[ind_w]
          w<-(s_G*w)/(s_G*w+1-w)
          num.vector<-rep(NA,length(lambda_vector))
          for(j in 1:length(lambda_vector)){
            mod<-ogClust_Surv(x=X.train,G=t(G.train),y=Y.train,y.ind=delta.train,c_center=center,lambda=lambda_vector[j],v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w)
            num.vector[j]<-sum(apply(mod$result_list$mu,1,function(x){length(unique(x))})!=1)
          }
          index<-which.min(abs(num.vector-50))
          lambda<-lambda_vector[index]

          mod.all<-ogClust_Surv(x=X.train,G=t(G.train),y=Y.train,y.ind=delta.train,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w)
          #cluster1<-predict.ogclust.test(mod.all,K=K,D.test=as.matrix(G.train),D.train=t(G.train),X1=X.train,p=ncol(G.train))
          cluster1<-predict_test(mod.all,K=K,D.test=as.matrix(G.train),X1=as.matrix(X.train),p=ncol(G.train))
          cluster1<-cluster1$clus
          if (length(table(cluster1))!=K) {
            temp.R2[ind_w]=NA
            next
          }
          centers<-cbind(apply(G.train[which(cluster1==1),,drop=F],2,mean),
                         apply(G.train[which(cluster1==2),,drop=F],2,mean),
                         apply(G.train[which(cluster1==3),,drop=F],2,mean))

          #print(table(cluster1))

          cluster<-rep(NA,nrow(G.train))
          for(ind_cv in 1:length(unique(index.sample))){
            ind_sample<-which(index.sample==ind_cv)
            G.test1<-G.train[ind_sample,]
            X.test1<-X.train[ind_sample,]
            Y.test1<-Y.train[ind_sample]
            delta.test1<-delta.train[ind_sample]
            G.train1<-G.train[(1:nrow(X.train))[-ind_sample],]
            X.train1<-X.train[(1:nrow(X.train))[-ind_sample],]
            Y.train1<-Y.train[(1:nrow(X.train))[-ind_sample]]
            delta.train1<-delta.train[(1:nrow(X.train))[-ind_sample]]

            label.test1<-label.train[ind_sample]
            label.train1<-label.train[(1:nrow(X.train))[-ind_sample]]
            num.vector<-rep(NA,length(lambda_vector))
            for(j in 1:length(lambda_vector)){
              mod<-ogClust_Surv(x=X.train1,G=t(G.train1),y=Y.train1,y.ind=delta.train1,c_center=center,lambda=lambda_vector[j],v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w_vector[ind_w],w_G=1-w_vector[ind_w])
              num.vector[j]<-sum(apply(mod$result_list$mu,1,function(x){length(unique(x))})!=1)
            }

            index<-which.min(abs(num.vector-50))
            lambda<-lambda_vector[index]
            #temp.num[ind_w,ind_cv]<-num.vector[index]
            mod<-ogClust_Surv(x=X.train1,G=t(G.train1),y=Y.train1,y.ind=delta.train1,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w)
            res1<-predict_test(mod,K=K,D.test=as.matrix(G.test1),X1=as.matrix(X.test1),p=ncol(G.test1))
            #res1<-predict.ogclust.test(mod,K=K,D.test=as.matrix(G.test1),D.train=t(G.train1),X1=X.test1,p=ncol(G.test1))

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

          }
          data1<-data.frame(cluster=cluster,Y=Y.train,Y.ind=delta.train,x1=X.train[,1],x2=X.train[,2])
          mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
          mod.cox2<-coxph(Surv(Y, Y.ind) ~ x2+x2, data = data1)
          mod.anova<-anova(mod.cox1,mod.cox2)
          temp.R2[ind_w]<-mod.anova$Chisq[2]
        }
        index<-which.max(temp.R2)
        end_time=Sys.time()
        tune_time=as_hms(end_time-start_time)

        w0<-w_vector[index]
        w<-(s_G*w0)/(s_G*w0+1-w0)
        print(paste0("Tuning w0=",w0," consumed time  ",tune_time))

        R2.outcome.tune<-rep(NA,length(lambda_vector1))
        num.vector<-rep(NA,length(lambda_vector1))
        R2.gene.tune<-rep(NA,length(lambda_vector1))
        for(j in 1:length(lambda_vector1)){
          mod<-ogClust_Surv(x=X.train,G=t(G.train),y=Y.train,y.ind=delta.train,c_center=center,lambda=lambda_vector1[j],v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w)
          select.feature<-as.numeric(apply(mod$result_list$mu,1,function(x){length(unique(x))})!=1)
          num.vector[j]<-sum(select.feature==1)

          cluster1<-predict_test(mod,K=K,D.test=as.matrix(G.train),X1=as.matrix(X.train),p=ncol(G.train))
          cluster1<-cluster1$clus

          if(length(table(cluster1))!=K){
            R2.outcome.tune[j]=NA
            R2.gene.tune[j]=NA
            next
          }

          centers<-cbind(apply(G.train[which(cluster1==1),,drop=F],2,mean),
                         apply(G.train[which(cluster1==2),,drop=F],2,mean),
                         apply(G.train[which(cluster1==3),,drop=F],2,mean))


          cluster<-rep(NA,nrow(G.train))
          for(ind_cv in 1:length(unique(index.sample))){
            ind_sample<-which(index.sample==ind_cv)
            G.test11<-G.train[ind_sample,]
            X.test1<-X.train[ind_sample,]
            Y.test1<-Y.train[ind_sample]
            delta.test1<-delta.train[ind_sample]
            delta.train1<-delta.train[(1:nrow(X.train))[-ind_sample]]
            G.train11<-G.train[(1:nrow(G.train))[-ind_sample],]
            X.train1<-X.train[(1:nrow(G.train))[-ind_sample],]
            Y.train1<-Y.train[(1:nrow(G.train))[-ind_sample]]

            mod<-ogClust_Surv(x=X.train1,G=t(G.train11),y=Y.train1,y.ind=delta.train1,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=3,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)

            res1<-predict_test(mod,K=K,D.test=as.matrix(G.test11),X1=as.matrix(X.test1),p=ncol(G.test11))

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

          data1<-data.frame(cluster=cluster,Y=Y.train,Y.ind=delta.train,x1=X.train[,1],x2=X.train[,2])
          mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
          mod.cox2<-coxph(Surv(Y, Y.ind) ~ x2+x2, data = data1)
          mod.anova<-anova(mod.cox1,mod.cox2)
          R2.outcome.tune[j]<-mod.anova$Chisq[2]

          Fstat.gene<-c()
          for(i in 1:ncol(G.train)){
            data1<-data.frame(gene=G.train[,i],cluster=as.factor(cluster))
            mod<-lm(gene~cluster,data=data1)
            mod<-summary(mod)
            Fstat.gene[i]<-abs(mod$coefficients[2,3])
          }
          Fstat.gene<-mean(Fstat.gene[which(select.feature==1)])
          R2.gene.tune[j]<-Fstat.gene
        }

        #---------------------------------------------
        #use criteria sqrt(R2.outcome.tune*R2.gene.tune) to tune lambda
        #---------------------------------------------
        index<-which.max(sqrt(R2.outcome.tune*R2.gene.tune))
        lambda=lambda_vector1[index]
        print(paste0("By sqrt(R2.outcome*R2.gene) criteria, lambda=",lambda," that selects  ",num.vector[index]," features."))

        mod.tune<-ogClust_Surv(x=X.train,G=t(G.train),y=Y.train,y.ind=delta.train,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w)
        select.feature<-as.numeric(apply(mod.tune$result_list$mu,1,function(x){length(unique(x))})>1)
        num.WJL.tune=sum(select.feature)

        res.test.tune<-predict_test(mod.tune,K=K,D.test=G.test,X1=as.matrix(X.test),p=ncol(G.test))
        res.train.tune<-predict.ogclust.test.select(mod.tune,K=K,D.test=as.matrix(G.train),D.train=t(G.train),X1=as.matrix(X.train),p=ncol(G.train),s_G = s_G,O.test = Y.train,w = w)

        cluster.test<-res.test.tune$clus
        cluster.train = res.train.tune$clus

        ari.test.WJL.tune<-adjustedRandIndex(cluster.test,label.test)
        ari.train.WJL.tune<-adjustedRandIndex(cluster.train,label.train)

        data1<-data.frame(cluster=cluster.train,Y=Y.train,Y.ind=delta.train,x1=X.train[,1],x2=X.train[,2])
        #library(survival)
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.train.tune<-mod.anova$Chisq[2]

        data1<-data.frame(cluster=cluster.test,Y=Y.test,Y.ind=delta.test,x1=X.test[,1],x2=X.test[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.test.tune<-mod.anova$Chisq[2]

        mod.train<-Rsquare(cluster =cluster.train,Y = Y.train,X = X.train,G = t(G.train))
        R2.gene.train.tune<-mean(mod.train$gene[which(select.feature==1)])
        mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test))
        R2.gene.test.tune<-mean(mod.test$gene[which(select.feature==1)])

        cindex.train.WJL.tune.hard<-concordance.index(x=1/as.numeric(res.train.tune$Y.hard), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
        cindex.test.WJL.tune.hard<-concordance.index(x=1/as.numeric(res.test.tune$Y.hard), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        cindex.train.random<-c()
        cindex.test.random<-c()
        cluster.test1<-cluster.test
        cluster.train1<-cluster.train
        for(ind in 1:50){
          cluster.test1<-sample(cluster.test1)
          cluster.train1<-sample(cluster.train1)
          intercept<-rep(NA,nrow(X.train))
          for(i in 1:K){
            intercept[which(cluster.train1==i)]<-mod.tune$result_list$int_coef[i]
          }

          pred.train<-X.train %*% mod.tune$result_list$int_coef[4:5]+intercept

          intercept<-rep(NA,nrow(X.test))
          for(i in 1:K){
            intercept[which(cluster.test1==i)]<-mod.tune$result_list$int_coef[i]
          }

          pred.test<-X.test %*% mod.tune$result_list$int_coef[4:5]+intercept


          cindex.train.random[ind]<-concordance.index(x=1/as.numeric(pred.train), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
          cindex.test.random[ind]<-concordance.index(x=1/as.numeric(pred.test), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        }

        cindex.train.WJL.tune.adjust<-(cindex.train.WJL.tune.hard-mean(cindex.train.random))/(1-mean(cindex.train.random))
        cindex.test.WJL.tune.adjust<-(cindex.test.WJL.tune.hard-mean(cindex.test.random))/(1-mean(cindex.test.random))

        jaccard.index.WJL.tune<-Jaccard.index(true.feature, select.feature)

        #---------------------------------------------
        #use num to select lambda
        #---------------------------------------------
        index<-max(which(num.vector>=50))
        lambda.num<-lambda_vector1[index]
        print(paste0("Tuning lambda=",lambda.num," that selects  ",num.vector[index]," features."))

        start_time=Sys.time()
        mod.num<-ogClust_Surv(x=X.train,G=t(G.train),y=Y.train,y.ind=delta.train,c_center=center,lambda=lambda.num,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
        end_time=Sys.time()
        train_time=as_hms(end_time-start_time)
        time.WJL=as_hms(train_time+tune_time)
        #res.test<-predict.ogclust.test(mod.num,K=K,D.test=G.test,D.train=t(G.train),X1=X.test,p=ncol(G.test))

        res.test<-predict_test(mod.num,K=K,D.test=G.test,X1=as.matrix(X.test),p=ncol(G.test))
        res.train<-predict.ogclust.test.select(mod.num,K=K,D.test=as.matrix(G.train),D.train=t(G.train),X1=X.train,p=ncol(G.train),s_G = s_G,O.test = Y.train,w = w)

        cluster.test<-res.test$clus
        #cluster.train<-apply(mod.num$result_list$z,1,which.max)
        cluster.train = res.train$clus

        ari.test.WJL.num<-adjustedRandIndex(cluster.test,label.test)
        ari.train.WJL.num<-adjustedRandIndex(cluster.train,label.train)


        select.feature<-as.numeric(apply(mod.num$result_list$mu,1,function(x){length(unique(x))})>1)
        num.WJL.num=sum(select.feature)
        jaccard.index.WJL.num<-Jaccard.index(true.feature, select.feature)

        data1<-data.frame(cluster=cluster.train,Y=Y.train,Y.ind=delta.train,x1=X.train[,1],x2=X.train[,2])
        #library(survival)
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.train.num<-mod.anova$Chisq[2]

        data1<-data.frame(cluster=cluster.test,Y=Y.test,Y.ind=delta.test,x1=X.test[,1],x2=X.test[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.test.num<-mod.anova$Chisq[2]

        mod.train<-Rsquare(cluster =cluster.train,Y = Y.train,X = X.train,G = t(G.train))
        R2.gene.train.num<-mean(mod.train$gene[which(select.feature==1)])
        mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test))
        R2.gene.test.num<-mean(mod.test$gene[which(select.feature==1)])


        cindex.train.WJL.num.hard<-concordance.index(x=1/as.numeric(res.train$Y.hard), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
        cindex.test.WJL.num.hard<-concordance.index(x=1/as.numeric(res.test$Y.hard), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        cindex.train.random<-c()
        cindex.test.random<-c()
        cluster.test1<-cluster.test
        cluster.train1<-cluster.train
        for(ind in 1:50){
          cluster.test1<-sample(cluster.test1)
          cluster.train1<-sample(cluster.train1)
          intercept<-rep(NA,nrow(X.train))
          for(i in 1:K){
            intercept[which(cluster.train1==i)]<-mod.num$result_list$int_coef[i]
          }

          pred.train<-X.train %*% mod.num$result_list$int_coef[4:5]+intercept

          intercept<-rep(NA,nrow(X.test))
          for(i in 1:K){
            intercept[which(cluster.test1==i)]<-mod.num$result_list$int_coef[i]
          }

          pred.test<-X.test %*% mod.num$result_list$int_coef[4:5]+intercept


          cindex.train.random[ind]<-concordance.index(x=1/as.numeric(pred.train), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
          cindex.test.random[ind]<-concordance.index(x=1/as.numeric(pred.test), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        }

        cindex.train.WJL.num.adjust<-(cindex.train.WJL.num.hard-mean(cindex.train.random))/(1-mean(cindex.train.random))
        cindex.test.WJL.num.adjust<-(cindex.test.WJL.num.hard-mean(cindex.test.random))/(1-mean(cindex.test.random))

        #-------------------------------------------------
        WJL.res<-list(mod.tune=mod.tune,mod.num=mod.num,
                      pred.train.tune=res.train.tune,pred.train.num=res.train,
                      pred.test.tune=res.test.tune,pred.test.num=res.test,
                      ari.train.tune=ari.train.WJL.tune,ari.train.num=ari.train.WJL.num,
                      ari.test.tune=ari.test.WJL.tune,ari.test.num=ari.test.WJL.num,
                      jaccard.index.tune=jaccard.index.WJL.tune,jaccard.index.num=jaccard.index.WJL.num,
                      F.train.tune=F.train.tune,F.test.tune=F.test.tune,
                      F.train.num=F.train.num,F.test.num=F.test.num,
                      num.tune=num.WJL.tune,num.num=num.WJL.num,
                      cindex.train.num=cindex.train.WJL.num.hard,
                      cindex.test.num=cindex.test.WJL.num.hard,
                      cindex.train.num.adjust=cindex.train.WJL.num.adjust,
                      cindex.test.num.adjust=cindex.test.WJL.num.adjust,
                      cindex.train.tune=cindex.train.WJL.tune.hard,
                      cindex.test.tune=cindex.test.WJL.tune.hard,
                      cindex.train.tune.adjust=cindex.train.WJL.tune.adjust,
                      cindex.test.tune.adjust=cindex.test.WJL.tune.adjust,
                      R2.gene.train.tune=R2.gene.train.tune,
                      R2.gene.test.tune=R2.gene.test.tune,
                      R2.gene.train.num=R2.gene.train.num,
                      R2.gene.test.num=R2.gene.test.num,
                      lambda.tune=lambda,lambda.num=lambda.num,
                      train_time=time.WJL,
                      w0.est=w0)
        print("Finish WJL method!")
        #------------------------------
        #Peng's method
        #------------------------------
        n1=nrow(G.train) # number of samples
        NG=ncol(G.train) # number of genes
        np=ncol(X.train) # number of covariates
        lambda_vector_Peng<-c(seq(0.005,0.01,0.001),seq(0.02,0.5,0.02))

        mod.kmeans<-kmeans(G.train,centers = 3,nstart = 50)
        cluster<-mod.kmeans$cluster

        index1<-which(cluster==1)
        data1<-data.frame(y=Y.train[index1],x1=X.train[index1,1],x2=X.train[index1,2])
        mod1<-lm(y~.,data=data1)
        mod1<-summary(mod1)

        index2<-which(cluster==2)
        data2<-data.frame(y=Y.train[index2],x1=X.train[index2,1],x2=X.train[index2,2])
        mod2<-lm(y~.,data=data2)
        mod2<-summary(mod2)

        index3<-which(cluster==3)
        data3<-data.frame(y=Y.train[index3],x1=X.train[index3,1],x2=X.train[index3,2])
        mod3<-lm(y~.,data=data3)
        mod3<-summary(mod3)

        beta0_int = c(mod1$coefficients[1,1],mod2$coefficients[1,1],mod3$coefficients[1,1])

        data<-data.frame(y=Y.train,X.train)
        mod<-lm(y~.,data=data)
        mod<-summary(mod)
        beta_int<-mod$coefficients[2:3,1]
        K<-3

        mod<-glmnet(G.train,as.factor(cluster),family = "multinomial",lambda=0)
        mod1<-coef(mod)
        gamma_int<-c(as.numeric(mod1$`1`),as.numeric(mod1$`2`))
        sigma2_int<-1
        theta_int = c(beta_int, gamma_int, beta0_int, sigma2_int)

        R2.outcome.tune<-rep(NA,length(lambda_vector_Peng))
        R2.gene.tune<-rep(NA,length(lambda_vector_Peng))
        num<-rep(NA,length(lambda_vector_Peng))
        for(ind_lambda in 1:length(lambda_vector_Peng)){
          #print(i)
          fit.res<-fit.ogClust.surv(n=n1, K=K, np=np, NG=NG, lambda=lambda_vector_Peng[ind_lambda],delta = delta.train,
                                    alpha=0.5,G=as.data.frame(G.train), Y=Y.train, X=X.train, theta_int=theta_int)

          # theta_est<-fit.res$res
          beta=fit.res$par$gamma
          beta<-as.data.frame(beta)
          select.feature<-as.numeric(apply(beta,1,function(x){sum(x==0)})!=K)[-1] #genes with non-zero weight of at least one clsuter
          num[ind_lambda]<-sum(select.feature==1)
          cluster1<-fit.res$grp_assign
          #print(table(cluster1))

          if(length(table(cluster1))!=K){
            R2.outcome.tune[ind_lambda]<-NA
            R2.gene.tune[ind_lambda]<-NA
            next
          }

          centers<-cbind(apply(G.train[which(cluster1==1),,drop=F],2,mean),
                         apply(G.train[which(cluster1==2),,drop=F],2,mean),
                         apply(G.train[which(cluster1==3),,drop=F],2,mean))
          #print(paste0("Begin the 10-fold CV with temporary lambda=",lambda))
          cluster<-rep(NA,nrow(G.train))
          for(ind_cv in 1:length(unique(index.sample))){
            ind_sample<-which(index.sample==ind_cv)
            G.test1<-G.train[ind_sample,]
            X.test1<-X.train[ind_sample,]
            Y.test1<-Y.train[ind_sample]
            delta.test1<-delta.train[ind_sample]
            delta.train1<-delta.train[(1:nrow(X.train))[-ind_sample]]
            G.train1<-G.train[(1:nrow(G.train))[-ind_sample],]
            X.train1<-X.train[(1:nrow(G.train))[-ind_sample],]
            Y.train1<-Y.train[(1:nrow(G.train))[-ind_sample]]
            n1sub=nrow(G.train1)
            NGsub=ncol(G.train1)
            np=ncol(X.train1)
            mod<-fit.ogClust.surv(n=n1sub, K=K, np=np, NG=NGsub, lambda=lambda_vector_Peng[ind_lambda],delta = delta.train1,
                                  alpha=0.5, G=as.data.frame(G.train1), Y=Y.train1, X=X.train1, theta_int=theta_int)


            cluster1<-mod$grp_assign
            gamma_est_matrix=mod$par$gamma

            G.test1.add = as.matrix(cbind(1, G.test1))
            pai_est.num = sapply(1:K, function(k) exp(G.test1.add %*% gamma_est_matrix[,  k, drop = F])/rowSums(exp(G.test1.add %*% gamma_est_matrix)))
            cluster.test<-apply(pai_est.num,1,which.max)

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

          if(length(table(cluster))!=K){
            R2.outcome.tune[ind_lambda]<-NA
            R2.gene.tune[ind_lambda]<-NA
          }else{
            data1<-data.frame(cluster=cluster,Y=Y.train,Y.ind=delta.train,x1=X.train[,1],x2=X.train[,2])
            mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
            mod.cox2<-coxph(Surv(Y, Y.ind) ~ x2+x2, data = data1)
            mod.anova<-anova(mod.cox1,mod.cox2)
            R2.outcome.tune[ind_lambda]<-mod.anova$Chisq[2]

            Fstat.gene<-c()
            for(i in 1:ncol(G.train)){
              data1<-data.frame(gene=G.train[,i],cluster=as.factor(cluster))
              mod<-lm(gene~cluster,data=data1)
              mod<-summary(mod)
              Fstat.gene[i]<-abs(mod$coefficients[2,3])
            }
            Fstat.gene<-mean(Fstat.gene[which(select.feature==1)])
            R2.gene.tune[ind_lambda]<-Fstat.gene
          }
        }

        #---------------------------------------------
        #use criteria sqrt(R2.outcome.tune*R2.gene.tune) to tune lambda
        #---------------------------------------------
        tune_metric=sqrt(R2.outcome.tune*R2.gene.tune)
        print(tune_metric)
        index<-which.max(tune_metric)
        lambda=lambda_vector_Peng[index]
        print(paste0("By sqrt(R2.outcome*R2.gene) criteria, lambda=",lambda," that selects  ",num[index]," features."))

        fit.res.tune<-fit.ogClust.surv(n=n1, K=K, np=np, NG=NG, lambda=lambda,delta = delta.train,
                                       alpha=0.5, G=as.data.frame(G.train), Y=Y.train, X=X.train, theta_int=theta_int)


        beta=fit.res.tune$par$gamma
        beta<-as.data.frame(beta)
        select.feature<-as.numeric(apply(beta,1,function(x){sum(x==0)})!=K)[-1] #genes with non-zero weight of at least one clsuter

        jaccard.index.Peng.tune<-Jaccard.index(true.feature, select.feature)
        num.Peng.tune<-sum(select.feature==1)
        ari.train.Peng.tune<-adjustedRandIndex(fit.res.tune$grp_assign,label.train)

        data1<-data.frame(cluster=fit.res.tune$grp_assign,Y=Y.train,Y.ind=delta.train,x1=X.train[,1],x2=X.train[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.train.tune<-mod.anova$Chisq[2]

        beta_est=fit.res.tune$par$beta
        gamma_est_matrix=fit.res.tune$par$gamma
        beta0_est=fit.res.tune$par$beta0
        sigma2_est=fit.res.tune$par$sigma2

        intercept<-rep(NA,nrow(X.train))
        for(i in 1:K){
          intercept[which(fit.res.tune$grp_assign==i)]<-beta0_est[i]
        }

        pred.train<-X.train %*% beta_est+intercept
        cindex.train.Peng.tune<-concordance.index(x=1/as.numeric(pred.train), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
        mod.train<-Rsquare(cluster =fit.res.tune$grp_assign,Y = Y.train,X = X.train,G = t(G.train))
        R2.gene.train.tune<-mean(mod.train$gene[which(select.feature==1)])

        G.test.add = cbind(1, G.test)
        pai_est.tune = sapply(1:K, function(k) exp(G.test.add %*% gamma_est_matrix[,  k, drop = F])/rowSums(exp(G.test.add %*% gamma_est_matrix)))
        cluster.test<-apply(pai_est.tune,1,which.max)
        ari.test.Peng.tune<-adjustedRandIndex(cluster.test,label.test)

        data1<-data.frame(cluster=cluster.test,Y=Y.test,Y.ind=delta.test,x1=X.test[,1],x2=X.test[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.test.tune<-mod.anova$Chisq[2]

        mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test))
        R2.gene.test.tune<-mean(mod.test$gene[which(select.feature==1)])

        intercept<-rep(NA,nrow(X.test))
        for(i in 1:K){
          intercept[which(cluster.test==i)]<-beta0_est[i]
        }

        pred.test<-X.test %*% beta_est+intercept
        cindex.test.Peng.tune.hard<-concordance.index(x=1/as.numeric(pred.test), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        cindex.train.random<-c()
        cindex.test.random<-c()
        cluster.test1<-cluster.test
        cluster.train1<-fit.res.tune$grp_assign
        for(ind in 1:50){
          cluster.test1<-sample(cluster.test1)
          cluster.train1<-sample(cluster.train1)
          intercept<-rep(NA,nrow(X.train))
          for(i in 1:K){
            intercept[which(cluster.train1==i)]<-beta0_est[i]
          }

          pred.train<-X.train %*% beta_est+intercept

          intercept<-rep(NA,nrow(X.test))
          for(i in 1:K){
            intercept[which(cluster.test1==i)]<-beta0_est[i]
          }

          pred.test<-X.test %*% beta_est+intercept


          cindex.train.random[ind]<-concordance.index(x=1/as.numeric(pred.train), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
          cindex.test.random[ind]<-concordance.index(x=1/as.numeric(pred.test), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        }

        cindex.train.Peng.tune.adjust<-(cindex.train.Peng.tune-mean(cindex.train.random))/(1-mean(cindex.train.random))
        cindex.test.Peng.tune.adjust<-(cindex.test.Peng.tune.hard-mean(cindex.test.random))/(1-mean(cindex.test.random))

        #---------------------------------------------
        #use num to select lambda
        #---------------------------------------------
        index<-which.min(abs(num-50))
        # index<-min(which(num>50))
        lambda.num=lambda_vector_Peng[index]
        print(paste0("Tuning lambda=",lambda.num," that selects  ",num[index]," features."))
        start_time=Sys.time()
        fit.res.num<-fit.ogClust.surv(n=n1, K=K, np=np, NG=NG, lambda=lambda.num,delta = delta.train,
                                      alpha=0.5, G=as.data.frame(G.train), Y=Y.train, X=X.train, theta_int=theta_int)
        end_time=Sys.time()
        time.GM=as_hms(end_time-start_time)

        beta=fit.res.num$par$gamma
        beta<-as.data.frame(beta)
        select.feature<-as.numeric(apply(beta,1,function(x){sum(x==0)})!=K)[-1]

        jaccard.index.Peng.num<-Jaccard.index(true.feature, select.feature)
        num.Peng.num<-sum(select.feature==1)
        ari.train.Peng.num<-adjustedRandIndex(fit.res.num$grp_assign,label.train)
        data1<-data.frame(cluster=fit.res.num$grp_assign,Y=Y.train,Y.ind=delta.train,x1=X.train[,1],x2=X.train[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.train.num<-mod.anova$Chisq[2]

        beta_est=fit.res.num$par$beta
        gamma_est_matrix=fit.res.num$par$gamma
        beta0_est=fit.res.num$par$beta0
        sigma2_est=fit.res.num$par$sigma2
        intercept<-rep(NA,nrow(X.train))
        for(i in 1:K){
          intercept[which(fit.res.num$grp_assign==i)]<-beta0_est[i]
        }

        pred.train<-X.train %*% beta_est+intercept
        cindex.train.Peng.num<-concordance.index(x=1/as.numeric(pred.train), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index

        mod.train<-Rsquare(cluster =fit.res.num$grp_assign,Y = Y.train,X = X.train,G = t(G.train))
        R2.gene.train.num<-mean(mod.train$gene[which(select.feature==1)])


        G.test.add = cbind(1, G.test)
        pai_est.num = sapply(1:K, function(k) exp(G.test.add %*% gamma_est_matrix[,  k, drop = F])/rowSums(exp(G.test.add %*% gamma_est_matrix)))

        cluster.test<-apply(pai_est.num,1,which.max)
        ari.test.Peng.num<-adjustedRandIndex(cluster.test,label.test)

        data1<-data.frame(cluster=cluster.test,Y=Y.test,Y.ind=delta.test,x1=X.test[,1],x2=X.test[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.test.num<-mod.anova$Chisq[2]

        mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test))
        R2.gene.test.num<-mean(mod.test$gene[which(select.feature==1)])

        intercept<-rep(NA,nrow(X.test))
        for(i in 1:K){
          intercept[which(cluster.test==i)]<-beta0_est[i]
        }
        pred.test<-X.test %*% beta_est+intercept
        cindex.test.Peng.num.hard<-concordance.index(x=1/as.numeric(pred.test), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        cindex.train.random<-c()
        cindex.test.random<-c()
        cluster.test1<-cluster.test
        cluster.train1<-fit.res.num$grp_assign
        for(ind in 1:50){
          cluster.test1<-sample(cluster.test1)
          cluster.train1<-sample(cluster.train1)
          intercept<-rep(NA,nrow(X.train))
          for(i in 1:K){
            intercept[which(cluster.train1==i)]<-beta0_est[i]
          }

          pred.train<-X.train %*% beta_est+intercept

          intercept<-rep(NA,nrow(X.test))
          for(i in 1:K){
            intercept[which(cluster.test1==i)]<-beta0_est[i]
          }

          pred.test<-X.test %*% beta_est+intercept


          cindex.train.random[ind]<-concordance.index(x=1/as.numeric(pred.train), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
          cindex.test.random[ind]<-concordance.index(x=1/as.numeric(pred.test), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        }

        cindex.train.Peng.num.adjust<-(cindex.train.Peng.num-mean(cindex.train.random))/(1-mean(cindex.train.random))
        cindex.test.Peng.num.adjust<-(cindex.test.Peng.num.hard-mean(cindex.test.random))/(1-mean(cindex.test.random))

        #-------------------------------------------------
        GM.res<-list(mod.tune=fit.res.tune,mod.num=fit.res.num,
                     ari.train.tune=ari.train.Peng.tune,ari.train.num=ari.train.Peng.num,
                     ari.test.tune=ari.test.Peng.tune,ari.test.num=ari.test.Peng.num,
                     jaccard.index.tune=jaccard.index.Peng.tune,jaccard.index.num=jaccard.index.Peng.num,
                     num.tune=num.Peng.tune,num.num=num.Peng.num,
                     cindex.train.tune=cindex.train.Peng.tune,
                     cindex.test.tune=cindex.test.Peng.tune.hard,
                     cindex.train.tune.adjust=cindex.train.Peng.tune.adjust,
                     cindex.test.tune.adjust=cindex.test.Peng.tune.adjust,
                     cindex.train.num=cindex.train.Peng.num,
                     cindex.test.num=cindex.test.Peng.num.hard,
                     cindex.train.num.adjust=cindex.train.Peng.num.adjust,
                     cindex.test.num.adjust=cindex.test.Peng.num.adjust,
                     F.train.tune=F.train.tune,
                     F.test.tune=F.test.tune,
                     F.train.num=F.train.num,
                     F.test.num=F.test.num,
                     R2.gene.test.tune=R2.gene.test.tune,
                     R2.gene.train.tune=R2.gene.train.tune,
                     R2.gene.test.num=R2.gene.test.num,
                     R2.gene.train.num=R2.gene.train.num,
                     lambda.tune=lambda,
                     lambda.num=lambda.num,
                     train_time=time.GM,
                     tune_metric=tune_metric,
                     num=num)
        #return(GM.res)
        print("Finish GM method!")

        #------------------------------
        #SKM method
        #------------------------------
        wbounds<-seq(1.1,25,0.5)
        #-------------------------------------------------
        #Permutation tune wbound
        #-------------------------------------------------

        optimal_wbound<-KMeansSparseCluster.permute(G.train,K=K,wbounds=wbounds,silent = T)$bestw

        start_time=Sys.time()
        mod.tune<-KMeansSparseCluster(G.train,K=K,wbounds=optimal_wbound,silent = T)
        end_time=Sys.time()
        time.SKM=as_hms(end_time-start_time)

        weight<-mod.tune[[1]]$ws
        select.feature<-as.numeric(weight!=0)
        num.skmeans.tune<-sum(select.feature==1)
        print(paste0("By permutation, wbound=",optimal_wbound," that selects  ",num.skmeans.tune," features."))

        ari.train.kmeans.tune<-adjustedRandIndex(mod.tune[[1]]$Cs,label.train)
        mod.train<-Rsquare(cluster =mod.tune[[1]]$Cs,Y = Y.train,X = X.train,G = t(G.train))
        R2.gene.train.tune<-mean(mod.train$gene[which(select.feature==1)])

        data1<-data.frame(cluster=mod.tune[[1]]$Cs,Y=Y.train,Y.ind=delta.train,x1=X.train[,1],x2=X.train[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.train.tune<-mod.anova$Chisq[2]

        centers<-cbind(apply(G.train[which(mod.tune[[1]]$Cs==1),,drop=F],2,mean),
                       apply(G.train[which(mod.tune[[1]]$Cs==2),,drop=F],2,mean),
                       apply(G.train[which(mod.tune[[1]]$Cs==3),,drop=F],2,mean))

        distance<-apply(G.test,1,function(x){
          distance1<-apply(centers,2,function(y){
            sum(weight*(x-y)^2)
          })
        })
        cluster.test<-apply(distance,2,which.min)
        ari.test.kmeans.tune<-adjustedRandIndex(cluster.test,label.test)

        mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test))
        R2.gene.test.tune<-mean(mod.test$gene[which(select.feature==1)])
        #Rsquare.outcome.test<-mod.test$outcome
        data1<-data.frame(cluster=cluster.test,Y=Y.test,Y.ind=delta.test,x1=X.test[,1],x2=X.test[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.test.tune<-mod.anova$Chisq[2]

        pred.train.kmeans<-rep(NA,length(Y.train))
        pred.test.kmeans<-rep(NA,length(Y.test))
        for(i in 1:K){
          data.temp<-data.frame(Y=Y.train[which(mod.tune[[1]]$Cs==i)],X.train[which(mod.tune[[1]]$Cs==i),])
          mod.lm<-lm(Y~.,data=data.temp)
          pred.train.kmeans[which(mod.tune[[1]]$Cs==i)]<-mod.lm$fitted.values
          pred.test.kmeans[which(cluster.test==i)]<-predict(mod.lm,newdata = as.data.frame(X.test[which(cluster.test==i),]))
        }
        cindex.train.kmeans.tune<-concordance.index(x=1/as.numeric(pred.train.kmeans), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
        cindex.test.kmeans.tune<-concordance.index(x=1/as.numeric(pred.test.kmeans), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        cindex.train.random<-c()
        cindex.test.random<-c()
        cluster.test1<-cluster.test
        cluster.train1<-mod.tune[[1]]$Cs
        for(ind in 1:50){
          cluster.test1<-sample(cluster.test1)
          cluster.train1<-sample(cluster.train1)
          pred.train.kmeans<-rep(NA,length(Y.train))
          pred.test.kmeans<-rep(NA,length(Y.test))
          for(i in 1:K){
            data.temp<-data.frame(Y=Y.train[which(cluster.train1==i)],X.train[which(cluster.train1==i),],delta=delta.train[which(cluster.train1==i)])
            mod.lm<-survreg(Surv(Y,delta)~x1+x2,dist="loglogistic",data=data.temp)
            pred.train.kmeans[which(cluster.train1==i)]<-predict(mod.lm,data.temp)
            pred.test.kmeans[which(cluster.test1==i)]<-predict(mod.lm,newdata = as.data.frame(X.test[which(cluster.test1==i),]))
          }
          cindex.train.random[ind]<-concordance.index(x=1/as.numeric(pred.train.kmeans), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
          cindex.test.random[ind]<-concordance.index(x=1/as.numeric(pred.test.kmeans), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        }

        cindex.train.kmeans.tune.adjust<-(cindex.train.kmeans.tune-mean(cindex.train.random))/(1-mean(cindex.train.random))
        cindex.test.kmeans.tune.adjust<-(cindex.test.kmeans.tune-mean(cindex.test.random))/(1-mean(cindex.test.random))

        jaccard.index.skmeans.tune<-Jaccard.index(true.feature,select.feature)

        #-------------------------------------------------
        #choose number of genes larger but cloest to 50
        #-------------------------------------------------

        mod.num<-KMeansSparseCluster(G.train,K=K,nstart=50)
        wbound.vector<-unlist(lapply(mod.num,function(x){x$wbound}))
        num.vector<-unlist(lapply(mod.num,function(x){sum(x$ws>0)}))
        # index<-min(which(num.vector>=50))
        index<-which.min(abs(num.vector-50))
        wbound.num=wbound.vector[index]
        mod.num<-mod.num[[index]]
        weight<-mod.num$ws
        select.feature<-as.numeric(weight!=0)
        num.skmeans.num<-sum(select.feature==1)
        print(paste0("Tuning wbound=",wbound.num," that selects  ",num.skmeans.num," features."))


        mod.train<-Rsquare(cluster =mod.num$Cs,Y = Y.train,X = X.train,G = t(G.train ))
        R2.gene.train.num<-mean(mod.train$gene[which(select.feature==1)])
        #mod.kmeans<-kmeans(G.train1,center=K,nstart=20)
        ari.train.kmeans.num<-adjustedRandIndex(mod.num$Cs,label.train)
        data1<-data.frame(cluster=mod.num$Cs,Y=Y.train,Y.ind=delta.train,x1=X.train[,1],x2=X.train[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.train.num<-mod.anova$Chisq[2]

        centers<-cbind(apply(G.train[which(mod.num$Cs==1),,drop=F],2,mean),
                       apply(G.train[which(mod.num$Cs==2),,drop=F],2,mean),
                       apply(G.train[which(mod.num$Cs==3),,drop=F],2,mean))

        distance<-apply(G.test,1,function(x){
          distance1<-apply(centers,2,function(y){
            sum(weight*(x-y)^2)
          })
        })
        cluster.test<-apply(distance,2,which.min)
        mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test))
        R2.gene.test.num<-mean(mod.test$gene[which(select.feature==1)])
        #Rsquare.outcome.test<-mod.test$outcome
        data1<-data.frame(cluster=cluster.test,Y=Y.test,Y.ind=delta.test,x1=X.test[,1],x2=X.test[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.test.num<-mod.anova$Chisq[2]

        ari.test.kmeans.num<-adjustedRandIndex(cluster.test,label.test)
        pred.train.kmeans<-rep(NA,length(Y.train))
        pred.test.kmeans<-rep(NA,length(Y.test))
        for(i in 1:K){
          data.temp<-data.frame(Y=Y.train[which(mod.num$Cs==i)],X.train[which(mod.num$Cs==i),])
          mod.lm<-lm(Y~.,data=data.temp)
          pred.train.kmeans[which(mod.num$Cs==i)]<-mod.lm$fitted.values
          pred.test.kmeans[which(cluster.test==i)]<-predict(mod.lm,newdata = as.data.frame(X.test[which(cluster.test==i),]))
        }
        cindex.train.kmeans.num<-concordance.index(x=1/as.numeric(pred.train.kmeans), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
        cindex.test.kmeans.num<-concordance.index(x=1/as.numeric(pred.test.kmeans), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        cindex.train.random<-c()
        cindex.test.random<-c()
        cluster.test1<-cluster.test
        cluster.train1<-mod.num$Cs
        for(ind in 1:50){
          cluster.test1<-sample(cluster.test1)
          cluster.train1<-sample(cluster.train1)
          pred.train.kmeans<-rep(NA,length(Y.train))
          pred.test.kmeans<-rep(NA,length(Y.test))
          for(i in 1:K){
            data.temp<-data.frame(Y=Y.train[which(cluster.train1==i)],X.train[which(cluster.train1==i),],delta=delta.train[which(cluster.train1==i)])
            mod.lm<-survreg(Surv(Y,delta)~x1+x2,dist="loglogistic",data=data.temp)
            pred.train.kmeans[which(cluster.train1==i)]<-predict(mod.lm,data.temp)
            pred.test.kmeans[which(cluster.test1==i)]<-predict(mod.lm,newdata = as.data.frame(X.test[which(cluster.test1==i),]))
          }
          cindex.train.random[ind]<-concordance.index(x=1/as.numeric(pred.train.kmeans), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
          cindex.test.random[ind]<-concordance.index(x=1/as.numeric(pred.test.kmeans), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        }

        cindex.train.kmeans.num.adjust<-(cindex.train.kmeans.num-mean(cindex.train.random))/(1-mean(cindex.train.random))
        cindex.test.kmeans.num.adjust<-(cindex.test.kmeans.num-mean(cindex.test.random))/(1-mean(cindex.test.random))

        jaccard.index.skmeans.num<-Jaccard.index(TrueFeature =true.feature,SelectedFeature1 =select.feature )

        #-------------------------------------------------
        SKM.res<-list(mod.tune=mod.tune,mod.num=mod.num,
                      ari.train.tune=ari.train.kmeans.tune,ari.train.num=ari.train.kmeans.num,
                      ari.test.tune=ari.test.kmeans.tune,ari.test.num=ari.test.kmeans.num,
                      num.tune=num.skmeans.tune,num.num=num.skmeans.num,
                      jaccard.index.tune=jaccard.index.skmeans.tune,jaccard.index.num=jaccard.index.skmeans.num,
                      cindex.train.tune=cindex.train.kmeans.tune,
                      cindex.test.tune=cindex.test.kmeans.tune,
                      cindex.train.tune.adjust=cindex.train.kmeans.tune.adjust,
                      cindex.test.tune.adjust=cindex.test.kmeans.tune.adjust,
                      cindex.train.num=cindex.train.kmeans.num,
                      cindex.test.num=cindex.test.kmeans.num,
                      cindex.train.num.adjust=cindex.train.kmeans.num.adjust,
                      cindex.test.num.adjust=cindex.test.kmeans.num.adjust,
                      F.train.tune=F.train.tune,
                      F.test.tune=F.test.tune,
                      F.train.num=F.train.num,
                      F.test.num=F.test.num,
                      R2.gene.test.tune=R2.gene.test.tune,
                      R2.gene.train.tune=R2.gene.train.tune,
                      R2.gene.test.num=R2.gene.test.num,
                      R2.gene.train.num=R2.gene.train.num,
                      wbound.tune=optimal_wbound,wbound.num=wbound.num,
                      train_time=time.SKM)

        print("Finish SKM method!")


        #------------------------------
        #PMBC method
        #------------------------------
        data.kmeans<-kmeans(as.matrix(G.train),K)
        c_size<-data.kmeans$size
        n=nrow(G.train)
        pi_int<-c_size/n
        miu_int<-sapply(1:K,function(x) apply(as.matrix(G.train[data.kmeans$cluster==x,]),2,mean))
        #dim(miu_int) #2000 x K
        sigma_c<-sapply(1:K,function(x) apply(as.matrix(G.train[data.kmeans$cluster==x,]),2,var))
        sum_squares<-sapply(1:K,function(x) c_size[x]*sigma_c[,x])
        sigma_int<-apply(sum_squares,1,sum)/n

        #lambda.PMBC=seq(20,100,2)
        lambda.PMBC=seq(8,70,2)
        #-------------------------------------------------
        #Tune lambda by BIC
        #-------------------------------------------------
        BIC_vector=rep(NA,length(lambda.PMBC))
        for(ind_lambda in 1:length(lambda.PMBC)) {
          lambda=lambda.PMBC[ind_lambda]
          BIC<-tryCatch({
            BIC<-em_mbc(t(G.train), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                        lambda=lambda, K=K, no_init=1, max_iter = 200)$BIC
          }, error=function(e){
            BIC=Inf
            return(BIC)
          })
          BIC_vector[ind_lambda]=BIC
        }

        lambda.tune=lambda.PMBC[which.min(BIC_vector)]

        start_time=Sys.time()
        mod.tune<-em_mbc(t(G.train), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                         lambda=lambda.tune, K=K, no_init=1, max_iter = 200)
        end_time=Sys.time()
        time.PMBC=as_hms(end_time-start_time)

        num.PMBC.tune=length(mod.tune$selected.index)
        print(paste0("By BIC, lambda=",lambda.tune," that selects  ",num.PMBC.tune," features."))
        select.feature=rep(0,ncol(G.train))
        select.feature[mod.tune$selected.index]=1

        pi_est_train<-mod.tune$optimal_result$pi
        mu_est_train<-mod.tune$optimal_result$mu
        sigma_est_train<-mod.tune$optimal_result$sigma
        grp_assign_train<-apply(mod.tune$optimal_result$z,1,which.max)
        cluster.train=grp_assign_train
        ari.train.PMBC.tune<-adjustedRandIndex(cluster.train,label.train)

        data1<-data.frame(cluster=cluster.train,Y=Y.train,Y.ind=delta.train,x1=X.train[,1],x2=X.train[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.train.tune<-mod.anova$Chisq[2]
        mod.train<-Rsquare(cluster =cluster.train,Y = Y.train,X = X.train,G = t(G.train))
        R2.gene.train.tune<-mean(mod.train$gene[which(select.feature==1)])

        mult_pdf=t(apply(G.test,1,function(x) {
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
        ari.test.PMBC.tune<-adjustedRandIndex(cluster.test,label.test)
        data1<-data.frame(cluster=cluster.test,Y=Y.test,Y.ind=delta.test,x1=X.test[,1],x2=X.test[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.test.tune<-mod.anova$Chisq[2]

        mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test))
        R2.gene.test.tune<-mean(mod.test$gene[which(select.feature==1)])


        pred.train.PMBC<-rep(NA,length(Y.train))
        pred.test.PMBC<-rep(NA,length(Y.test))
        for(i in 1:K){
          data.temp<-data.frame(Y=Y.train[which(cluster.train==i)],X.train[which(cluster.train==i),])
          mod.lm<-lm(Y~.,data=data.temp)
          pred.train.PMBC[which(cluster.train==i)]<-mod.lm$fitted.values
          pred.test.PMBC[which(cluster.test==i)]<-predict(mod.lm,newdata = as.data.frame(X.test[which(cluster.test==i),]))
        }

        cindex.train.tune<-concordance.index(x=1/as.numeric(pred.train.PMBC), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
        cindex.test.tune<-concordance.index(x=1/as.numeric(pred.test.PMBC), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        cindex.train.random<-c()
        cindex.test.random<-c()
        cluster.test1<-cluster.test
        cluster.train1<-cluster.train
        for(ind in 1:50){
          cluster.test1<-sample(cluster.test1)
          cluster.train1<-sample(cluster.train1)
          pred.train.PMBC<-rep(NA,length(Y.train))
          pred.test.PMBC<-rep(NA,length(Y.test))
          for(i in 1:K){
            data.temp<-data.frame(Y=Y.train[which(cluster.train1==i)],X.train[which(cluster.train1==i),],delta=delta.train[which(cluster.train1==i)])
            mod.lm<-survreg(Surv(Y,delta)~x1+x2,dist="loglogistic",data=data.temp)
            pred.train.PMBC[which(cluster.train1==i)]<-predict(mod.lm,data.temp)
            pred.test.PMBC[which(cluster.test1==i)]<-predict(mod.lm,newdata = as.data.frame(X.test[which(cluster.test1==i),]))
          }
          cindex.train.random[ind]<-concordance.index(x=1/as.numeric(pred.train.PMBC), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
          cindex.test.random[ind]<-concordance.index(x=1/as.numeric(pred.test.PMBC), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        }

        cindex.train.tune.adjust<-(cindex.train.tune-mean(cindex.train.random))/(1-mean(cindex.train.random))
        cindex.test.tune.adjust<-(cindex.test.tune-mean(cindex.test.random))/(1-mean(cindex.test.random))

        jaccard.index.PMBC.tune<-Jaccard.index(true.feature,select.feature)

        #-------------------------------------------------
        #choose number of genes cloest to 50
        #-------------------------------------------------

        selected.all=rep(NA, length(lambda.PMBC))
        for(m in 1:length(lambda.PMBC)){
          lambda=lambda.PMBC[m]
          ngenes<-tryCatch({
            ngenes<-length(em_mbc(t(G.train), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                                  lambda=lambda, K=K, no_init=1, max_iter = 200)$selected.index)
          }, error=function(e){
            ngenes=Inf
            return(ngenes)
          })
          selected.all[m]<-ngenes
        }

        lambda.num=lambda.PMBC[which.min(abs(selected.all-50))]

        mod.num<-em_mbc(t(G.train), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                        lambda=lambda.num, K=K, no_init=1, max_iter = 200)

        num.PMBC.num=length(mod.num$selected.index)
        print(paste0("Tuning lambda=",lambda.num," that selects  ",num.PMBC.num," features."))

        pi_est_train<-mod.num$optimal_result$pi
        mu_est_train<-mod.num$optimal_result$mu
        sigma_est_train<-mod.num$optimal_result$sigma
        grp_assign_train<-apply(mod.num$optimal_result$z,1,which.max)
        cluster.train=grp_assign_train
        ari.train.PMBC.num<-adjustedRandIndex(cluster.train,label.train)
        mod.train<-Rsquare(cluster =cluster.train,Y = Y.train,X = X.train,G = t(G.train ))
        R2.gene.train.num<-mean(mod.train$gene[mod.num$selected.index])

        data1<-data.frame(cluster=cluster.train,Y=Y.train,Y.ind=delta.train,x1=X.train[,1],x2=X.train[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.train.num<-mod.anova$Chisq[2]


        mult_pdf=t(apply(G.test,1,function(x) {
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
        ari.test.PMBC.num<-adjustedRandIndex(cluster.test,label.test)
        mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test))
        R2.gene.test.num<-mean(mod.test$gene[mod.num$selected.index])
        data1<-data.frame(cluster=cluster.test,Y=Y.test,Y.ind=delta.test,x1=X.test[,1],x2=X.test[,2])
        mod.cox1<-coxph(Surv(Y, Y.ind) ~ as.factor(cluster)+x1+x2, data = data1)
        mod.cox2<-coxph(Surv(Y, Y.ind) ~ x1+x2, data = data1)
        mod.anova<-anova(mod.cox1,mod.cox2)
        F.test.num<-mod.anova$Chisq[2]


        pred.train.PMBC<-rep(NA,length(Y.train))
        pred.test.PMBC<-rep(NA,length(Y.test))
        for(i in 1:K){
          data.temp<-data.frame(Y=Y.train[which(cluster.train==i)],X.train[which(cluster.train==i),])
          mod.lm<-lm(Y~.,data=data.temp)
          pred.train.PMBC[which(cluster.train==i)]<-mod.lm$fitted.values
          pred.test.PMBC[which(cluster.test==i)]<-predict(mod.lm,newdata = as.data.frame(X.test[which(cluster.test==i),]))
        }

        cindex.train.num<-concordance.index(x=1/as.numeric(pred.train.PMBC), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
        cindex.test.num<-concordance.index(x=1/as.numeric(pred.test.PMBC), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        cindex.train.random<-c()
        cindex.test.random<-c()
        cluster.test1<-cluster.test
        cluster.train1<-cluster.train
        for(ind in 1:50){
          cluster.test1<-sample(cluster.test1)
          cluster.train1<-sample(cluster.train1)
          pred.train.PMBC<-rep(NA,length(Y.train))
          pred.test.PMBC<-rep(NA,length(Y.test))
          for(i in 1:K){
            data.temp<-data.frame(Y=Y.train[which(cluster.train1==i)],X.train[which(cluster.train1==i),],delta=delta.train[which(cluster.train1==i)])
            mod.lm<-survreg(Surv(Y,delta)~x1+x2,dist="loglogistic",data=data.temp)
            pred.train.PMBC[which(cluster.train1==i)]<-predict(mod.lm,data.temp)
            pred.test.PMBC[which(cluster.test1==i)]<-predict(mod.lm,newdata = as.data.frame(X.test[which(cluster.test1==i),]))
          }
          cindex.train.random[ind]<-concordance.index(x=1/as.numeric(pred.train.PMBC), surv.time=Y.train, surv.event=delta.train, method="noether")$c.index
          cindex.test.random[ind]<-concordance.index(x=1/as.numeric(pred.test.PMBC), surv.time=Y.test, surv.event=delta.test, method="noether")$c.index

        }

        cindex.train.num.adjust<-(cindex.train.num-mean(cindex.train.random))/(1-mean(cindex.train.random))
        cindex.test.num.adjust<-(cindex.test.num-mean(cindex.test.random))/(1-mean(cindex.test.random))

        select.feature=rep(0,ncol(G.train))
        select.feature[mod.num$selected.index]=1
        jaccard.index.PMBC.num<-Jaccard.index(true.feature,select.feature)

        #-------------------------------------------------
        PMBC.res<-list(mod.tune=mod.tune,mod.num=mod.num,
                       ari.train.tune=ari.train.PMBC.tune,ari.train.num=ari.train.PMBC.num,
                       ari.test.tune=ari.test.PMBC.tune,ari.test.num=ari.test.PMBC.num,
                       num.tune=num.PMBC.tune,num.num=num.PMBC.num,
                       jaccard.index.tune=jaccard.index.PMBC.tune,jaccard.index.num=jaccard.index.PMBC.num,
                       cindex.train.tune=cindex.train.tune,
                       cindex.test.tune=cindex.test.tune,
                       cindex.train.tune.adjust=cindex.train.tune.adjust,
                       cindex.test.tune.adjust=cindex.test.tune.adjust,
                       cindex.train.num=cindex.train.num,
                       cindex.test.num=cindex.test.num,
                       cindex.train.num.adjust=cindex.train.num.adjust,
                       cindex.test.num.adjust=cindex.test.num.adjust,
                       F.train.tune=F.train.tune,
                       F.test.tune=F.test.tune,
                       F.train.num=F.train.num,
                       F.test.num=F.test.num,
                       R2.gene.test.tune=R2.gene.test.tune,
                       R2.gene.train.tune=R2.gene.train.tune,
                       R2.gene.test.num=R2.gene.test.num,
                       R2.gene.train.num=R2.gene.train.num,
                       lambda.tune=lambda.tune,lambda.num=lambda.num,
                       train_time=time.PMBC)
        return(list(WJL.res=WJL.res,GM.res=GM.res,SKM.res=SKM.res,PMBC.res=PMBC.res,index.test=index.test))
      },mc.cores=150)
      if (SC) {
        file.res<-paste("SimII_RES_c1=",c1,"_mu=",mu,".SC.rds",sep="")
      }
      else {
        file.res<-paste("SimII_RES_c1=",c1,"_mu=",mu,".rds",sep="")
      }
      saveRDS(RES,file.res)
    }
  }
}

res.SC=resSurv.gen(K=3,SC=TRUE,num.sim=50,mu_vector=c(0.9,1.2,1.5,1.8),c1_vector_all=c(2))
res=resSurv.gen(K=3,SC=FALSE,num.sim=50,mu_vector=c(0.9,1.2,1.5,1.8),c1_vector_all=c(2))


library(gridExtra)
library(tidyr)
library(ggplot2)
summarySurv_gen=function(SC=TRUE,tune=FALSE,mu_vector=c(0.9,1.2,1.5,1.8),c1_vector_all=c(2),figure.name){
  ari.train.yujia.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.yujia.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.yujia.matrix.hard.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.yujia.matrix.hard.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.yujia.matrix.hard.num.adjust<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.yujia.matrix.hard.num.adjust<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))

  Jaccard.yujia.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  num.yujia.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.train<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.test<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))

  ari.train.yujia.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.yujia.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.yujia.matrix.num.hard.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.yujia.matrix.num.hard.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.yujia.matrix.num.hard.adjust.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.yujia.matrix.num.hard.adjust.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.yujia.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.train.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.test.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))

  w.yujia.matrix<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  w.yujia.matrix.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))

  ari.train.Peng.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.Peng.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.Peng.matrix.num.hard<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.Peng.matrix.num.hard<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.Peng.matrix.num.hard.adjust<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.Peng.matrix.num.hard.adjust<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.Peng.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  num.Peng.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.train.Peng<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.Peng<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.test.Peng<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.Peng<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))

  ari.train.Peng.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.Peng.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.Peng.matrix.num.hard.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.Peng.matrix.num.hard.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.Peng.matrix.num.hard.adjust.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.Peng.matrix.num.hard.adjust.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.Peng.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.train.Peng.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.Peng.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.test.Peng.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.Peng.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))

  ari.train.skmeans.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.skmeans.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.skmeans.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.skmeans.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.skmeans.matrix.num.adjust<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.skmeans.matrix.num.adjust<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.skmeans.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  num.skmeans.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.train.skmeans<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.skmeans<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.test.skmeans<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.skmeans<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))

  ari.train.skmeans.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.skmeans.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.skmeans.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.skmeans.matrix.num.adjust.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.skmeans.matrix.num.adjust.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.skmeans.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.skmeans.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.train.skmeans.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.skmeans.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.test.skmeans.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.skmeans.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))

  ari.train.PMBC.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.PMBC.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.PMBC.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.PMBC.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.PMBC.matrix.num.adjust<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.PMBC.matrix.num.adjust<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.PMBC.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  num.PMBC.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.train.PMBC<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.PMBC<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.test.PMBC<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.PMBC<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))

  ari.train.PMBC.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.PMBC.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.PMBC.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.PMBC.matrix.num.adjust.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.train.PMBC.matrix.num.adjust.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  cindex.test.PMBC.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.PMBC.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.train.PMBC.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.PMBC.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  F.outcome.matrix.test.PMBC.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.PMBC.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))


  for(ind_c1 in 1:length(c1_vector_all)) {
    c1<-c1_vector_all[ind_c1]
    for (ind_mu in 1:length(mu_vector)) {
      mu=mu_vector[ind_mu]
      if (SC) {
        file.name= paste("SimII_RES_c1=",c1,"_mu=",mu,".SC.rds",sep="")
      }
      else {
        file.name= paste("SimII_RES_c1=",c1,"_mu=",mu,".rds",sep="")
      }
      RES=readRDS(file.name)
      RES=RES[sapply(RES,length)==5]
      if (tune) {
        ari.train.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$ari.train.tune})),na.rm=T)
        ari.test.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$ari.test.tune})),na.rm=T)
        cindex.train.yujia.matrix.hard.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$cindex.train.tune})),na.rm=T)
        cindex.test.yujia.matrix.hard.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$cindex.test.tune})),na.rm=T)
        cindex.train.yujia.matrix.hard.num.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$cindex.train.tune.adjust})),na.rm=T)
        cindex.test.yujia.matrix.hard.num.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$cindex.test.tune.adjust})),na.rm=T)
        Jaccard.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$jaccard.index.tune})),na.rm=T)
        num.yujia.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$WJL.res$num.tune})))
        F.outcome.matrix.train[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$F.train.tune})),na.rm=T)
        F.outcome.matrix.test[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$F.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.test.tune})),na.rm=T)

        ari.train.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$ari.train.tune})),na.rm=T)
        ari.test.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$ari.test.tune})),na.rm=T)
        cindex.train.yujia.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$cindex.train.tune})),na.rm=T)
        cindex.test.yujia.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$cindex.test.tune})),na.rm=T)
        cindex.train.yujia.matrix.num.hard.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$cindex.train.tune.adjust})),na.rm=T)
        cindex.test.yujia.matrix.num.hard.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$cindex.test.tune.adjust})),na.rm=T)
        Jaccard.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$jaccard.index.tune})),na.rm=T)
        F.outcome.matrix.train.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$F.train.tune})),na.rm=T)
        F.outcome.matrix.test.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$F.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$Rsquare.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.test.tune})),na.rm=T)

        w.yujia.matrix[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$w0.est})),na.rm=T)
        w.yujia.matrix.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$w0.est})),na.rm=T)/sqrt(50)

        ari.train.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$ari.train.tune})),na.rm=T)
        ari.test.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$ari.test.tune})),na.rm=T)
        cindex.train.Peng.matrix.num.hard[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$cindex.train.tune})),na.rm=T)
        cindex.test.Peng.matrix.num.hard[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$cindex.test.tune})),na.rm=T)
        cindex.train.Peng.matrix.num.hard.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$cindex.train.tune.adjust})),na.rm=T)
        cindex.test.Peng.matrix.num.hard.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$cindex.test.tune.adjust})),na.rm=T)
        Jaccard.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$jaccard.index.tune})),na.rm=T)
        num.Peng.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$GM.res$num.tune})))
        F.outcome.matrix.train.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$F.train.tune})),na.rm=T)
        F.outcome.matrix.test.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$F.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$R2.gene.test.tune})),na.rm=T)

        ari.train.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$ari.train.tune})),na.rm=T)
        ari.test.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$ari.test.tune})),na.rm=T)
        cindex.train.Peng.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$cindex.train.tune})),na.rm=T)
        cindex.test.Peng.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$cindex.test.tune})),na.rm=T)
        cindex.train.Peng.matrix.num.hard.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$cindex.train.tune.adjust})),na.rm=T)
        cindex.test.Peng.matrix.num.hard.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$cindex.test.tune.adjust})),na.rm=T)
        Jaccard.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$jaccard.index.tune})),na.rm=T)
        F.outcome.matrix.train.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$F.train.tune})),na.rm=T)
        F.outcome.matrix.test.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$F.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$R2.gene.test.tune})),na.rm=T)

        ari.train.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$ari.train.tune})),na.rm=T)
        ari.test.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$ari.test.tune})),na.rm=T)
        cindex.train.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$cindex.train.tune})),na.rm=T)
        cindex.test.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$cindex.test.tune})),na.rm=T)
        cindex.train.skmeans.matrix.num.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$cindex.train.tune.adjust})),na.rm=T)
        cindex.test.skmeans.matrix.num.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$cindex.test.tune.adjust})),na.rm=T)
        Jaccard.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$jaccard.index.tune})),na.rm=T)
        num.skmeans.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$SKM.res$num.tune})))
        F.outcome.matrix.train.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$F.train.tune})),na.rm=T)
        F.outcome.matrix.test.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$F.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.test.tune})),na.rm=T)

        ari.train.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$ari.train.tune})),na.rm=T)
        ari.test.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$ari.test.tune})),na.rm=T)
        cindex.train.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$cindex.train.tune})),na.rm=T)
        cindex.test.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$cindex.test.tune})),na.rm=T)
        cindex.train.skmeans.matrix.num.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$cindex.train.tune.adjust})),na.rm=T)
        cindex.test.skmeans.matrix.num.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$cindex.test.tune.adjust})),na.rm=T)
        Jaccard.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$jaccard.index.tune})),na.rm=T)
        F.outcome.matrix.train.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$F.train.tune})),na.rm=T)
        F.outcome.matrix.test.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$F.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.test.tune})),na.rm=T)

        ari.train.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$ari.train.tune})),na.rm=T)
        ari.test.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$ari.test.tune})),na.rm=T)
        cindex.train.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$cindex.train.tune})),na.rm=T)
        cindex.test.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$cindex.test.tune})),na.rm=T)
        cindex.train.PMBC.matrix.num.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$cindex.train.tune.adjust})),na.rm=T)
        cindex.test.PMBC.matrix.num.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$cindex.test.tune.adjust})),na.rm=T)
        Jaccard.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$jaccard.index.tune})),na.rm=T)
        num.PMBC.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$PMBC.res$num.tune})))
        F.outcome.matrix.train.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$F.train.tune})),na.rm=T)
        F.outcome.matrix.test.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$F.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.test.tune})),na.rm=T)

        ari.train.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$ari.train.tune})),na.rm=T)
        ari.test.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$ari.test.tune})),na.rm=T)
        cindex.train.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$cindex.train.tune})),na.rm=T)
        cindex.test.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$cindex.test.tune})),na.rm=T)
        cindex.train.PMBC.matrix.num.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$cindex.train.tune.adjust})),na.rm=T)
        cindex.test.PMBC.matrix.num.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$cindex.test.tune.adjust})),na.rm=T)
        Jaccard.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$jaccard.index.tune})),na.rm=T)
        F.outcome.matrix.train.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$F.train.tune})),na.rm=T)
        F.outcome.matrix.test.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$F.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.test.tune})),na.rm=T)

      }
      else {
        ari.train.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$ari.train.num})),na.rm=T)
        ari.test.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$ari.test.num})),na.rm=T)
        cindex.train.yujia.matrix.hard.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$cindex.train.num})),na.rm=T)
        cindex.test.yujia.matrix.hard.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$cindex.test.num})),na.rm=T)
        cindex.train.yujia.matrix.hard.num.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$cindex.train.num.adjust})),na.rm=T)
        cindex.test.yujia.matrix.hard.num.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$cindex.test.num.adjust})),na.rm=T)
        Jaccard.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$jaccard.index.num})),na.rm=T)
        num.yujia.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$WJL.res$num.num})))
        F.outcome.matrix.train[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$F.train.num})),na.rm=T)
        F.outcome.matrix.test[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$F.test.num})),na.rm=T)
        Rsquare.gene.matrix.train[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.train.num})),na.rm=T)
        Rsquare.gene.matrix.test[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.test.num})),na.rm=T)

        ari.train.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$ari.train.num})),na.rm=T)
        ari.test.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$ari.test.num})),na.rm=T)
        cindex.train.yujia.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$cindex.train.num})),na.rm=T)
        cindex.test.yujia.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$cindex.test.num})),na.rm=T)
        cindex.train.yujia.matrix.num.hard.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$cindex.train.num.adjust})),na.rm=T)
        cindex.test.yujia.matrix.num.hard.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$cindex.test.num.adjust})),na.rm=T)
        Jaccard.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$jaccard.index.num})),na.rm=T)
        F.outcome.matrix.train.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$F.train.num})),na.rm=T)
        F.outcome.matrix.test.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$F.test.num})),na.rm=T)
        Rsquare.gene.matrix.train.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$Rsquare.gene.train.num})),na.rm=T)
        Rsquare.gene.matrix.test.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.test.num})),na.rm=T)

        w.yujia.matrix[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$w0.est})),na.rm=T)
        w.yujia.matrix.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$w0.est})),na.rm=T)/sqrt(50)

        ari.train.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$ari.train.num})),na.rm=T)
        ari.test.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$ari.test.num})),na.rm=T)
        cindex.train.Peng.matrix.num.hard[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$cindex.train.num})),na.rm=T)
        cindex.test.Peng.matrix.num.hard[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$cindex.test.num})),na.rm=T)
        cindex.train.Peng.matrix.num.hard.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$cindex.train.num.adjust})),na.rm=T)
        cindex.test.Peng.matrix.num.hard.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$cindex.test.num.adjust})),na.rm=T)
        Jaccard.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$jaccard.index.num})),na.rm=T)
        num.Peng.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$GM.res$num.num})),na.rm=T)
        F.outcome.matrix.train.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$F.train.num})),na.rm=T)
        F.outcome.matrix.test.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$F.test.num})),na.rm=T)
        Rsquare.gene.matrix.train.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$R2.gene.train.num})),na.rm=T)
        Rsquare.gene.matrix.test.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$R2.gene.test.num})),na.rm=T)

        ari.train.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$ari.train.num})),na.rm=T)
        ari.test.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$ari.test.num})),na.rm=T)
        cindex.train.Peng.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$cindex.train.num})),na.rm=T)
        cindex.test.Peng.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$cindex.test.num})),na.rm=T)
        cindex.train.Peng.matrix.num.hard.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$cindex.train.num.adjust})),na.rm=T)
        cindex.test.Peng.matrix.num.hard.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$cindex.test.num.adjust})),na.rm=T)
        Jaccard.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$jaccard.index.num})),na.rm=T)
        F.outcome.matrix.train.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$F.train.num})),na.rm=T)
        F.outcome.matrix.test.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$F.test.num})),na.rm=T)
        Rsquare.gene.matrix.train.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$R2.gene.train.num})),na.rm=T)
        Rsquare.gene.matrix.test.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$R2.gene.test.num})),na.rm=T)

        ari.train.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$ari.train.num})),na.rm=T)
        ari.test.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$ari.test.num})),na.rm=T)
        cindex.train.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$cindex.train.num})),na.rm=T)
        cindex.test.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$cindex.test.num})),na.rm=T)
        cindex.train.skmeans.matrix.num.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$cindex.train.num.adjust})),na.rm=T)
        cindex.test.skmeans.matrix.num.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$cindex.test.num.adjust})),na.rm=T)
        Jaccard.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$jaccard.index.num})),na.rm=T)
        num.skmeans.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$SKM.res$num.num})))
        F.outcome.matrix.train.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$F.train.num})),na.rm=T)
        F.outcome.matrix.test.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$F.test.num})),na.rm=T)
        Rsquare.gene.matrix.train.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.train.num})),na.rm=T)
        Rsquare.gene.matrix.test.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.test.num})),na.rm=T)

        ari.train.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$ari.train.num})),na.rm=T)
        ari.test.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$ari.test.num})),na.rm=T)
        cindex.train.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$cindex.train.num})),na.rm=T)
        cindex.test.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$cindex.test.num})),na.rm=T)
        cindex.train.skmeans.matrix.num.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$cindex.train.num.adjust})),na.rm=T)
        cindex.test.skmeans.matrix.num.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$cindex.test.num.adjust})),na.rm=T)
        Jaccard.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$jaccard.index.num})),na.rm=T)
        F.outcome.matrix.train.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$F.train.num})),na.rm=T)
        F.outcome.matrix.test.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$F.test.num})),na.rm=T)
        Rsquare.gene.matrix.train.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.train.num})),na.rm=T)
        Rsquare.gene.matrix.test.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.test.num})),na.rm=T)

        ari.train.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$ari.train.num})),na.rm=T)
        ari.test.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$ari.test.num})),na.rm=T)
        cindex.train.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$cindex.train.num})),na.rm=T)
        cindex.test.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$cindex.test.num})),na.rm=T)
        cindex.train.PMBC.matrix.num.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$cindex.train.num.adjust})),na.rm=T)
        cindex.test.PMBC.matrix.num.adjust[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$cindex.test.num.adjust})),na.rm=T)
        Jaccard.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$jaccard.index.num})),na.rm=T)
        num.PMBC.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$PMBC.res$num.num})))
        F.outcome.matrix.train.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$F.train.num})),na.rm=T)
        F.outcome.matrix.test.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$F.test.num})),na.rm=T)
        Rsquare.gene.matrix.train.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.train.num})),na.rm=T)
        Rsquare.gene.matrix.test.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.test.num})),na.rm=T)

        ari.train.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$ari.train.num})),na.rm=T)
        ari.test.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$ari.test.num})),na.rm=T)
        cindex.train.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$cindex.train.num})),na.rm=T)
        cindex.test.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$cindex.test.num})),na.rm=T)
        cindex.train.PMBC.matrix.num.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$cindex.train.num.adjust})),na.rm=T)
        cindex.test.PMBC.matrix.num.adjust.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$cindex.test.num.adjust})),na.rm=T)
        Jaccard.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$jaccard.index.num})),na.rm=T)
        F.outcome.matrix.train.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$F.train.num})),na.rm=T)
        F.outcome.matrix.test.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$F.test.num})),na.rm=T)
        Rsquare.gene.matrix.train.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.train.num})),na.rm=T)
        Rsquare.gene.matrix.test.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.test.num})),na.rm=T)

      }

    }

    ari.train1<-rbind(mu_vector,
                      ari.train.yujia.matrix.num[ind_c1,],ari.train.Peng.matrix.num[ind_c1,],
                      ari.train.skmeans.matrix.num[ind_c1,],ari.train.PMBC.matrix.num[ind_c1,])
    ari.test1<-rbind(mu_vector,
                     ari.test.yujia.matrix.num[ind_c1,],ari.test.Peng.matrix.num[ind_c1,],
                     ari.test.skmeans.matrix.num[ind_c1,],ari.test.PMBC.matrix.num[ind_c1,])

    cindex.train1<-rbind(mu_vector,
                         cindex.train.yujia.matrix.hard.num.adjust[ind_c1,],cindex.train.Peng.matrix.num.hard.adjust[ind_c1,],
                         cindex.train.skmeans.matrix.num.adjust[ind_c1,],cindex.train.PMBC.matrix.num.adjust[ind_c1,])
    cindex.test1<-rbind(mu_vector,
                        cindex.test.yujia.matrix.hard.num.adjust[ind_c1,],cindex.test.Peng.matrix.num.hard.adjust[ind_c1,],
                        cindex.test.skmeans.matrix.num.adjust[ind_c1,],cindex.test.PMBC.matrix.num.adjust[ind_c1,])
    jaccard1<-rbind(mu_vector,
                    Jaccard.yujia.matrix.num[ind_c1,],Jaccard.Peng.matrix.num[ind_c1,],
                    Jaccard.skmeans.matrix.num[ind_c1,],Jaccard.PMBC.matrix.num[ind_c1,])
    # num1<-rbind(mu_vector,
    #             num.yujia.matrix.num[ind_c1,],num.Peng.matrix.num[ind_c1,],
    #             num.skmeans.matrix.num[ind_c1,],num.PMBC.matrix.num[ind_c1,])

    F.outcome.train1<-rbind(mu_vector,
                            F.outcome.matrix.train[ind_c1,],F.outcome.matrix.train.Peng[ind_c1,],
                            F.outcome.matrix.train.skmeans[ind_c1,],F.outcome.matrix.train.PMBC[ind_c1,])

    Rsquare.gene.train1<-rbind(mu_vector,
                               Rsquare.gene.matrix.train[ind_c1,],Rsquare.gene.matrix.train.Peng[ind_c1,],
                               Rsquare.gene.matrix.train.skmeans[ind_c1,],Rsquare.gene.matrix.train.PMBC[ind_c1,])

    F.outcome.test1<-rbind(mu_vector,
                           F.outcome.matrix.test[ind_c1,],F.outcome.matrix.test.Peng[ind_c1,],
                           F.outcome.matrix.test.skmeans[ind_c1,],F.outcome.matrix.test.PMBC[ind_c1,])

    Rsquare.gene.test1<-rbind(mu_vector,
                              Rsquare.gene.matrix.test[ind_c1,],Rsquare.gene.matrix.test.Peng[ind_c1,],
                              Rsquare.gene.matrix.test.skmeans[ind_c1,],Rsquare.gene.matrix.test.PMBC[ind_c1,])

    ari.train1.se<-rbind(mu_vector,
                         ari.train.yujia.matrix.num.se[ind_c1,],ari.train.Peng.matrix.num.se[ind_c1,],
                         ari.train.skmeans.matrix.num.se[ind_c1,],ari.train.PMBC.matrix.num.se[ind_c1,])
    ari.test1.se<-rbind(mu_vector,
                        ari.test.yujia.matrix.num.se[ind_c1,],ari.test.Peng.matrix.num.se[ind_c1,],
                        ari.test.skmeans.matrix.num.se[ind_c1,],ari.test.PMBC.matrix.num.se[ind_c1,])

    cindex.train1.se<-rbind(mu_vector,
                            cindex.train.yujia.matrix.num.hard.adjust.se[ind_c1,],cindex.train.Peng.matrix.num.hard.adjust.se[ind_c1,],
                            cindex.train.skmeans.matrix.num.adjust.se[ind_c1,],cindex.train.PMBC.matrix.num.adjust.se[ind_c1,])
    cindex.test1.se<-rbind(mu_vector,
                           cindex.test.yujia.matrix.num.hard.adjust.se[ind_c1,],cindex.test.Peng.matrix.num.hard.adjust.se[ind_c1,],
                           cindex.test.skmeans.matrix.num.adjust.se[ind_c1,],cindex.test.PMBC.matrix.num.adjust.se[ind_c1,])
    jaccard1.se<-rbind(mu_vector,
                       Jaccard.yujia.matrix.num.se[ind_c1,],Jaccard.Peng.matrix.num.se[ind_c1,],
                       Jaccard.skmeans.matrix.num.se[ind_c1,],Jaccard.PMBC.matrix.num.se[ind_c1,])

    F.outcome.train1.se<-rbind(mu_vector,
                               F.outcome.matrix.train.se[ind_c1,],F.outcome.matrix.train.Peng.se[ind_c1,],
                               F.outcome.matrix.train.skmeans.se[ind_c1,],F.outcome.matrix.train.PMBC.se[ind_c1,])

    Rsquare.gene.train1.se<-rbind(mu_vector,
                                  Rsquare.gene.matrix.train.se[ind_c1,],Rsquare.gene.matrix.train.Peng.se[ind_c1,],
                                  Rsquare.gene.matrix.train.skmeans.se[ind_c1,],Rsquare.gene.matrix.train.PMBC.se[ind_c1,])

    F.outcome.test1.se<-rbind(mu_vector,
                              F.outcome.matrix.test.se[ind_c1,],F.outcome.matrix.test.Peng.se[ind_c1,],
                              F.outcome.matrix.test.skmeans.se[ind_c1,],F.outcome.matrix.test.PMBC.se[ind_c1,])

    Rsquare.gene.test1.se<-rbind(mu_vector,
                                 Rsquare.gene.matrix.test.se[ind_c1,],Rsquare.gene.matrix.test.Peng.se[ind_c1,],
                                 Rsquare.gene.matrix.test.skmeans.se[ind_c1,],Rsquare.gene.matrix.test.PMBC.se[ind_c1,])


    rownames(ari.train1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(ari.test1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(cindex.train1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(cindex.test1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(jaccard1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(F.outcome.train1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rsquare.gene.train1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(F.outcome.test1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rsquare.gene.test1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")

    rownames(ari.train1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(ari.test1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(cindex.train1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(cindex.test1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(jaccard1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(F.outcome.train1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rsquare.gene.train1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(F.outcome.test1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rsquare.gene.test1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")

    ari.train1.se<-pivot_longer(data = as.data.frame(t(ari.train1.se)),cols=!Mu,names_to="method",values_to="value")
    ari.test1.se<-pivot_longer(data = as.data.frame(t(ari.test1.se)),cols=!Mu,names_to="method",values_to="value")
    cindex.train1.se<-pivot_longer(data = as.data.frame(t(cindex.train1.se)),cols=!Mu,names_to="method",values_to="value")
    cindex.test1.se<-pivot_longer(data = as.data.frame(t(cindex.test1.se)),cols=!Mu,names_to="method",values_to="value")
    jaccard1.se<-pivot_longer(data = as.data.frame(t(jaccard1.se)),cols=!Mu,names_to="method",values_to="value")
    F.outcome.train1.se<-pivot_longer(data = as.data.frame(t(F.outcome.train1.se)),cols=!Mu,names_to="method",values_to="value")
    Rsquare.gene.train1.se<-pivot_longer(data = as.data.frame(t(Rsquare.gene.train1.se)),cols=!Mu,names_to="method",values_to="value")
    F.outcome.test1.se<-pivot_longer(data = as.data.frame(t(F.outcome.test1.se)),cols=!Mu,names_to="method",values_to="value")
    Rsquare.gene.test1.se<-pivot_longer(data = as.data.frame(t(Rsquare.gene.test1.se)),cols=!Mu,names_to="method",values_to="value")


    ari.train1<-pivot_longer(data = as.data.frame(t(ari.train1)),cols=!Mu,names_to="method",values_to="value")
    ari.test1<-pivot_longer(data = as.data.frame(t(ari.test1)),cols=!Mu,names_to="method",values_to="value")
    cindex.train1<-pivot_longer(data = as.data.frame(t(cindex.train1)),cols=!Mu,names_to="method",values_to="value")
    cindex.test1<-pivot_longer(data = as.data.frame(t(cindex.test1)),cols=!Mu,names_to="method",values_to="value")
    jaccard1<-pivot_longer(data = as.data.frame(t(jaccard1)),cols=!Mu,names_to="method",values_to="value")
    F.outcome.train1<-pivot_longer(data = as.data.frame(t(F.outcome.train1)),cols=!Mu,names_to="method",values_to="value")
    Rsquare.gene.train1<-pivot_longer(data = as.data.frame(t(Rsquare.gene.train1)),cols=!Mu,names_to="method",values_to="value")
    F.outcome.test1<-pivot_longer(data = as.data.frame(t(F.outcome.test1)),cols=!Mu,names_to="method",values_to="value")
    Rsquare.gene.test1<-pivot_longer(data = as.data.frame(t(Rsquare.gene.test1)),cols=!Mu,names_to="method",values_to="value")

    ari.train1$se<-ari.train1.se$value/sqrt(50)
    ari.test1$se<-ari.test1.se$value/sqrt(50)
    cindex.train1$se<-cindex.train1.se$value/sqrt(50)
    cindex.test1$se<-cindex.test1.se$value/sqrt(50)
    jaccard1$se<-jaccard1.se$value/sqrt(50)
    F.outcome.train1$se<-F.outcome.train1.se$value/sqrt(50)
    F.outcome.test1$se<-F.outcome.test1.se$value/sqrt(50)
    Rsquare.gene.train1$se<-Rsquare.gene.train1.se$value/sqrt(50)
    Rsquare.gene.test1$se<-Rsquare.gene.test1.se$value/sqrt(50)


    Group<-c(rep("ARI.train",nrow(ari.train1)),rep("ARI.test",nrow(ari.train1)),rep("cindex.train",nrow(ari.train1)),
             rep("cindex.test",nrow(ari.train1)),rep("Jaccard",nrow(ari.train1)),
             rep("F.outcome.train",nrow(ari.train1)),rep("R2.gene.train",nrow(ari.train1)),
             rep("F.outcome.test",nrow(ari.train1)),rep("R2.gene.test",nrow(ari.train1)))
    data1<-rbind(ari.train1,ari.test1,cindex.train1,cindex.test1,jaccard1,
                 F.outcome.train1,Rsquare.gene.train1,F.outcome.test1,Rsquare.gene.test1)
    data1<-cbind(data1,Group)
    data1$Group<-factor(data1$Group,levels = c("ARI.train","cindex.train","Jaccard",
                                               "ARI.test","cindex.test",
                                               "F.outcome.train","R2.gene.train",
                                               "F.outcome.test","R2.gene.test"))
    data1$method=factor(data1$method,levels = c("ogClust_GM",  "ogClust_WJL", "SKM","PMBC"))
    data<-data1[which(data1$Group=="ARI.train"),]
    p<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=10, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+facet_wrap(~Group,scales = "free")+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"),
            strip.text.x = element_text(size = 20),legend.text = element_text(size=25),
            legend.title = element_text(size=25),
            legend.key.size = unit(1.5, 'cm'))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+
      scale_color_discrete(name = "method", labels = c(expression(ogClust[GM]), expression(ogClust[WJL]), expression(paste("sparse ",italic(K),"-means")),expression("PMBC")))+
      scale_shape_discrete(name = "method", labels = c(expression(ogClust[GM]), expression(ogClust[WJL]), expression(paste("sparse ",italic(K),"-means")),expression("PMBC")))+
      theme(legend.text.align = 0)
    p
    g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)}
    mylegend<-g_legend(p)

    p1<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab("ARI")
    p1

    data<-data1[which(data1$Group=="cindex.train"),]
    p2<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab("Adjusted C-index")
    p2

    data<-data1[which(data1$Group=="F.outcome.train"),]
    p3<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab("log-rank test statistics")
    p3

    data<-data1[which(data1$Group=="R2.gene.train"),]
    p4<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab(expression(paste(italic(R)^2, "(genes)",sep="")))
    p4


    data<-data1[which(data1$Group=="ARI.test"),]
    p5<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab("ARI")
    p5

    data<-data1[which(data1$Group=="cindex.test"),]
    p6<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab("Adjusted C-index")
    p6

    data<-data1[which(data1$Group=="F.outcome.test"),]
    p7<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab("log-rank test statistics")
    p7

    data<-data1[which(data1$Group=="R2.gene.test"),]
    p8<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab(expression(paste(italic(R)^2, "(genes)",sep="")))
    p8

    data<-data1[which(data1$Group=="Jaccard"),]
    p9<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab("Jaccard")
    p9

    p10<-as_ggplot(mylegend)+theme(plot.margin = unit(c(5,5,5,5), "cm"))
    p11<-ggplot() +theme_void()
    p12<-ggplot() +theme_void()
    # p1 <- annotate_figure(p1,top = text_grob("(A1) clustering accuracy",face = "bold", size = 20 ))
    # p2 <- annotate_figure(p2,top = text_grob("(A2) outcome prediction",face = "bold", size = 20 ))
    # p3 <- annotate_figure(p3,top = text_grob("(A3) outcome prediction",face = "bold", size = 20 ))
    # p4 <- annotate_figure(p4,top = text_grob("(A4) gene separation",face = "bold", size = 20 ))
    # p5 <- annotate_figure(p5,top = text_grob("(B1) clustering accuracy",face = "bold", size = 20 ))
    # p6 <- annotate_figure(p6,top = text_grob("(B2) outcome prediction",face = "bold", size = 20 ))
    # p7 <- annotate_figure(p7,top = text_grob("(B3) outcome prediction",face = "bold", size = 20 ))
    # p8 <- annotate_figure(p8,top = text_grob("(B4) gene separation",face = "bold", size = 20 ))
    # p9 <- annotate_figure(p9,top = text_grob("(A5) feature selection",face = "bold", size = 20 ))
    #
    # plot1<-ggarrange(p1, p3, p9,nrow=1)
    # plot2<-ggarrange(p2, p4, p10,nrow=1)
    # plot3<-ggarrange(p5, p7, p11,nrow=1)
    # plot4<-ggarrange(p6, p8, p12,nrow=1)
    # plot.add<-ggarrange(p11, p11, p11,nrow=1)
    #
    # plot_combined<-grid.arrange(plot1,plot2,plot.add,plot3,plot4,
    #                             layout_matrix = cbind(c(1,1,2,2,3,4,4,5,5)))
    #
    #
    # if(SC & tune) {
    #   file.name<-paste("SimII_Figure_c1=",c1_vector_all[ind_c1],".SC.unknown.png",sep="")
    # }
    # else if (SC & !tune) {
    #   file.name<-paste("SimII_Figure_c1=",c1_vector_all[ind_c1],".SC.known.png",sep="")
    # }
    # else if (!SC & tune) {
    #   file.name<-paste("SimII_Figure_c1=",c1_vector_all[ind_c1],".unknown.png",sep="")
    # }
    # else {
    #   file.name<-paste("SimII_Figure_c1=",c1_vector_all[ind_c1],".known.png",sep="")
    # }
    # ggsave(plot_combined,filename = file.name, width = 30, height = 40, units = "cm")

    p1 <- annotate_figure(p1,top = text_grob("(C1) Train: clustering accuracy",face = "bold", size = 20 ))
    p2 <- annotate_figure(p2,top = text_grob("(F1) Train: outcome prediction",face = "bold", size = 20 ))
    p3 <- annotate_figure(p3,top = text_grob("(E1) Train: outcome prediction",face = "bold", size = 20 ))
    p4 <- annotate_figure(p4,top = text_grob("(D1) Train: gene separation",face = "bold", size = 20 ))
    p5 <- annotate_figure(p5,top = text_grob("(C2) Test: clustering accuracy",face = "bold", size = 20 ))
    p6 <- annotate_figure(p6,top = text_grob("(F2) Test: outcome prediction",face = "bold", size = 20 ))
    p7 <- annotate_figure(p7,top = text_grob("(E2) Test: outcome prediction",face = "bold", size = 20 ))
    p8 <- annotate_figure(p8,top = text_grob("(D2) Test: gene separation",face = "bold", size = 20 ))
    p9 <- annotate_figure(p9,top = text_grob("(B) feature selection",face = "bold", size = 20 ))
    plot1<-ggarrange(p11, p11,p9,nrow=1)
    plot2<-ggarrange(p11, p11,p10,nrow=1)
    plot3<-ggarrange(p1, p4,p3, p2,nrow=1)
    plot4<-ggarrange(p5, p8,p7, p6,nrow=1)

    plot_combined<-grid.arrange(plot1,plot2,plot3,plot4,
                                layout_matrix = cbind(c(1,1,2,2,3,3,3,4,4,4)))
    ggsave(plot_combined,filename = figure.name, width = 50, height = 40, units = "cm")

  }
  print("WJL # features:")
  print(num.yujia.matrix.num)
  print("GM # features:")
  print(num.Peng.matrix.num)
  print("SKM # features:")
  print(num.skmeans.matrix.num)
  print("PMBC # features:")
  print(num.PMBC.matrix.num)
}

#generate Figrue II
summarySurv_gen(SC=FALSE,tune=TRUE,c1_vector_all=c(2),figure.name = "FigureII.png")
#generate Supplementary Figrue V
summarySurv_gen(SC=FALSE,tune=FALSE,c1_vector_all=c(2),figure.name="SupplementFigrueV.png")

