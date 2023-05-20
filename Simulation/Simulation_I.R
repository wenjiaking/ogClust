rm(list=ls()) 
library(sparcl) #to implement sparse K means
library(adaHuber)
library(parallel)
library(Brobdingnag)
library(glmnet)
library(mclust)
library(psych)
library(hms)


####Generate simulated data####
mu_vector<-c(0.9,1.2,1.5,1.8)
mu1_vector<-mu_vector+0.2
c1_vector_all<-c(3,2,1)
beta_vector<-2
sigma_y<-1
var_g<-2
n<-99
q<-50
q1<-50
q2<-50
q3<-1850
K<-3
num.sim<-50

for(ind_c1 in 1:length(c1_vector_all)) {
  c1<-c1_vector_all[ind_c1]
  for (ind_mu in 1:length(mu_vector)) {
    mu<-mu_vector[ind_mu]
    mu1<-mu1_vector[ind_mu]
    file.name<-paste("SimI_Data_c1=",c1,"_mu=",mu,".rds",sep="")
    data.temp=mclapply(1:num.sim,function(ind.data) {
      set.seed(ind.data)
      data<-Sim3(n=n,beta1=beta_vector,beta2=beta_vector,q=q,q1=q1,q2=q2,q3=q3,c1=c1,var_g=var_g,mu=mu,mu1=mu1,sigma_y=sigma_y)
      X<-data$x
      Y<-data$y
      G<-data$G
      true.label<-c(rep(1,33),rep(2,33),rep(3,33))
      true.feature<-c(rep(1,q),rep(0,2000-q))
      return(list(X=X,Y=Y,G=G,true.label=true.label,true.feature=true.feature))
    },mc.cores=num.sim)
    saveRDS(data.temp,file.name)
  }
}

num.sim<-50
mu_vector<-c(0.9,1.2,1.5,1.8)
mu1_vector<-mu_vector+0.2
c1_vector_all<-c(3,1)

res.gen=function(K=3,SC=TRUE,num.sim=50,mu_vector=c(0.9,1.2,1.5,1.8),c1_vector_all=c(3,1)) {
  
  for(ind_c1 in 1:length(c1_vector_all)) {
    c1<-c1_vector_all[ind_c1]
    for (ind_mu in 1:length(mu_vector)) {
      mu=mu_vector[ind_mu]
      file.name<-paste("SimI_Data_c1=",c1,"_mu=",mu,".rds",sep="")
      print(paste0("Begin ",file.name))
      data.temp=readRDS(file.name)
      RES=mclapply(1:num.sim,function(ind.data) {
        data=data.temp[[ind.data]]
        X<-data$X
        Y<-data$Y
        G<-data$G
        true.label<-data$true.label
        true.feature<-data$true.feature
        set.seed(ind.data)
        index.test<-sample(1:nrow(X),round(nrow(X)/3))
        index.train<-(1:nrow(X))[-index.test]
        X.train<-X[index.train,]
        G.train<-G[index.train,]
        G.train<-scale(G.train)
        Y.train<-Y[index.train]
        X.test<-X[index.test,]
        G.test<-G[index.test,]
        G.test<-scale(G.test)
        Y.test<-Y[index.test]
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
            mod<-ogclust_con_Yujia(x=X.train,G=t(G.train),y=Y.train,c_center=center,lambda=lambda_vector[j],v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
            num.vector[j]<-sum(apply(mod$result_list$mu,1,function(x){length(unique(x))})!=1)
          }
          index<-which.min(abs(num.vector-50))
          lambda<-lambda_vector[index]
          
          mod.all<-ogclust_con_Yujia(x=X.train,G=t(G.train),y=Y.train,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
          #cluster1<-predict.ogclust.test(mod.all,K=K,D.test=as.matrix(G.train),D.train=t(G.train),X1=X.train,p=ncol(G.train))
          cluster1<-predict_test(mod.all,K=K,D.test=as.matrix(G.train),X1=as.matrix(X.train),p=ncol(G.train))
          cluster1<-cluster1$clus
          if (length(table(cluster1))!=K) {
            temp.R2[ind_w]=NA
            next
          }
          centers<-cbind(apply(G.train[which(cluster1==1),],2,mean),
                         apply(G.train[which(cluster1==2),],2,mean),
                         apply(G.train[which(cluster1==3),],2,mean))
          
          #print(table(cluster1))
          
          cluster<-rep(NA,nrow(G.train))
          for(ind_cv in 1:length(unique(index.sample))){
            ind_sample<-which(index.sample==ind_cv)
            G.test1<-G.train[ind_sample,]
            X.test1<-X.train[ind_sample,]
            Y.test1<-Y.train[ind_sample]
            G.train1<-G.train[(1:nrow(X.train))[-ind_sample],]
            X.train1<-X.train[(1:nrow(X.train))[-ind_sample],]
            Y.train1<-Y.train[(1:nrow(X.train))[-ind_sample]]
            
            label.test1<-label.train[ind_sample]
            label.train1<-label.train[(1:nrow(X.train))[-ind_sample]]
            num.vector<-rep(NA,length(lambda_vector))
            for(j in 1:length(lambda_vector)){
              mod<-ogclust_con_Yujia(x=X.train1,G=t(G.train1),y=Y.train1,c_center=center,lambda=lambda_vector[j],v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w_vector[ind_w],w_G=1-w_vector[ind_w],z_int=NULL)
              num.vector[j]<-sum(apply(mod$result_list$mu,1,function(x){length(unique(x))})!=1)
            }
            
            index<-which.min(abs(num.vector-50))
            lambda<-lambda_vector[index]
            #temp.num[ind_w,ind_cv]<-num.vector[index]
            mod<-ogclust_con_Yujia(x=X.train1,G=t(G.train1),y=Y.train1,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
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
          res1<-Rsquare(cluster = cluster,Y = Y.train,X = X.train,G = t(G.train))
          temp.R2[ind_w]<-res1$outcome
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
          mod<-ogclust_con_Yujia(x=X.train,G=t(G.train),y=Y.train,c_center=center,lambda=lambda_vector1[j],v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
          select.feature<-as.numeric(apply(mod$result_list$mu,1,function(x){length(unique(x))})!=1)
          num.vector[j]<-sum(select.feature==1)
          
          cluster1<-predict_test(mod,K=K,D.test=as.matrix(G.train),X1=as.matrix(X.train),p=ncol(G.train))
          cluster1<-cluster1$clus
          
          if(length(table(cluster1))!=K){
            R2.outcome.tune[j]=NA
            R2.gene.tune[j]=NA
            next
          }
          
          centers<-cbind(apply(G.train[which(cluster1==1),],2,mean),
                         apply(G.train[which(cluster1==2),],2,mean),
                         apply(G.train[which(cluster1==3),],2,mean))
          
          
          cluster<-rep(NA,nrow(G.train))
          for(ind_cv in 1:length(unique(index.sample))){
            ind_sample<-which(index.sample==ind_cv)
            G.test11<-G.train[ind_sample,]
            X.test1<-X.train[ind_sample,]
            Y.test1<-Y.train[ind_sample]
            G.train11<-G.train[(1:nrow(G.train))[-ind_sample],]
            X.train1<-X.train[(1:nrow(G.train))[-ind_sample],]
            Y.train1<-Y.train[(1:nrow(G.train))[-ind_sample]]
            
            mod<-ogclust_con_Yujia(x=X.train1,G=t(G.train11),y=Y.train1,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=3,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
            
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
          res1<-Rsquare(cluster = cluster,Y = Y.train,X =X.train,G = t(G.train))
          R2.outcome.tune[j]<-res1$outcome
          R2.gene.tune[j]<-mean(res1$gene[which(select.feature==1)])
        }
        
        #---------------------------------------------
        #use criteria sqrt(R2.outcome.tune*R2.gene.tune) to tune lambda
        #---------------------------------------------
        index<-which.max(sqrt(R2.outcome.tune*R2.gene.tune))
        lambda=lambda_vector1[index]
        print(paste0("By sqrt(R2.outcome*R2.gene) criteria, lambda=",lambda," that selects  ",num.vector[index]," features."))
        
        mod.tune<-ogclust_con_Yujia(x=X.train,G=t(G.train),y=Y.train,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
        
        res.test.tune<-predict_test(mod.tune,K=K,D.test=G.test,X1=as.matrix(X.test),p=ncol(G.test))
        res.train.tune<-predict.ogclust.test.select(mod.tune,K=K,D.test=as.matrix(G.train),D.train=t(G.train),X1=as.matrix(X.train),p=ncol(G.train),s_G = s_G,O.test = Y.train,w = w)
        
        cluster.test<-res.test.tune$clus
        cluster.train = res.train.tune$clus
        
        ari.test.WJL.tune<-adjustedRandIndex(cluster.test,label.test)
        ari.train.WJL.tune<-adjustedRandIndex(cluster.train,label.train)
        Rmse.train.WJL.tune<-sqrt(sum((res.train.tune$Y-Y.train)^2)/length(Y.train))
        Rmse.train.WJL.tune.hard<-sqrt(sum((res.train.tune$Y.hard-Y.train)^2)/length(Y.train))
        Rmse.test.WJL.tune<-sqrt(sum((res.test.tune$Y-Y.test)^2)/length(Y.test))
        Rmse.test.WJL.tune.hard<-sqrt(sum((res.test.tune$Y.hard-Y.test)^2)/length(Y.test))
        select.feature<-as.numeric(apply(mod.tune$result_list$mu,1,function(x){length(unique(x))})>1)
        num.WJL.tune=sum(select.feature)
        jaccard.index.WJL.tune<-Jaccard.index(true.feature, select.feature)
        
        mod.train<-Rsquare(cluster =cluster.train,Y = Y.train,X = X.train,G = t(G.train))
        Rsquare.gene.WJL.tune<-mean(mod.train$gene[which(select.feature==1)])
        Rsquare.outcome.WJL.tune<-mod.train$outcome
        mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test))
        Rsquare.gene.test.WJL.tune<-mean(mod.test$gene[which(select.feature==1)])
        Rsquare.outcome.test.WJL.tune<-mod.test$outcome
        
        #---------------------------------------------
        #use num to select lambda
        #---------------------------------------------
        index<-max(which(num.vector>=50))
        lambda.num<-lambda_vector1[index]
        print(paste0("Tuning lambda=",lambda.num," that selects  ",num.vector[index]," features."))
        
        start_time=Sys.time()
        mod.num<-ogclust_con_Yujia(x=X.train,G=t(G.train),y=Y.train,c_center=center,lambda=lambda.num,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
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
        Rmse.train.WJL.num<-sqrt(sum((res.train$Y-Y.train)^2)/length(Y.train))
        Rmse.train.WJL.num.hard<-sqrt(sum((res.train$Y.hard-Y.train)^2)/length(Y.train))
        Rmse.test.WJL.num<-sqrt(sum((res.test$Y-Y.test)^2)/length(Y.test))
        Rmse.test.WJL.num.hard<-sqrt(sum((res.test$Y.hard-Y.test)^2)/length(Y.test))
        
        select.feature<-as.numeric(apply(mod.num$result_list$mu,1,function(x){length(unique(x))})>1)
        jaccard.index.WJL.num<-Jaccard.index(true.feature, select.feature)
        
        mod.train<-Rsquare(cluster =cluster.train,Y = Y.train,X = X.train,G = t(G.train))
        Rsquare.gene.WJL<-mean(mod.train$gene[which(select.feature==1)])
        Rsquare.outcome.WJL<-mod.train$outcome
        mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test))
        Rsquare.gene.test.WJL<-mean(mod.test$gene[which(select.feature==1)])
        Rsquare.outcome.test.WJL<-mod.test$outcome
        num.WJL.num<-sum(select.feature==1)
        
        #w.est<-w
        G.train11<-G.train[,which(select.feature==1)]
        G.test11<-G.test[,which(select.feature==1)]
        mod.G<-ogclust_con_Yujia(x=X.train,G=t(G.train11),y=Y.train,c_center=center[which(select.feature==1),],lambda=0,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=0,w_G=1,z_int=NULL)
        mod.O<-ogclust_con_Yujia(x=X.train,G=t(G.train11),y=Y.train,c_center=center[which(select.feature==1),],lambda=0,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=1,w_G=0,z_int=NULL)
        cluster.G<-apply(mod.G$result_list$z,1,which.max)
        cluster.O<-apply(mod.O$result_list$z,1,which.max)
        ari.consistent.WJL<-adjustedRandIndex(cluster.G,cluster.O)
        # res.G<-predict.ogclust.test(mod.G,K=K,D.test=G.test11,D.train=t(G.train11),X1=X.test,p=ncol(G.test11))
        # res.O<-predict.ogclust.test(mod.O,K=K,D.test=G.test11,D.train=t(G.train11),X1=X.test,p=ncol(G.test11))
        res.G<-predict_test(mod.G,K=K,D.test=G.test11,X1=X.test,p=ncol(G.test11))
        res.O<-predict_test(mod.O,K=K,D.test=G.test11,X1=X.test,p=ncol(G.test11))
        ari.consistent.test.WJL<-adjustedRandIndex(res.G$clus,res.O$clus)
        #res.test.num<-res.test
        #-------------------------------------------------
        WJL.res<-list(mod.tune=mod.tune,mod.num=mod.num,
                      pred.train.tune=res.train.tune,pred.train.num=res.train,
                      pred.test.tune=res.test.tune,pred.test.num=res.test,
                      ari.train.tune=ari.train.WJL.tune,ari.train.num=ari.train.WJL.num,
                      ari.test.tune=ari.test.WJL.tune,ari.test.num=ari.test.WJL.num,
                      Rmse.train.tune=Rmse.train.WJL.tune,Rmse.train.num=Rmse.train.WJL.num,
                      Rmse.test.tune=Rmse.test.WJL.tune,Rmse.test.num=Rmse.test.WJL.num,
                      Rmse.test.tune.hard=Rmse.test.WJL.tune.hard,Rmse.test.num.hard=Rmse.test.WJL.num.hard,
                      Rmse.train.tune.hard=Rmse.train.WJL.tune.hard,Rmse.train.num.hard=Rmse.train.WJL.num.hard,
                      jaccard.index.tune=jaccard.index.WJL.tune,jaccard.index.num=jaccard.index.WJL.num,
                      num.tune=num.WJL.tune,num.num=num.WJL.num,
                      R2.gene.train.tune=Rsquare.gene.WJL.tune,
                      R2.outcome.train.tune=Rsquare.outcome.WJL.tune,
                      R2.gene.test.tune=Rsquare.gene.test.WJL.tune,
                      R2.outcome.test.tune=Rsquare.outcome.test.WJL.tune,
                      R2.gene.train=Rsquare.gene.WJL,
                      R2.outcome.train=Rsquare.outcome.WJL,
                      R2.gene.test=Rsquare.gene.test.WJL,
                      R2.outcome.test=Rsquare.outcome.test.WJL,
                      lambda.tune=lambda,lambda.num=lambda.num,
                      train_time=time.WJL,
                      w0.est=w0,
                      ari.consistent.train=ari.consistent.WJL,
                      ari.consistent.test=ari.consistent.test.WJL)
        print("Finish WJL method!")
        #------------------------------
        #Peng's method
        #------------------------------
        n1=nrow(G.train) # number of samples
        NG=ncol(G.train) # number of genes
        np=ncol(X.train) # number of covariates
        #lambda_vector_Peng<-c(0.005,0.01,0.05,seq(0.1,1,0.1))
        lambda_vector_Peng<-c(seq(0.005,0.01,0.001),seq(0.02,0.5,0.02))
        
        mod.kmeans<-kmeans(G.train,centers = K,nstart = 50)
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
          fit.res<-fit.ogClust(n=n1, K=K, np=np, NG=NG, lambda=lambda_vector_Peng[ind_lambda],
                               alpha=0.5, G=as.data.frame(G.train), Y=Y.train, X=X.train, theta_int=theta_int)
          
          beta=fit.res$par$gamma 
          beta<-as.data.frame(beta)
          select.feature<-as.numeric(apply(beta,1,function(x){sum(x==0)})!=K)[-1] #genes with non-zero weight of at least one clsuter
          num[ind_lambda]<-sum(select.feature==1)
          print(sum(select.feature==1))
          print(num)
          cluster1<-fit.res$grp_assign 
          
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
            G.train1<-G.train[(1:nrow(G.train))[-ind_sample],]
            X.train1<-X.train[(1:nrow(G.train))[-ind_sample],]
            Y.train1<-Y.train[(1:nrow(G.train))[-ind_sample]]
            n1sub=nrow(G.train1) 
            NGsub=ncol(G.train1) 
            np=ncol(X.train1)
            mod<-fit.ogClust(n=n1sub, K=K, np=np, NG=NGsub, lambda=lambda_vector_Peng[ind_lambda],
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
            res1<-Rsquare(cluster = cluster,Y = Y.train,X = X.train,G = t(G.train))
            R2.outcome.tune[ind_lambda]<-res1$outcome
            R2.gene.tune[ind_lambda]<-mean(res1$gene[which(select.feature==1)])
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
        
        fit.res.tune<-fit.ogClust(n=n1, K=K, np=np, NG=NG, lambda=lambda,
                                  alpha=0.5, G=as.data.frame(G.train), Y=Y.train, X=X.train, theta_int=theta_int)
        
        
        beta=fit.res.tune$par$gamma 
        beta<-as.data.frame(beta)
        select.feature<-as.numeric(apply(beta,1,function(x){sum(x==0)})!=K)[-1] #genes with non-zero weight of at least one clsuter
        
        jaccard.index.Peng.tune<-Jaccard.index(true.feature, select.feature)
        num.Peng.tune<-sum(select.feature==1)
        ari.train.Peng.tune<-adjustedRandIndex(fit.res.tune$grp_assign,label.train)
        mod.train<-Rsquare(cluster =fit.res.tune$grp_assign,Y = Y.train,X = X.train,G = t(G.train))
        Rsquare.gene.train.tune<-mean(mod.train$gene[which(select.feature==1)])
        Rsquare.outcome.train.tune<-mod.train$outcome
        
        beta_est=fit.res.tune$par$beta
        gamma_est_matrix=fit.res.tune$par$gamma
        beta0_est=fit.res.tune$par$beta0
        sigma2_est=fit.res.tune$par$sigma2
        
        intercept<-rep(NA,nrow(X.train))
        for(i in 1:K){
          intercept[which(fit.res.tune$grp_assign==i)]<-beta0_est[i]
        }
        
        pred.train<-X.train %*% beta_est+intercept
        Rmse.train.Peng.tune.hard<-sqrt(sum((pred.train-Y.train)^2)/length(Y.train))
        Rmse.train.Peng.tune<-sqrt(sum((fit.res.tune$Y_prd-Y.train)^2)/length(Y.train))
        
        
        G.test.add = cbind(1, G.test)
        pai_est.tune = sapply(1:K, function(k) exp(G.test.add %*% gamma_est_matrix[,  k, drop = F])/rowSums(exp(G.test.add %*% gamma_est_matrix)))
        cluster.test<-apply(pai_est.tune,1,which.max)
        if(length(unique(cluster.test))>1){
          mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test[,-1] ))
          Rsquare.gene.test.tune<-mean(mod.test$gene[which(select.feature==1)])
          Rsquare.outcome.test.tune<-mod.test$outcome
        }else{
          Rsquare.gene.test.tune<-NA
          Rsquare.outcome.test.tune<-NA
        }
        
        intercept<-rep(NA,nrow(X.test))
        for(i in 1:K){
          intercept[which(cluster.test==i)]<-beta0_est[i]
        }
        
        pred.test<-X.test %*% beta_est+intercept
        Rmse.test.Peng.tune.hard<-sqrt(sum((pred.test-Y.test)^2)/length(Y.test))
        y.pred.test = apply(sapply(1:K, function(x) pai_est.tune[, x] * (beta0_est[x] + X.test %*% beta_est)), 1, sum)
        Rmse.test.Peng.tune<-sqrt(sum((y.pred.test-Y.test)^2)/length(Y.test))
        ari.test.Peng.tune<-adjustedRandIndex(cluster.test,label.test)
        
        #---------------------------------------------
        #use num to select lambda
        #---------------------------------------------
        index<-which.min(abs(num-50))
        # index<-min(which(num>50))
        lambda.num=lambda_vector_Peng[index]
        print(paste0("Tuning lambda=",lambda.num," that selects  ",num[index]," features."))
        start_time=Sys.time()
        fit.res.num<-fit.ogClust(n=n1, K=K, np=np, NG=NG, lambda=lambda.num,
                                 alpha=0.5, G=as.data.frame(G.train), Y=Y.train, X=X.train, theta_int=theta_int)
        end_time=Sys.time()
        time.GM=as_hms(end_time-start_time)
        
        beta=fit.res.num$par$gamma 
        beta<-as.data.frame(beta)
        select.feature<-as.numeric(apply(beta,1,function(x){sum(x==0)})!=K)[-1]
        
        jaccard.index.Peng.num<-Jaccard.index(true.feature, select.feature)
        num.Peng.num<-sum(select.feature==1)
        
        mod.train<-Rsquare(cluster =fit.res.num$grp_assign,Y = Y.train,X = X.train,G = t(G.train))
        Rsquare.gene.train<-mean(mod.train$gene[which(select.feature==1)])
        Rsquare.outcome.train<-mod.train$outcome
        
        ari.train.Peng.num<-adjustedRandIndex(fit.res.num$grp_assign,label.train)
        beta_est=fit.res.num$par$beta
        gamma_est_matrix=fit.res.num$par$gamma
        beta0_est=fit.res.num$par$beta0
        sigma2_est=fit.res.num$par$sigma2
        intercept<-rep(NA,nrow(X.train))
        for(i in 1:K){
          intercept[which(fit.res.num$grp_assign==i)]<-beta0_est[i]
        }
        
        pred.train<-X.train %*% beta_est+intercept
        Rmse.train.Peng.num.hard<-sqrt(sum((pred.train-Y.train)^2)/length(Y.train))
        Rmse.train.Peng.num<-sqrt(sum((fit.res.num$Y_prd-Y.train)^2)/length(Y.train))
        
        G.test.add = cbind(1, G.test)
        pai_est.num = sapply(1:K, function(k) exp(G.test.add %*% gamma_est_matrix[,  k, drop = F])/rowSums(exp(G.test.add %*% gamma_est_matrix)))
        
        cluster.test<-apply(pai_est.num,1,which.max)
        if(length(unique(cluster.test))>1){
          mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test[,-1] ))
          Rsquare.gene.test<-mean(mod.test$gene[which(select.feature==1)])
          Rsquare.outcome.test<-mod.test$outcome
        }else{
          Rsquare.gene.test<-NA
          Rsquare.outcome.test<-NA
        }
        
        
        intercept<-rep(NA,nrow(X.test))
        for(i in 1:K){
          intercept[which(cluster.test==i)]<-beta0_est[i]
        }
        
        pred.test<-X.test %*% beta_est+intercept
        Rmse.test.Peng.num.hard<-sqrt(sum((pred.test-Y.test)^2)/length(Y.test))
        
        y.pred.test = apply(sapply(1:K, function(x) pai_est.num[, x] * (beta0_est[x] + X.test %*% beta_est)), 1, sum)
        Rmse.test.Peng.num<-sqrt(sum((y.pred.test-Y.test)^2)/length(Y.test))
        ari.test.Peng.num<-adjustedRandIndex(cluster.test,label.test)
        #-------------------------------------------------
        GM.res<-list(mod.tune=fit.res.tune,mod.num=fit.res.num,
                     ari.train.tune=ari.train.Peng.tune,ari.train.num=ari.train.Peng.num,
                     ari.test.tune=ari.test.Peng.tune,ari.test.num=ari.test.Peng.num,
                     Rmse.train.tune=Rmse.train.Peng.tune,Rmse.train.num=Rmse.train.Peng.num,
                     Rmse.train.tune.hard=Rmse.train.Peng.tune.hard,Rmse.train.num.hard=Rmse.train.Peng.num.hard,
                     Rmse.test.tune=Rmse.test.Peng.tune,Rmse.test.num=Rmse.test.Peng.num,
                     Rmse.test.tune.hard=Rmse.test.Peng.tune.hard,Rmse.test.num.hard=Rmse.test.Peng.num.hard,
                     jaccard.index.tune=jaccard.index.Peng.tune,jaccard.index.num=jaccard.index.Peng.num,
                     num.tune=num.Peng.tune,num.num=num.Peng.num,
                     R2.gene.test.tune=Rsquare.gene.test.tune,
                     R2.gene.train.tune=Rsquare.gene.train.tune,
                     R2.outcome.test.tune=Rsquare.outcome.test.tune,
                     R2.outcome.train.tune=Rsquare.outcome.train.tune,
                     R2.gene.test=Rsquare.gene.test,
                     R2.gene.train=Rsquare.gene.train,
                     R2.outcome.test=Rsquare.outcome.test,
                     R2.outcome.train=Rsquare.outcome.train,
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
        mod.train<-Rsquare(cluster =mod.tune[[1]]$Cs,Y = Y.train,X = X.train,G = t(G.train ))
        Rsquare.gene.train.tune<-mean(mod.train$gene[which(select.feature==1)])
        Rsquare.outcome.train.tune<-mod.train$outcome
        
        ari.train.kmeans.tune<-adjustedRandIndex(mod.tune[[1]]$Cs,label.train)
        centers<-cbind(apply(G.train[which(mod.tune[[1]]$Cs==1),],2,mean),
                       apply(G.train[which(mod.tune[[1]]$Cs==2),],2,mean),
                       apply(G.train[which(mod.tune[[1]]$Cs==3),],2,mean))
        
        distance<-apply(G.test,1,function(x){
          distance1<-apply(centers,2,function(y){
            sum(weight*(x-y)^2)
          })
        })
        cluster.test<-apply(distance,2,which.min)
        mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test))
        Rsquare.gene.test.tune<-mean(mod.test$gene[which(select.feature==1)])
        Rsquare.outcome.test.tune<-mod.test$outcome
        
        ari.test.kmeans.tune<-adjustedRandIndex(cluster.test,label.test)
        pred.train.kmeans<-rep(NA,length(Y.train))
        pred.test.kmeans<-rep(NA,length(Y.test))
        for(i in 1:K){
          data.temp<-data.frame(Y=Y.train[which(mod.tune[[1]]$Cs==i)],X.train[which(mod.tune[[1]]$Cs==i),])
          mod.lm<-lm(Y~.,data=data.temp)
          pred.train.kmeans[which(mod.tune[[1]]$Cs==i)]<-mod.lm$fitted.values
          pred.test.kmeans[which(cluster.test==i)]<-predict(mod.lm,newdata = as.data.frame(X.test[which(cluster.test==i),]))
        }
        Rmse.train.kmeans.tune<-sqrt(sum((pred.train.kmeans-Y.train)^2)/length(Y.train))
        Rmse.test.kmeans.tune<-sqrt(sum((pred.test.kmeans-Y.test)^2)/length(Y.test))
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
        Rsquare.gene.train<-mean(mod.train$gene[which(select.feature==1)])
        Rsquare.outcome.train<-mod.train$outcome
        #mod.kmeans<-kmeans(G.train1,center=K,nstart=20)
        ari.train.kmeans.num<-adjustedRandIndex(mod.num$Cs,label.train)
        centers<-cbind(apply(G.train[which(mod.num$Cs==1),],2,mean),
                       apply(G.train[which(mod.num$Cs==2),],2,mean),
                       apply(G.train[which(mod.num$Cs==3),],2,mean))
        
        distance<-apply(G.test,1,function(x){
          distance1<-apply(centers,2,function(y){
            sum(weight*(x-y)^2)
          })
        })
        cluster.test.num<-apply(distance,2,which.min)
        mod.test<-Rsquare(cluster =cluster.test.num,Y = Y.test,X = X.test,G = t(G.test))
        Rsquare.gene.test<-mean(mod.test$gene[which(select.feature==1)])
        Rsquare.outcome.test<-mod.test$outcome
        
        ari.test.kmeans.num<-adjustedRandIndex(cluster.test.num,label.test)
        pred.train.kmeans<-rep(NA,length(Y.train))
        pred.test.kmeans<-rep(NA,length(Y.test))
        for(i in 1:K){
          data.temp<-data.frame(Y=Y.train[which(mod.num$Cs==i)],X.train[which(mod.num$Cs==i),])
          mod.lm<-lm(Y~.,data=data.temp)
          pred.train.kmeans[which(mod.num$Cs==i)]<-mod.lm$fitted.values
          pred.test.kmeans[which(cluster.test.num==i)]<-predict(mod.lm,newdata = as.data.frame(X.test[which(cluster.test.num==i),]))
        }
        Rmse.train.kmeans.num<-sqrt(sum((pred.train.kmeans-Y.train)^2)/length(Y.train))
        Rmse.test.kmeans.num<-sqrt(sum((pred.test.kmeans-Y.test)^2)/length(Y.test))
        jaccard.index.skmeans.num<-Jaccard.index(TrueFeature =true.feature,SelectedFeature1 =select.feature )
        
        #-------------------------------------------------
        SKM.res<-list(mod.tune=mod.tune,mod.num=mod.num,
                      ari.train.tune=ari.train.kmeans.tune,ari.train.num=ari.train.kmeans.num,
                      ari.test.tune=ari.test.kmeans.tune,ari.test.num=ari.test.kmeans.num,
                      Rmse.train.tune=Rmse.train.kmeans.tune,Rmse.train.num=Rmse.train.kmeans.num,
                      Rmse.test.tune=Rmse.test.kmeans.tune,Rmse.test.num=Rmse.test.kmeans.num,
                      num.tune=num.skmeans.tune,num.num=num.skmeans.num,
                      jaccard.index.tune=jaccard.index.skmeans.tune,jaccard.index.num=jaccard.index.skmeans.num,
                      R2.gene.test.tune=Rsquare.gene.test.tune,
                      R2.gene.train.tune=Rsquare.gene.train.tune,
                      R2.outcome.test.tune=Rsquare.outcome.test.tune,
                      R2.outcome.train.tune=Rsquare.outcome.train.tune,
                      R2.gene.test=Rsquare.gene.test,
                      R2.gene.train=Rsquare.gene.train,
                      R2.outcome.test=Rsquare.outcome.test,
                      R2.outcome.train=Rsquare.outcome.train,
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
        
        pi_est_train<-mod.tune$optimal_result$pi
        mu_est_train<-mod.tune$optimal_result$mu
        sigma_est_train<-mod.tune$optimal_result$sigma
        grp_assign_train<-apply(mod.tune$optimal_result$z,1,which.max)
        cluster.train=grp_assign_train
        ari.train.PMBC.tune<-adjustedRandIndex(cluster.train,label.train)
        mod.train<-Rsquare(cluster =cluster.train,Y = Y.train,X = X.train,G = t(G.train ))
        Rsquare.gene.train.tune<-mean(mod.train$gene[mod.tune$selected.index])
        Rsquare.outcome.train.tune<-mod.train$outcome
        
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
        
        mod.test<-Rsquare(cluster =cluster.test,Y = Y.test,X = X.test,G = t(G.test))
        Rsquare.gene.test.tune<-mean(mod.test$gene[mod.tune$selected.index])
        Rsquare.outcome.test.tune<-mod.test$outcome
        
        pred.train.PMBC<-rep(NA,length(Y.train))
        pred.test.PMBC<-rep(NA,length(Y.test))
        for(i in 1:K){
          data.temp<-data.frame(Y=Y.train[which(cluster.train==i)],X.train[which(cluster.train==i),])
          mod.lm<-lm(Y~.,data=data.temp)
          pred.train.PMBC[which(cluster.train==i)]<-mod.lm$fitted.values
          pred.test.PMBC[which(cluster.test==i)]<-predict(mod.lm,newdata = as.data.frame(X.test[which(cluster.test==i),]))
        }
        Rmse.train.PMBC.tune<-sqrt(sum((pred.train.PMBC-Y.train)^2)/length(Y.train))
        Rmse.test.PMBC.tune<-sqrt(sum((pred.test.PMBC-Y.test)^2)/length(Y.test))
        select.feature=rep(0,ncol(G.train))
        select.feature[mod.tune$selected.index]=1
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
        Rsquare.gene.train<-mean(mod.train$gene[mod.num$selected.index])
        Rsquare.outcome.train<-mod.train$outcome
        
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
        Rsquare.gene.test<-mean(mod.test$gene[mod.num$selected.index])
        Rsquare.outcome.test<-mod.test$outcome
        
        pred.train.PMBC<-rep(NA,length(Y.train))
        pred.test.PMBC<-rep(NA,length(Y.test))
        for(i in 1:K){
          data.temp<-data.frame(Y=Y.train[which(cluster.train==i)],X.train[which(cluster.train==i),])
          mod.lm<-lm(Y~.,data=data.temp)
          pred.train.PMBC[which(cluster.train==i)]<-mod.lm$fitted.values
          pred.test.PMBC[which(cluster.test==i)]<-predict(mod.lm,newdata = as.data.frame(X.test[which(cluster.test==i),]))
        }
        Rmse.train.PMBC.num<-sqrt(sum((pred.train.PMBC-Y.train)^2)/length(Y.train))
        Rmse.test.PMBC.num<-sqrt(sum((pred.test.PMBC-Y.test)^2)/length(Y.test))
        select.feature=rep(0,ncol(G.train))
        select.feature[mod.num$selected.index]=1
        jaccard.index.PMBC.num<-Jaccard.index(true.feature,select.feature)
        
        #-------------------------------------------------
        PMBC.res<-list(mod.tune=mod.tune,mod.num=mod.num,
                       ari.train.tune=ari.train.PMBC.tune,ari.train.num=ari.train.PMBC.num,
                       ari.test.tune=ari.test.PMBC.tune,ari.test.num=ari.test.PMBC.num,
                       Rmse.train.tune=Rmse.train.PMBC.tune,Rmse.train.num=Rmse.train.PMBC.num,
                       Rmse.test.tune=Rmse.test.PMBC.tune,Rmse.test.num=Rmse.test.PMBC.num,
                       num.tune=num.PMBC.tune,num.num=num.PMBC.num,
                       jaccard.index.tune=jaccard.index.PMBC.tune,jaccard.index.num=jaccard.index.PMBC.num,
                       R2.gene.test.tune=Rsquare.gene.test.tune,
                       R2.gene.train.tune=Rsquare.gene.train.tune,
                       R2.outcome.test.tune=Rsquare.outcome.test.tune,
                       R2.outcome.train.tune=Rsquare.outcome.train.tune,
                       R2.gene.test=Rsquare.gene.test,
                       R2.gene.train=Rsquare.gene.train,
                       R2.outcome.test=Rsquare.outcome.test,
                       R2.outcome.train=Rsquare.outcome.train,
                       lambda.tune=lambda.tune,lambda.num=lambda.num,
                       train_time=time.PMBC)
        return(list(WJL.res=WJL.res,GM.res=GM.res,SKM.res=SKM.res,PMBC.res=PMBC.res,index.test=index.test))
      },mc.cores=150)
      if (SC) {
        file.res<-paste("SimI_RES_c1=",c1,"_mu=",mu,".SC.rds",sep="")
      }
      else {
        file.res<-paste("SimI_RES_c1=",c1,"_mu=",mu,".rds",sep="")
      }
      
      saveRDS(RES,file.res)
    }
  }
}

res.SC=res.gen(K=3,SC=TRUE,num.sim=50,mu_vector=c(0.9,1.2,1.5,1.8),c1_vector_all=c(3,1)) 
res=res.gen(K=3,SC=FALSE,num.sim=50,mu_vector=c(0.9,1.2,1.5,1.8),c1_vector_all=c(3,1))

library(gridExtra)
library(ggpubr)
library(tidyr)
library(ggplot2)
summary_gen=function(SC=TRUE,tune=FALSE,mu_vector=c(0.9,1.2,1.5,1.8),c1_vector_all=c(3,1),figure.name="SupplementFigrueI.png"){
  ari.train.yujia.matrix<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.yujia.matrix<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.yujia.matrix<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.yujia.matrix<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.yujia.matrix<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  #num.yujia.matrix<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  
  ari.train.yujia.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.yujia.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.yujia.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.yujia.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.yujia.matrix.num.hard<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.yujia.matrix.num.hard<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.yujia.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  num.yujia.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.train<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.test<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  
  ari.train.yujia.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.yujia.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.yujia.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.yujia.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.yujia.matrix.num.hard.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.yujia.matrix.num.hard.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.yujia.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.train.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.test.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  
  w.yujia.matrix<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  w.yujia.matrix.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  
  ari.train.Peng.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.Peng.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.Peng.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.Peng.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.Peng.matrix.num.hard<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.Peng.matrix.num.hard<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.Peng.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  num.Peng.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.train.Peng<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.Peng<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.test.Peng<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.Peng<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  
  ari.train.Peng.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.Peng.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.Peng.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.Peng.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.Peng.matrix.num.hard.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.Peng.matrix.num.hard.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.Peng.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.train.Peng.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.Peng.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.test.Peng.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.Peng.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  
  
  
  ari.train.skmeans.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.skmeans.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.skmeans.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.skmeans.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.skmeans.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  num.skmeans.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.train.skmeans<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.skmeans<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.test.skmeans<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.skmeans<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  
  ari.train.skmeans.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.skmeans.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.skmeans.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.skmeans.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.skmeans.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.train.skmeans.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.skmeans.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.test.skmeans.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.skmeans.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  
  ari.train.PMBC.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.PMBC.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.PMBC.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.PMBC.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.PMBC.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  num.PMBC.matrix.num<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.train.PMBC<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.PMBC<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.test.PMBC<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.PMBC<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  
  ari.train.PMBC.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  ari.test.PMBC.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.train.PMBC.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rmse.test.PMBC.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Jaccard.PMBC.matrix.num.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.train.PMBC.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.train.PMBC.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.outcome.matrix.test.PMBC.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  Rsquare.gene.matrix.test.PMBC.se<-matrix(NA,nrow=length(c1_vector_all),ncol=length(mu_vector),dimnames = list(paste("c1=",c1_vector_all,sep=""),1:length(mu_vector)))
  
  
  for(ind_c1 in 1:length(c1_vector_all)) {
    c1<-c1_vector_all[ind_c1]
    for (ind_mu in 1:length(mu_vector)) {
      mu=mu_vector[ind_mu]
      if (SC) {
        file.name= paste("SimI_RES_c1=",c1,"_mu=",mu,".SC.rds",sep="")
      }
      else {
        file.name= paste("SimI_RES_c1=",c1,"_mu=",mu,".rds",sep="")
      }
      RES=readRDS(file.name)
      RES=RES[sapply(RES,length)==5]
      if (tune) {
        ari.train.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$ari.train.tune})),na.rm=T)
        ari.test.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$ari.test.tune})),na.rm=T)
        Rmse.train.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$Rmse.train.tune})),na.rm=T)
        Rmse.test.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$Rmse.test.tune})),na.rm=T)
        Rmse.train.yujia.matrix.num.hard[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$Rmse.train.tune.hard})),na.rm=T)
        Rmse.test.yujia.matrix.num.hard[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$Rmse.test.tune.hard})),na.rm=T)
        Jaccard.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$jaccard.index.tune})),na.rm=T)
        num.yujia.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$WJL.res$num.tune})))
        Rsquare.outcome.matrix.train[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$R2.outcome.train.tune})),na.rm=T)
        Rsquare.outcome.matrix.test[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$R2.outcome.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.test.tune})),na.rm=T)
        
        ari.train.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$ari.train.tune})),na.rm=T)
        ari.test.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$ari.test.tune})),na.rm=T)
        Rmse.train.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$Rmse.train.tune})),na.rm=T)
        Rmse.test.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$Rmse.test.tune})),na.rm=T)
        Rmse.train.yujia.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$Rmse.train.tune.hard})),na.rm=T)
        Rmse.test.yujia.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$Rmse.test.tune.hard})),na.rm=T)
        Jaccard.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$jaccard.index.tune})),na.rm=T)
        Rsquare.outcome.matrix.train.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$R2.outcome.train.tune})),na.rm=T)
        Rsquare.outcome.matrix.test.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$R2.outcome.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.test.tune})),na.rm=T)
        
        w.yujia.matrix[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$w0.est})),na.rm=T)
        w.yujia.matrix.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$w0.est})))/sqrt(50)
        
        ari.train.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$ari.train.tune})),na.rm=T)
        ari.test.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$ari.test.tune})),na.rm=T)
        Rmse.train.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$Rmse.train.tune})),na.rm=T)
        Rmse.test.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$Rmse.test.tune})),na.rm=T)
        Rmse.train.Peng.matrix.num.hard[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$Rmse.train.tune.hard})),na.rm=T)
        Rmse.test.Peng.matrix.num.hard[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$Rmse.test.tune.hard})),na.rm=T)
        Jaccard.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$jaccard.index.tune})),na.rm=T)
        num.Peng.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$GM.res$num.tune})))
        Rsquare.outcome.matrix.train.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$R2.outcome.train.tune})),na.rm=T)
        Rsquare.outcome.matrix.test.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$R2.outcome.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$R2.gene.test.tune})),na.rm=T)
        
        ari.train.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$ari.train.tune})),na.rm=T)
        ari.test.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$ari.test.tune})),na.rm=T)
        Rmse.train.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$Rmse.train.tune})),na.rm=T)
        Rmse.test.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$Rmse.test.tune})),na.rm=T)
        Rmse.train.Peng.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$Rmse.train.tune.hard})),na.rm=T)
        Rmse.test.Peng.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$Rmse.test.tune.hard})),na.rm=T)
        Jaccard.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$jaccard.index.tune})),na.rm=T)
        Rsquare.outcome.matrix.train.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$R2.outcome.train.tune})),na.rm=T)
        Rsquare.outcome.matrix.test.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$R2.outcome.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$R2.gene.test.tune})),na.rm=TRUE)
        
        
        ari.train.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$ari.train.tune})),na.rm=T)
        ari.test.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$ari.test.tune})),na.rm=T)
        Rmse.train.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$Rmse.train.tune})),na.rm=T)
        Rmse.test.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$Rmse.test.tune})),na.rm=T)
        Jaccard.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$jaccard.index.tune})),na.rm=T)
        num.skmeans.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$SKM.res$num.tune})))
        Rsquare.outcome.matrix.train.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$R2.outcome.train.tune})),na.rm=T)
        Rsquare.outcome.matrix.test.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$R2.outcome.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.test.tune})),na.rm=T)
        
        ari.train.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$ari.train.tune})),na.rm=T)
        ari.test.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$ari.test.tune})),na.rm=T)
        Rmse.train.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$Rmse.train.tune})),na.rm=T)
        Rmse.test.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$Rmse.test.tune})),na.rm=T)
        Jaccard.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$jaccard.index.tune})),na.rm=T)
        Rsquare.outcome.matrix.train.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$R2.outcome.train.tune})),na.rm=T)
        Rsquare.outcome.matrix.test.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$R2.outcome.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.test.tune})),na.rm=T)
        
        ari.train.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$ari.train.tune})),na.rm=T)
        ari.test.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$ari.test.tune})),na.rm=T)
        Rmse.train.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$Rmse.train.tune})),na.rm=T)
        Rmse.test.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$Rmse.test.tune})),na.rm=T)
        Jaccard.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$jaccard.index.tune})),na.rm=T)
        num.PMBC.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$PMBC.res$num.tune})))
        Rsquare.outcome.matrix.train.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$R2.outcome.train.tune})),na.rm=T)
        Rsquare.outcome.matrix.test.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$R2.outcome.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.test.tune})),na.rm=T)
        
        ari.train.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$ari.train.tune})),na.rm=T)
        ari.test.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$ari.test.tune})),na.rm=T)
        Rmse.train.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$Rmse.train.tune})),na.rm=T)
        Rmse.test.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$Rmse.test.tune})),na.rm=T)
        Jaccard.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$jaccard.index.tune})),na.rm=T)
        Rsquare.outcome.matrix.train.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$R2.outcome.train.tune})),na.rm=T)
        Rsquare.outcome.matrix.test.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$R2.outcome.test.tune})),na.rm=T)
        Rsquare.gene.matrix.train.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.train.tune})),na.rm=T)
        Rsquare.gene.matrix.test.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.test.tune})),na.rm=T)
        
      }
      else {
        ari.train.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$ari.train.num})),na.rm=T)
        ari.test.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$ari.test.num})),na.rm=T)
        Rmse.train.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$Rmse.train.num})),na.rm=T)
        Rmse.test.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$Rmse.test.num})),na.rm=T)
        Rmse.train.yujia.matrix.num.hard[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$Rmse.train.num.hard})),na.rm=T)
        Rmse.test.yujia.matrix.num.hard[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$Rmse.test.num.hard})),na.rm=T)
        Jaccard.yujia.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$jaccard.index.num})),na.rm=T)
        num.yujia.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$WJL.res$num.num})))
        Rsquare.outcome.matrix.train[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$R2.outcome.train})),na.rm=T)
        Rsquare.outcome.matrix.test[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$R2.outcome.test})),na.rm=T)
        Rsquare.gene.matrix.train[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.train})),na.rm=T)
        Rsquare.gene.matrix.test[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.test})),na.rm=T)
        
        ari.train.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$ari.train.num})),na.rm=T)
        ari.test.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$ari.test.num})),na.rm=T)
        Rmse.train.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$Rmse.train.num})),na.rm=T)
        Rmse.test.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$Rmse.test.num})),na.rm=T)
        Rmse.train.yujia.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$Rmse.train.num.hard})),na.rm=T)
        Rmse.test.yujia.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$Rmse.test.num.hard})),na.rm=T)
        Jaccard.yujia.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$jaccard.index.num})),na.rm=T)
        Rsquare.outcome.matrix.train.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$R2.outcome.train})),na.rm=T)
        Rsquare.outcome.matrix.test.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$R2.outcome.test})),na.rm=T)
        Rsquare.gene.matrix.train.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.train})),na.rm=T)
        Rsquare.gene.matrix.test.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$R2.gene.test})),na.rm=T)
        
        w.yujia.matrix[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$WJL.res$w0.est})),na.rm=T)
        w.yujia.matrix.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$WJL.res$w0.est})))/sqrt(50)
        
        ari.train.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$ari.train.num})),na.rm=T)
        ari.test.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$ari.test.num})),na.rm=T)
        Rmse.train.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$Rmse.train.num})),na.rm=T)
        Rmse.test.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$Rmse.test.num})),na.rm=T)
        Rmse.train.Peng.matrix.num.hard[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$Rmse.train.num.hard})),na.rm=T)
        Rmse.test.Peng.matrix.num.hard[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$Rmse.test.num.hard})),na.rm=T)
        Jaccard.Peng.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$jaccard.index.num})),na.rm=T)
        num.Peng.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$GM.res$num.num})))
        Rsquare.outcome.matrix.train.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$R2.outcome.train})),na.rm=T)
        Rsquare.outcome.matrix.test.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$R2.outcome.test})),na.rm=T)
        Rsquare.gene.matrix.train.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$R2.gene.train})),na.rm=T)
        Rsquare.gene.matrix.test.Peng[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$GM.res$R2.gene.test})),na.rm=T)
        
        ari.train.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$ari.train.num})),na.rm=T)
        ari.test.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$ari.test.num})),na.rm=T)
        Rmse.train.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$Rmse.train.num})),na.rm=T)
        Rmse.test.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$Rmse.test.num})),na.rm=T)
        Rmse.train.Peng.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$Rmse.train.num.hard})),na.rm=T)
        Rmse.test.Peng.matrix.num.hard.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$Rmse.test.num.hard})),na.rm=T)
        Jaccard.Peng.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$jaccard.index.num})),na.rm=T)
        Rsquare.outcome.matrix.train.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$R2.outcome.train})),na.rm=T)
        Rsquare.outcome.matrix.test.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$R2.outcome.test})),na.rm=T)
        Rsquare.gene.matrix.train.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$R2.gene.train})),na.rm=T)
        Rsquare.gene.matrix.test.Peng.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$GM.res$R2.gene.test})),na.rm=TRUE)
        
        
        ari.train.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$ari.train.num})),na.rm=T)
        ari.test.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$ari.test.num})),na.rm=T)
        Rmse.train.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$Rmse.train.num})),na.rm=T)
        Rmse.test.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$Rmse.test.num})),na.rm=T)
        Jaccard.skmeans.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$jaccard.index.num})),na.rm=T)
        num.skmeans.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$SKM.res$num.num})))
        Rsquare.outcome.matrix.train.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$R2.outcome.train})),na.rm=T)
        Rsquare.outcome.matrix.test.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$R2.outcome.test})),na.rm=T)
        Rsquare.gene.matrix.train.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.train})),na.rm=T)
        Rsquare.gene.matrix.test.skmeans[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.test})),na.rm=T)
        
        ari.train.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$ari.train.num})),na.rm=T)
        ari.test.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$ari.test.num})),na.rm=T)
        Rmse.train.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$Rmse.train.num})),na.rm=T)
        Rmse.test.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$Rmse.test.num})),na.rm=T)
        Jaccard.skmeans.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$jaccard.index.num})),na.rm=T)
        Rsquare.outcome.matrix.train.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$R2.outcome.train})),na.rm=T)
        Rsquare.outcome.matrix.test.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$R2.outcome.test})),na.rm=T)
        Rsquare.gene.matrix.train.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.train})),na.rm=T)
        Rsquare.gene.matrix.test.skmeans.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$SKM.res$R2.gene.test})),na.rm=T)
        
        ari.train.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$ari.train.num})),na.rm=T)
        ari.test.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$ari.test.num})),na.rm=T)
        Rmse.train.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$Rmse.train.num})),na.rm=T)
        Rmse.test.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$Rmse.test.num})),na.rm=T)
        Jaccard.PMBC.matrix.num[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$jaccard.index.num})),na.rm=T)
        num.PMBC.matrix.num[ind_c1,ind_mu]<-median(unlist(lapply(RES,function(x){x$PMBC.res$num.num})))
        Rsquare.outcome.matrix.train.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$R2.outcome.train})),na.rm=T)
        Rsquare.outcome.matrix.test.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$R2.outcome.test})),na.rm=T)
        Rsquare.gene.matrix.train.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.train})),na.rm=T)
        Rsquare.gene.matrix.test.PMBC[ind_c1,ind_mu]<-mean(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.test})),na.rm=T)
        
        ari.train.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$ari.train.num})),na.rm=T)
        ari.test.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$ari.test.num})),na.rm=T)
        Rmse.train.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$Rmse.train.num})),na.rm=T)
        Rmse.test.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$Rmse.test.num})),na.rm=T)
        Jaccard.PMBC.matrix.num.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$jaccard.index.num})),na.rm=T)
        Rsquare.outcome.matrix.train.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$R2.outcome.train})),na.rm=T)
        Rsquare.outcome.matrix.test.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$R2.outcome.test})),na.rm=T)
        Rsquare.gene.matrix.train.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.train})),na.rm=T)
        Rsquare.gene.matrix.test.PMBC.se[ind_c1,ind_mu]<-sd(unlist(lapply(RES,function(x){x$PMBC.res$R2.gene.test})),na.rm=T)
        
      }
      
    }
    
    ari.train1<-rbind(mu_vector,
                      ari.train.yujia.matrix.num[ind_c1,],ari.train.Peng.matrix.num[ind_c1,],
                      ari.train.skmeans.matrix.num[ind_c1,],ari.train.PMBC.matrix.num[ind_c1,])
    ari.test1<-rbind(mu_vector,
                     ari.test.yujia.matrix.num[ind_c1,],ari.test.Peng.matrix.num[ind_c1,],
                     ari.test.skmeans.matrix.num[ind_c1,],ari.test.PMBC.matrix.num[ind_c1,])
    Rmse.train1<-rbind(mu_vector,
                       Rmse.train.yujia.matrix.num.hard[ind_c1,],Rmse.train.Peng.matrix.num.hard[ind_c1,],
                       Rmse.train.skmeans.matrix.num[ind_c1,],Rmse.train.PMBC.matrix.num[ind_c1,])
    Rmse.test1<-rbind(mu_vector,
                      Rmse.test.yujia.matrix.num.hard[ind_c1,],Rmse.test.Peng.matrix.num.hard[ind_c1,],
                      Rmse.test.skmeans.matrix.num[ind_c1,],Rmse.test.PMBC.matrix.num[ind_c1,])
    jaccard1<-rbind(mu_vector,
                    Jaccard.yujia.matrix.num[ind_c1,],Jaccard.Peng.matrix.num[ind_c1,],
                    Jaccard.skmeans.matrix.num[ind_c1,],Jaccard.PMBC.matrix.num[ind_c1,])
    Rsquare.outcome.train1<-rbind(mu_vector,
                                  Rsquare.outcome.matrix.train[ind_c1,],Rsquare.outcome.matrix.train.Peng[ind_c1,],
                                  Rsquare.outcome.matrix.train.skmeans[ind_c1,],Rsquare.outcome.matrix.train.PMBC[ind_c1,])
    
    Rsquare.gene.train1<-rbind(mu_vector,
                               Rsquare.gene.matrix.train[ind_c1,],Rsquare.gene.matrix.train.Peng[ind_c1,],
                               Rsquare.gene.matrix.train.skmeans[ind_c1,],Rsquare.gene.matrix.train.PMBC[ind_c1,])
    
    Rsquare.outcome.test1<-rbind(mu_vector,
                                 Rsquare.outcome.matrix.test[ind_c1,],Rsquare.outcome.matrix.test.Peng[ind_c1,],
                                 Rsquare.outcome.matrix.test.skmeans[ind_c1,],Rsquare.outcome.matrix.test.PMBC[ind_c1,])
    
    Rsquare.gene.test1<-rbind(mu_vector,
                              Rsquare.gene.matrix.test[ind_c1,],Rsquare.gene.matrix.test.Peng[ind_c1,],
                              Rsquare.gene.matrix.test.skmeans[ind_c1,],Rsquare.gene.matrix.test.PMBC[ind_c1,])
    
    ari.train1.se<-rbind(mu_vector,
                         ari.train.yujia.matrix.num.se[ind_c1,],ari.train.Peng.matrix.num.se[ind_c1,],
                         ari.train.skmeans.matrix.num.se[ind_c1,],Rsquare.gene.matrix.train.PMBC[ind_c1,])
    ari.test1.se<-rbind(mu_vector,
                        ari.test.yujia.matrix.num.se[ind_c1,],ari.test.Peng.matrix.num.se[ind_c1,],
                        ari.test.skmeans.matrix.num.se[ind_c1,],Rsquare.gene.matrix.test.PMBC[ind_c1,])
    
    Rmse.train1.se<-rbind(mu_vector,
                          Rmse.train.yujia.matrix.num.hard.se[ind_c1,],Rmse.train.Peng.matrix.num.hard.se[ind_c1,],
                          Rmse.train.skmeans.matrix.num.se[ind_c1,],Rmse.train.PMBC.matrix.num.se[ind_c1,])
    Rmse.test1.se<-rbind(mu_vector,
                         Rmse.test.yujia.matrix.num.hard.se[ind_c1,],Rmse.test.Peng.matrix.num.hard.se[ind_c1,],
                         Rmse.test.skmeans.matrix.num.se[ind_c1,],Rmse.test.PMBC.matrix.num.se[ind_c1,])
    jaccard1.se<-rbind(mu_vector,
                       Jaccard.yujia.matrix.num.se[ind_c1,],Jaccard.Peng.matrix.num.se[ind_c1,],
                       Jaccard.skmeans.matrix.num.se[ind_c1,],Jaccard.PMBC.matrix.num.se[ind_c1,])
    
    Rsquare.outcome.train1.se<-rbind(mu_vector,
                                     Rsquare.outcome.matrix.train.se[ind_c1,],Rsquare.outcome.matrix.train.Peng.se[ind_c1,],
                                     Rsquare.outcome.matrix.train.skmeans.se[ind_c1,],Rsquare.outcome.matrix.train.PMBC.se[ind_c1,])
    
    Rsquare.gene.train1.se<-rbind(mu_vector,
                                  Rsquare.gene.matrix.train.se[ind_c1,],Rsquare.gene.matrix.train.Peng.se[ind_c1,],
                                  Rsquare.gene.matrix.train.skmeans.se[ind_c1,], Rsquare.gene.matrix.train.PMBC.se[ind_c1,])
    
    Rsquare.outcome.test1.se<-rbind(mu_vector,
                                    Rsquare.outcome.matrix.test.se[ind_c1,],Rsquare.outcome.matrix.test.Peng.se[ind_c1,],
                                    Rsquare.outcome.matrix.test.skmeans.se[ind_c1,],Rsquare.outcome.matrix.test.PMBC.se[ind_c1,])
    
    Rsquare.gene.test1.se<-rbind(mu_vector,
                                 Rsquare.gene.matrix.test.se[ind_c1,],Rsquare.gene.matrix.test.Peng.se[ind_c1,],
                                 Rsquare.gene.matrix.test.skmeans.se[ind_c1,],Rsquare.gene.matrix.test.PMBC.se[ind_c1,])
    
    
    rownames(ari.train1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(ari.test1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rmse.train1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rmse.test1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(jaccard1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rsquare.outcome.train1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rsquare.gene.train1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rsquare.outcome.test1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rsquare.gene.test1)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    
    rownames(ari.train1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(ari.test1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rmse.train1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rmse.test1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(jaccard1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rsquare.outcome.train1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rsquare.gene.train1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rsquare.outcome.test1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    rownames(Rsquare.gene.test1.se)<-c("Mu","ogClust_WJL","ogClust_GM","SKM","PMBC")
    
    
    ari.train1.se<-pivot_longer(data = as.data.frame(t(ari.train1.se)),cols=!Mu,names_to="method",values_to="value")
    ari.test1.se<-pivot_longer(data = as.data.frame(t(ari.test1.se)),cols=!Mu,names_to="method",values_to="value")
    Rmse.train1.se<-pivot_longer(data = as.data.frame(t(Rmse.train1.se)),cols=!Mu,names_to="method",values_to="value")
    Rmse.test1.se<-pivot_longer(data = as.data.frame(t(Rmse.test1.se)),cols=!Mu,names_to="method",values_to="value")
    jaccard1.se<-pivot_longer(data = as.data.frame(t(jaccard1.se)),cols=!Mu,names_to="method",values_to="value")
    Rsquare.outcome.train1.se<-pivot_longer(data = as.data.frame(t(Rsquare.outcome.train1.se)),cols=!Mu,names_to="method",values_to="value")
    Rsquare.gene.train1.se<-pivot_longer(data = as.data.frame(t(Rsquare.gene.train1.se)),cols=!Mu,names_to="method",values_to="value")
    Rsquare.outcome.test1.se<-pivot_longer(data = as.data.frame(t(Rsquare.outcome.test1.se)),cols=!Mu,names_to="method",values_to="value")
    Rsquare.gene.test1.se<-pivot_longer(data = as.data.frame(t(Rsquare.gene.test1.se)),cols=!Mu,names_to="method",values_to="value")
    
    
    ari.train1<-pivot_longer(data = as.data.frame(t(ari.train1)),cols=!Mu,names_to="method",values_to="value")
    ari.test1<-pivot_longer(data = as.data.frame(t(ari.test1)),cols=!Mu,names_to="method",values_to="value")
    Rmse.train1<-pivot_longer(data = as.data.frame(t(Rmse.train1)),cols=!Mu,names_to="method",values_to="value")
    Rmse.test1<-pivot_longer(data = as.data.frame(t(Rmse.test1)),cols=!Mu,names_to="method",values_to="value")
    jaccard1<-pivot_longer(data = as.data.frame(t(jaccard1)),cols=!Mu,names_to="method",values_to="value")
    Rsquare.outcome.train1<-pivot_longer(data = as.data.frame(t(Rsquare.outcome.train1)),cols=!Mu,names_to="method",values_to="value")
    Rsquare.gene.train1<-pivot_longer(data = as.data.frame(t(Rsquare.gene.train1)),cols=!Mu,names_to="method",values_to="value")
    Rsquare.outcome.test1<-pivot_longer(data = as.data.frame(t(Rsquare.outcome.test1)),cols=!Mu,names_to="method",values_to="value")
    Rsquare.gene.test1<-pivot_longer(data = as.data.frame(t(Rsquare.gene.test1)),cols=!Mu,names_to="method",values_to="value")
    
    ari.train1$se<-ari.train1.se$value/sqrt(50)
    ari.test1$se<-ari.test1.se$value/sqrt(50)
    Rmse.train1$se<-Rmse.train1.se$value/sqrt(50)
    Rmse.test1$se<-Rmse.test1.se$value/sqrt(50)
    jaccard1$se<-jaccard1.se$value/sqrt(50)
    Rsquare.outcome.train1$se<-Rsquare.outcome.train1.se$value/sqrt(50)
    Rsquare.outcome.test1$se<-Rsquare.outcome.test1.se$value/sqrt(50)
    Rsquare.gene.train1$se<-Rsquare.gene.train1.se$value/sqrt(50)
    Rsquare.gene.test1$se<-Rsquare.gene.test1.se$value/sqrt(50)
    
    
    Group<-c(rep("ARI.train",nrow(ari.train1)),rep("ARI.test",nrow(ari.train1)),rep("Rmse.train",nrow(ari.train1)),
             rep("Rmse.test",nrow(ari.train1)),rep("Jaccard",nrow(ari.train1)),
             rep("R2.outcome.train",nrow(ari.train1)),rep("R2.gene.train",nrow(ari.train1)),
             rep("R2.outcome.test",nrow(ari.train1)),rep("R2.gene.test",nrow(ari.train1)))
    data1<-rbind(ari.train1,ari.test1,Rmse.train1,Rmse.test1,jaccard1,
                 Rsquare.outcome.train1,Rsquare.gene.train1,Rsquare.outcome.test1,Rsquare.gene.test1)  
    data1<-cbind(data1,Group)
    data1$Group<-factor(data1$Group,levels = c("ARI.train","Rmse.train","Jaccard",
                                               "ARI.test","Rmse.test",
                                               "R2.outcome.train","R2.gene.train",
                                               "R2.outcome.test","R2.gene.test"))
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
    
    data<-data1[which(data1$Group=="Rmse.train"),]
    p2<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab("RMSE")
    p2
    
    data<-data1[which(data1$Group=="R2.outcome.train"),]
    p3<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab(expression(paste(italic(R)^2, "(outcome)",sep="")))
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
    
    data<-data1[which(data1$Group=="Rmse.test"),]
    p6<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab("RMSE")
    p6
    
    data<-data1[which(data1$Group=="R2.outcome.test"),]
    p7<-ggplot(data)+aes(x=Mu,y=value,color=method,shape=method)+geom_point(size=2.5, position = position_dodge(0.1))+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.6,width =.2,lwd=0.3,position=position_dodge(0.1))+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
            strip.text.x = element_text(size = 20),legend.position = 'none',
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      scale_x_continuous(breaks=mu_vector)+xlab(expression(paste('effect size of gene ',"(",mu,")")))+ylab(expression(paste(italic(R)^2, "(outcome)",sep="")))
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

#generate Figrue I
summary_gen(SC=FALSE,tune=TRUE,c1_vector_all=c(3),figure.name = "FigureI.png")
#generate Supplementary Figrue II
summary_gen(SC=FALSE,tune=TRUE,c1_vector_all=c(1),figure.name="SupplementFigrueII.png")
#generate Supplementary Figrue III
summary_gen(SC=FALSE,tune=FALSE,c1_vector_all=c(3),figure.name="SupplementFigrueIII.png")
#generate Supplementary Figrue IV
summary_gen(SC=TRUE,tune=TRUE,c1_vector_all=c(3),figure.name="SupplementFigrueIV.png")


