
###ogClust_WJL method #####
ogclust_con_Yujia<-function(x,G,y,c_center=NULL,lambda,v_int=NULL,pi_int=NULL,K=3,max_iter=200,w_outcome,w_G,z_int=NULL){
  x.origin<-x
  mult_density1<-function(x,mu,sigma){

    lik<-dnorm(x,mean=mu,sd=sqrt(sigma),log=TRUE)
    return(lik)
  }

  #Column refers to samples and row refers to genes

  p<-dim(G)[1] # number of variables
  n<-dim(G)[2] # number of subjects

  #========== E-step: initialization===========#

  mu_int=c_center



  if(length(v_int)!=p){
    v_int<-rep(1,p)
  }
  if(length(pi_int)!=p){
    pi_int<-rep(1/K,K)
  }
  pi<-pi_int
  v<-v_int
  mu<-mu_int


  if(is.null(z_int)){
    #-----------------------------------Set initial clustering assignment
    z_int<-matrix(,nrow=n,ncol=K) # initial value of prob in cluster k for each subject
    mult_pdf<-matrix(,nrow=n,ncol=K)
    #------------------------------each sample in each cluster (sum of log likelihood)
    for(j in 1:K){
      temp<-mult_density1(as.numeric(G),mu=rep(mu_int[,j],times=n),sigma=rep(v_int,times=n))
      temp_matrix<-matrix(temp,nrow=p,ncol=n,byrow=FALSE)
      mult_pdf[,j]<-t(temp_matrix)%*%rep(1,p)
    }

    #---------------------The initial clustering assignment is by Gene expression only
    max_pdf<-apply(mult_pdf,1,max)
    max_pdf_matrix<-matrix(rep(max_pdf,times=K),ncol=K)
    mult_pdf1<-mult_pdf-max_pdf_matrix
    mult_pdf2<-exp(mult_pdf1)*matrix(rep(pi_int,times=n),byrow=T,ncol=K)
    sum_mult_pdf2<-mult_pdf2%*%rep(1,K)
    z_int<-mult_pdf2/matrix(rep(sum_mult_pdf2,times=K),ncol=K)

    #initialize coefficient in outcome association
    #x<-cbind(1,x)
    #----------------------------------------
    x<-matrix(NA,nrow=nrow(x.origin)*K,ncol=ncol(x.origin))
    intercept<-matrix(0,nrow=n*K,ncol=K)
    for(i in 1:ncol(x.origin)){
      x[,i]<-rep(x.origin[,i],K)
    }
    for(i in 1:K){
      index<-((i-1)*n+1):(i*n)
      intercept[index,i]<-1
    }
    colnames(x)<-colnames(x.origin)
    colnames(intercept)<-paste("intercept",1:K,sep="")
    x<-cbind(intercept,x)
    W<-as.numeric(z_int)
    W[which(W == 0)] = 10^(-100)
    y1<-rep(y,K)
    data1<-cbind(y1,x)
    data1<-as.data.frame(data1)
    mod<-lm(y1~ -1+., data1,weights = W)
    #mod$residuals[1]
    int_coef<-mod$coefficients
    #solve(t(x)%*%diag(W)%*%x)%*%t(x)%*%diag(W)%*%y1
    #(y1-x%*%int_coef)[1]
    #int_sigma_coef<-mod$scale
    #mod$fitted.values

    int_sigma_coef<-t(y1-x%*%int_coef)%*%diag(W)%*%(y1-x%*%int_coef)/n
    int_sigma_coef<-sqrt(int_sigma_coef)
    x<-x.origin
    #----------------------------------------

    #update the z matrix
    z_outcome<-matrix(,nrow=n,ncol=K)
    for(j in 1:K){
      z_outcome[,j]<-dnorm(y,mean=int_coef[j]+x.origin%*%int_coef[(K+1):length(int_coef)],sd=int_sigma_coef,log=T)
      #z_outcome[,j]<-dnorm(y,mean=int_coef[j]+int_coef[4]*x.origin[,1]+int_coef[5]*x.origin[,2],sd=int_sigma_coef,log=T)
    }
    #z_outcome<-log(z_outcome)
    #mult_pdf<-mult_pdf*s_G
    #z_outcome<-z_outcome/s_G
    mult_pdf.add<-w_G*mult_pdf+w_outcome*z_outcome

    max_pdf<-apply(mult_pdf.add,1,max)
    max_pdf_matrix<-matrix(rep(max_pdf,times=K),ncol=K)
    mult_pdf1<-mult_pdf.add-max_pdf_matrix
    mult_pdf2<-exp(mult_pdf1)*matrix(rep(pi_int,times=n),byrow=T,ncol=K)
    sum_mult_pdf2<-mult_pdf2%*%rep(1,K)
    z<-mult_pdf2/matrix(rep(sum_mult_pdf2,times=K),ncol=K)

  }else{
    z<-z_int
  }

  pi_up<-apply(z,2,sum)/n
  v_up_matrix<-matrix(NA,ncol=K,nrow=p)
  for(j in 1:K){
    temp1<-(G-matrix(rep(mu[,j],times=n),ncol=n))^2
    v_up_matrix[,j]<-temp1%*%z[,j]/n
  }
  v_up<-v_up_matrix%*%rep(1,K)
  #update mu
  mu_tu<-matrix(,nrow=p,ncol=K)
  mu_up<-matrix(,nrow=p,ncol=K)
  for(i in 1:K){
    #temp<-sapply(1:n,function(x) z[x,i]*data[,x])
    G1<-G*matrix(rep(z[,i],times=p),ncol=n,byrow=T)
    mu_tu[,i]<-G1%*%rep(1,n)/sum(z[,i])
    #mu_tu[,i]<-apply(temp,1,sum)/sum(z[,i])
    #mu_up[,i]<-ifelse(lambda<=abs((G1%*%rep(1,n))/v),sign(mu_tu[,i])*(abs(mu_tu[,i])-(lambda/sum(z[,i]))*v),0)
    mu_up[,i]<-ifelse(lambda<=abs((G1%*%rep(1,n))/v),sign(mu_tu[,i])*(abs(mu_tu[,i])-(lambda/sum(z[,i]))*v),0)
  }

  mult_pdf<-matrix(,nrow=n,ncol=K)
  for(j in 1:K){
    temp<-mult_density1(as.numeric(G),mu=rep(mu_up[,j],times=n),sigma=rep(v_up,times=n))
    temp_matrix<-matrix(temp,nrow=p,ncol=n,byrow=FALSE)
    mult_pdf[,j]<-t(temp_matrix)%*%rep(1,p)
  }

  x<-matrix(NA,nrow=nrow(x.origin)*K,ncol=ncol(x.origin))
  intercept<-matrix(0,nrow=n*K,ncol=K)
  for(i in 1:ncol(x)){
    x[,i]<-rep(x.origin[,i],K)
  }
  for(i in 1:K){
    index<-((i-1)*n+1):(i*n)
    intercept[index,i]<-1
  }
  colnames(x)<-colnames(x.origin)
  colnames(intercept)<-paste("intercept",1:K,sep="")
  x<-cbind(intercept,x)
  W<-as.numeric(z)
  W[which(W == 0)] = 10^(-100)
  y1<-rep(y,K)
  data1<-cbind(y1,x)
  data1<-as.data.frame(data1)
  mod<-lm(y1~ -1+., data1,weights = W)
  int_coef<-mod$coefficients
  #int_sigma_coef<-mod$scale
  #mod$fitted.values
  #sqrt(sum((y1[1:33]-x[1:33,4:5]%*%int_coef[4:5]-int_coef[2])^2)/33)
  #sqrt(sum((y1[34:66]-x[34:66,4:5]%*%int_coef[4:5]-int_coef[3])^2)/33)
  #sqrt(sum((y1[67:99]-x[67:99,4:5]%*%int_coef[4:5]-int_coef[1])^2)/33)
  #sqrt(sum(W*(mod$residuals^2))/(900))
  int_sigma_coef<-t(y1-x%*%int_coef)%*%diag(W)%*%(y1-x%*%int_coef)/(n)
  int_sigma_coef<-sqrt(int_sigma_coef)
  x<-x.origin
  #----------------------------------------

  #update the z matrix
  z_outcome<-matrix(,nrow=n,ncol=K)
  for(j in 1:K){
    z_outcome[,j]<-dnorm(y,mean=int_coef[j]+x.origin%*%int_coef[(K+1):length(int_coef)],sd=int_sigma_coef,log=T)
    #z_outcome[,j]<-dnorm(y,mean=int_coef[j]+int_coef[4]*x.origin[,1]+int_coef[5]*x.origin[,2],sd=int_sigma_coef,log=T)
  }
  #z_outcome[114,j]<-dnorm(y[114],mean=int_coef[j]+x.origin[114,]%*%int_coef[(K+1):length(int_coef)],sd=int_sigma_coef,log=T)
  #z_outcome<-log(z_outcome)
  #mult_pdf<-mult_pdf*s_G
  #z_outcome<-z_outcome/s_G
  mult_pdf.add<-w_G*mult_pdf+w_outcome*z_outcome

  max_pdf<-apply(mult_pdf.add,1,max)
  max_pdf_matrix<-matrix(rep(max_pdf,times=K),ncol=K)
  mult_pdf1<-mult_pdf.add-max_pdf_matrix
  mult_pdf2<-exp(mult_pdf1)*matrix(rep(pi_int,times=n),byrow=T,ncol=K)
  sum_mult_pdf2<-mult_pdf2%*%rep(1,K)
  z_up<-mult_pdf2/matrix(rep(sum_mult_pdf2,times=K),ncol=K)

  # update parameters values
  pi<-pi_up
  v<-v_up
  mu<-mu_up
  z<-z_up
  #========= M step ==========#
  iter=1
  log_lik<-1 # initialize
  log_lik_up<-0 #initialize

  while(abs(log_lik_up-log_lik)>10^(-7) & iter <=max_iter){
    #print(iter)
    if(sum(pi==0)!=0){
      pi[which(pi==0)]<-10^(-6)
    }
    log_lik<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf.add[x,])))-lambda*sum(abs(mu))
    pi_up<-apply(z,2,sum)/n
    v_up_matrix<-matrix(NA,ncol=K,nrow=p)
    for(j in 1:K){
      temp1<-(G-matrix(rep(mu[,j],times=n),ncol=n))^2
      v_up_matrix[,j]<-temp1%*%z[,j]/n
    }
    v_up<-v_up_matrix%*%rep(1,K)
    #update mu
    mu_tu<-matrix(,nrow=p,ncol=K)
    mu_up<-matrix(,nrow=p,ncol=K)
    for(i in 1:K){
      #temp<-sapply(1:n,function(x) z[x,i]*data[,x])
      G1<-G*matrix(rep(z[,i],times=p),ncol=n,byrow=T)
      mu_tu[,i]<-G1%*%rep(1,n)/sum(z[,i])
      #mu_tu[,i]<-apply(temp,1,sum)/sum(z[,i])
      #mu_up[,i]<-ifelse(lambda<=abs((G1%*%rep(1,n))/v),sign(mu_tu[,i])*(abs(mu_tu[,i])-(lambda/sum(z[,i]))*v),0)
      mu_up[,i]<-ifelse(lambda<=abs((G1%*%rep(1,n))/v),sign(mu_tu[,i])*(abs(mu_tu[,i])-(lambda/sum(z[,i]))*v),0)
    }

    mult_pdf<-matrix(,nrow=n,ncol=K)
    for(j in 1:K){
      temp<-mult_density1(as.numeric(G),mu=rep(mu_up[,j],times=n),sigma=rep(v_up,times=n))
      temp_matrix<-matrix(temp,nrow=p,ncol=n,byrow=FALSE)
      mult_pdf[,j]<-t(temp_matrix)%*%rep(1,p)
    }

    x<-matrix(NA,nrow=nrow(x.origin)*K,ncol=ncol(x.origin))
    intercept<-matrix(0,nrow=n*K,ncol=K)
    for(i in 1:ncol(x)){
      x[,i]<-rep(x.origin[,i],K)
    }
    for(i in 1:K){
      index<-((i-1)*n+1):(i*n)
      intercept[index,i]<-1
    }
    colnames(x)<-colnames(x.origin)
    colnames(intercept)<-paste("intercept",1:K,sep="")
    x<-cbind(intercept,x)
    W<-as.numeric(z)
    W[which(W == 0)] = 10^(-100)
    y1<-rep(y,K)
    data1<-cbind(y1,x)
    data1<-as.data.frame(data1)
    mod<-lm(y1~ -1+., data1,weights = W)
    int_coef<-mod$coefficients
    #int_sigma_coef<-mod$scale
    #mod$fitted.values

    int_sigma_coef<-t(y1-x%*%int_coef)%*%diag(W)%*%(y1-x%*%int_coef)/n
    int_sigma_coef<-sqrt(int_sigma_coef)
    x<-x.origin
    #----------------------------------------

    #update the z matrix
    z_outcome<-matrix(,nrow=n,ncol=K)
    for(j in 1:K){
      z_outcome[,j]<-dnorm(y,mean=int_coef[j]+x.origin%*%int_coef[(K+1):length(int_coef)],sd=int_sigma_coef,log=T)
      #z_outcome[,j]<-dnorm(y,mean=int_coef[j]+int_coef[4]*x.origin[,1]+int_coef[5]*x.origin[,2],sd=int_sigma_coef,log=T)
    }

    #z_outcome<-log(z_outcome)
    #z_outcome<-log(z_outcome)
    #mult_pdf<-mult_pdf*s_G
    #z_outcome<-z_outcome/s_G
    mult_pdf.add<-w_G*mult_pdf+w_outcome*z_outcome

    max_pdf<-apply(mult_pdf.add,1,max)
    max_pdf_matrix<-matrix(rep(max_pdf,times=K),ncol=K)
    mult_pdf1<-mult_pdf.add-max_pdf_matrix
    mult_pdf2<-exp(mult_pdf1)*matrix(rep(pi_int,times=n),byrow=T,ncol=K)
    sum_mult_pdf2<-mult_pdf2%*%rep(1,K)
    z_up<-mult_pdf2/matrix(rep(sum_mult_pdf2,times=K),ncol=K)

    # update parameters values
    pi<-pi_up
    v<-v_up
    mu<-mu_up
    z<-z_up

    # if(sum(pi==0)!=0){
    #   pi[which(pi==0)]<-10^(-100)
    # }
    log_lik_up<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf.add[x,])))-lambda*sum(abs(mu))
    unpen_lik_up<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf.add[x,])))
    log_lik_uncomplete<-0
    max_mult<-apply(mult_pdf.add,1,max)
    for(i in 1:nrow(mult_pdf)){
      log_lik_uncomplete<-log_lik_uncomplete+sum(pi*exp(mult_pdf.add[i,]-max_mult[i]))+max_mult[i]
      #log_lik<-log_lik+log(sum(pi*exp(mult_pdf[i,]-max_mult[i])))+max_mult[i]
    }
    iter=iter+1
  }
  result_list <- list('mu'=mu,'sigma'=v,'z'=z,'pi'=pi,
                      'int_coef'=int_coef,'int_sigma_coef'=int_sigma_coef)
  log_lik<-log_lik_up
  unpen_lik<-unpen_lik_up
  max_lik<-log_lik_uncomplete
  max_mu<-result_list$mu
  s<-apply(max_mu,1,function(x) sum(x==0))
  q<-sum(s)
  #BIC
  if(lambda!=0) {
    BIC<--2*max_lik+log(n)*(K+p+K*p-q)
  } else if(lambda==0) {
    BIC<--2*max_lik+log(n)*(K+p+K*p)
  }

  res<-list('result_list'=result_list,'BIC'=BIC,'lik'=max_lik,'iter'=iter)
  return(res)

}

Fstat<-function(cluster,Y,G){
  #cluster<-apply(mod$result_list$z,1,which.max)
  F.gene<-c()
  for(i in 1:nrow(G)){
    data1<-data.frame(gene=G[i,],cluster=as.factor(cluster))
    res.aov <- aov(gene ~ cluster, data = data1)
    res.aov<-summary(res.aov)
    F.gene[i]<-res.aov[[1]]$`F value`[1]
  }
  data1<-data.frame(outcome=Y,cluster=as.factor(cluster))
  res.aov <- aov(outcome ~ cluster, data = data1)
  res.aov<-summary(res.aov)
  F.outcome<-res.aov[[1]]$`F value`[1]
  res<-list(outcome=F.outcome,gene=F.gene)
  return(res)
}

region_lambda<-function (lambda1 =18,lambda2=0,iteration = 10, Y, G,X,center,w,K){
  lambda.vector <- c(lambda1, lambda2)
  num.vector <- rep(-1, length(lambda.vector))
  R2.gene.vector<- rep(-1, length(lambda.vector))
  R2.outcome.vector<-rep(-1, length(lambda.vector))

  for (i in 1:length(num.vector)) {
    mod1 <- ogclust_con_Yujia(x=X,G=t(G),y=Y,c_center=center,lambda=lambda.vector[i],v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
    select.feature<-as.numeric(apply(mod1$result_list$mu,1,function(x){length(unique(x))})!=1)
    num.vector[i]<-sum(select.feature)

    cluster.train<-apply(mod1$result_list$z,1,which.max)
    mod.train<-Rsquare(cluster =cluster.train,Y = Y,X = X,G = t(G))
    R2.gene.vector[i]<-mean(mod.train$gene[which(select.feature==1)])
    R2.outcome.vector[i]<-mod.train$outcome
  }
  num.vector_old <- num.vector
  R2.gene.vector_old<-R2.gene.vector
  R2.outcome.vector_old<-R2.outcome.vector

  lambda <- (lambda.vector[1]+lambda.vector[2])/2
  iter<-2
  while (iter <= iteration) {
    #print(iter)
    lambda.vector <- sort(c(lambda.vector, lambda),decreasing = T)
    num.vector <- rep(-1, length(lambda.vector))
    R2.gene.vector<- rep(-1, length(lambda.vector))
    R2.outcome.vector<-rep(-1, length(lambda.vector))
    num.vector[-which(lambda.vector == lambda)] <- num.vector_old
    R2.gene.vector[-which(lambda.vector == lambda)] <- R2.gene.vector_old
    R2.outcome.vector[-which(lambda.vector == lambda)] <- R2.outcome.vector_old

    mod_new<-ogclust_con_Yujia(x=X,G=t(G),y=Y,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
    select.feature<-as.numeric(apply(mod_new$result_list$mu,1,function(x){length(unique(x))})!=1)
    num.vector[which(lambda.vector == lambda)] <-sum(select.feature)

    cluster.train<-apply(mod_new$result_list$z,1,which.max)
    mod.train<-Rsquare(cluster =cluster.train,Y = Y,X = X,G = t(G))
    R2.gene.vector[which(lambda.vector == lambda)] <- mean(mod.train$gene[which(select.feature==1)])
    R2.outcome.vector[which(lambda.vector == lambda)] <- mod.train$outcome

    d <- rep(-1, (length(lambda.vector) - 1))
    for (i in 1:length(d)) {
      if (num.vector[i] == num.vector[i + 1]) {
        d[i] <- 0
      }else {
        d[i] <- abs(log2(num.vector[i + 1]/num.vector[i]))
      }
    }
    index <- which.max(d)
    lambda <- (lambda.vector[index]+lambda.vector[index + 1])/2
    #print(lambda)
    iter <- iter + 1
    num.vector_old <- num.vector
    R2.gene.vector_old<-R2.gene.vector
    R2.outcome.vector_old<-R2.outcome.vector
  }
  lambda.vector <- sort(c(lambda.vector, lambda),decreasing =T)
  num.vector <- rep(-1, length(lambda.vector))
  num.vector[-which(lambda.vector == lambda)] <- num.vector_old
  R2.gene.vector<- rep(-1, length(lambda.vector))
  R2.outcome.vector<-rep(-1, length(lambda.vector))
  R2.gene.vector[-which(lambda.vector == lambda)] <- R2.gene.vector_old
  R2.outcome.vector[-which(lambda.vector == lambda)] <- R2.outcome.vector_old

  mod_new<-ogclust_con_Yujia(x=X,G=t(G),y=Y,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
  select.feature<-as.numeric(apply(mod_new$result_list$mu,1,function(x){length(unique(x))})!=1)
  num.vector[which(lambda.vector == lambda)] <-sum(select.feature)

  cluster.train<-apply(mod_new$result_list$z,1,which.max)
  mod.train<-Rsquare(cluster =cluster.train,Y = Y,X = X,G = t(G))
  R2.gene.vector[which(lambda.vector == lambda)] <- mean(mod.train$gene[which(select.feature==1)])
  R2.outcome.vector[which(lambda.vector == lambda)] <- mod.train$outcome

  res<-list(lambda=lambda.vector,num=num.vector,R2.gene=R2.gene.vector,R2.outcome=R2.outcome.vector)
  return(res)
}

mult_density1<-function(x,mu,sigma){

  lik<-dnorm(x,mean=mu,sd=sqrt(sigma),log=TRUE)
  return(lik)
}

predict.ogclust.test<-function(mod,K,D.test,D.train,X1,p=p){
  mu_est<-mod$result_list$mu
  pi_est<-mod$result_list$pi
  sigma_est<-mod$result_list$sigma
  coef_est<-mod$result_list$int_coef
  if(ncol(D.test)>1){
    D.test<-t(D.test)
  }
  #p<-ncol(D.test)
  #D.test<-t(D.test)
  #D.test<-as.matrix(D.test)
  mult_pdf<-matrix(,nrow=ncol(D.test),ncol=K)
  #------------------------------each sample in each cluster (sum of log likelihood)
  for(j1 in 1:K){
    temp<-mult_density1(as.numeric(D.test),mu=rep(mu_est[,j1],times=ncol(D.test)),sigma=rep(sigma_est,times=ncol(D.test)))
    #assume genes are uncorrelated, i.e. independently sample from normal distribution
    temp_matrix<-matrix(temp,nrow=p,ncol=ncol(D.test),byrow=FALSE)
    mult_pdf[,j1]<-t(temp_matrix)%*%rep(1,p) #the gene expression probability from cluster j1 for each sample
  }
  #---------------------get the predicted probability for each sample
  max_pdf<-apply(mult_pdf,1,max) #just based on the gene pribability?
  max_pdf_matrix<-matrix(rep(max_pdf,times=K),ncol=K)
  mult_pdf1<-mult_pdf-max_pdf_matrix
  mult_pdf2<-exp(mult_pdf1)*matrix(rep(pi_est,times=ncol(D.test)),byrow=T,ncol=K)
  sum_mult_pdf2<-mult_pdf2%*%rep(1,K)
  z_est<-mult_pdf2/matrix(rep(sum_mult_pdf2,times=K),ncol=K)
  beta_vector<-coef_est[1:K] #beta_0, i.e. intercepts for differnet clusters
  if(nrow(z_est)==1){
    y.pred<-t(as.matrix(X1))%*%as.matrix(coef_est[(K+1):length(coef_est)])+ z_est%*%beta_vector
  }else{
    y.pred<-as.matrix(X1)%*%as.matrix(coef_est[(K+1):length(coef_est)])+ z_est%*%beta_vector
  }

  cluster<-apply(mod$result_list$z,1,which.max) #the cluster for train data
  G.mean<-matrix(NA,nrow=ncol(D.train),ncol=K)
  for(i in 1:K){
    G.mean[,i]<-as.numeric(apply(D.train[which(cluster==i),],2,mean)) #calculate the center
  }

  index<-which(apply(mod$result_list$mu,1,function(x){length(unique(x))})!=1) #select genes
  G.pred<-as.matrix(G.mean)%*%as.matrix(t(z_est)) #predict the gene expression level based on the assigned clusters for test samples
  if(nrow(z_est)==1){ #if only one test sample
    clus<-which.max(z_est)
  }else{
    clus<-apply(z_est,1,which.max)
  }
  #clus<-which.max(z_est)
  res<-list(Y=y.pred,G=G.pred,clus=clus,index=index)
  #y.pred the predicted outcome Y, and G.pred is the predicted gene expression, index is the selected genes of the input trained mod.
  return(res)
}

predict_test<-function(mod,K,D.test,X1,p=p){
  mult_density1<-function(x,mu,sigma){
    lik<-dnorm(x,mean=mu,sd=sqrt(sigma),log=TRUE)
    return(lik)
  }
  mu_est<-mod$result_list$mu
  pi_est<-mod$result_list$pi
  sigma_est<-mod$result_list$sigma
  coef_est<-mod$result_list$int_coef
  if(ncol(D.test)>1){
    D.test<-t(D.test)
  }
  #p<-ncol(D.test)
  #D.test<-t(D.test)
  #D.test<-as.matrix(D.test)
  mult_pdf<-matrix(,nrow=ncol(D.test),ncol=K)
  #------------------------------each sample in each cluster (sum of log likelihood)
  for(j1 in 1:K){
    temp<-mult_density1(as.numeric(D.test),mu=rep(mu_est[,j1],times=ncol(D.test)),sigma=rep(sigma_est,times=ncol(D.test)))
    temp_matrix<-matrix(temp,nrow=p,ncol=ncol(D.test),byrow=FALSE)
    mult_pdf[,j1]<-t(temp_matrix)%*%rep(1,p)
  }
  #---------------------get the predicted probability for each sample
  max_pdf<-apply(mult_pdf,1,max)
  max_pdf_matrix<-matrix(rep(max_pdf,times=K),ncol=K)
  mult_pdf1<-mult_pdf-max_pdf_matrix
  mult_pdf2<-exp(mult_pdf1)*matrix(rep(pi_est,times=ncol(D.test)),byrow=T,ncol=K)
  sum_mult_pdf2<-mult_pdf2%*%rep(1,K)
  z_est<-mult_pdf2/matrix(rep(sum_mult_pdf2,times=K),ncol=K)
  beta_vector<-coef_est[1:K]

  if(nrow(z_est)==1){
    y.pred<-t(as.matrix(X1))%*%as.matrix(coef_est[(K+1):length(coef_est)])+ z_est%*%beta_vector
  }else{
    y.pred<-as.matrix(X1)%*%as.matrix(coef_est[(K+1):length(coef_est)])+ z_est%*%beta_vector
  }
  cluster<-apply(z_est,1,which.max)
  intercept<-rep(NA,nrow(X1))
  for(i in 1:K){
    intercept[which(cluster==i)]<-beta_vector[i]
  }
  y.pred.hard<-as.matrix(X1)%*%as.matrix(coef_est[(K+1):length(coef_est)])+intercept

  index<-which(apply(mod$result_list$mu,1,function(x){length(unique(x))})!=1)
  # G.pred<-as.matrix(G.mean)%*%as.matrix(t(z_est))
  if(nrow(z_est)==1){
    clus<-which.max(z_est)
  }else{
    clus<-apply(z_est,1,which.max)
  }
  #clus<-which.max(z_est)
  y.pred.hard
  #res<-list(Y=y.pred,G=G.pred,clus=clus,index=index)
  res<-list(Y=y.pred,clus=clus,index=index,Y.hard=y.pred.hard,z_est=z_est)
  return(res)
}

Rsquare<-function(cluster,Y,G,X){
  #cluster<-apply(mod$result_list$z,1,which.max)
  R2.gene<-c()
  for(i in 1:nrow(G)){
    data1<-data.frame(gene=G[i,],cluster=as.factor(cluster))
    res.lm <- lm(gene ~ cluster, data = data1)
    RSS<-sum(res.lm$residuals^2)
    res.lm <- lm(gene ~ 1, data = data1)
    TSS<-sum(res.lm$residuals^2)
    R2.gene[i]<-1-RSS/TSS
  }
  data1<-data.frame(Y=Y,cluster=as.factor(cluster),X=X)
  res.lm <- lm(Y~ cluster+X, data = data1)
  RSS<-sum(res.lm$residuals^2)
  res.lm <- lm(Y ~ 1+X, data = data1)
  TSS<-sum(res.lm$residuals^2)
  R2.outcome<-1-RSS/TSS
  res<-list(outcome=R2.outcome,gene=R2.gene)
  return(res)
}

# library(Brobdingnag)# for avoid the small number equals 0 in calculation
# calculate log likelihood for subject j and cluster k
mult_density<-function(x,mu,sigma){

  sum<-sum(dnorm(x,mean=mu,sd=sqrt(sigma),log=TRUE))
  return(sum)
}


em_mbc<-function(data,mu=0,sigma=1,lambda,c_center=NULL,v_int=NULL,pi_int=NULL,K=2,no_init=1,max_iter=200){
  result_list<-list()
  log_lik_vector<-c()
  unpen_lik_vector<-c()
  time<-c()
  #Column refers to samples and row refers to genes
  p<-dim(data)[1] # number of variables
  n<-dim(data)[2] # number of subjects

  for(init in 1:no_init){
    start<-Sys.time()
    #========== E-step: initialization===========#

    if(length(c_center)==0){
      # assume the means for each cluster follows the distribution N(0,1)
      mu_int<-matrix(rnorm(p*K,mu,sigma),nrow=p,ncol=K)
      #mu[which(mu==0)]<-100
    } else{
      mu_int=c_center
    }

    v_int<-ifelse(length(v_int)!=p,rep(sigma,p),v_int) # initial value of variance of each variables
    pi_int<-ifelse(length(pi_int)!=K,rep(1/K,K),pi_int) # initial value of pi
    z_int<-matrix(,nrow=n,ncol=K) # initial value of prob in cluster k for each subject

    mult_pdf<-matrix(,nrow=n,ncol=K)
    for(i in 1:n){
      for(j in 1:K){
        mult_pdf[i,j]<-as.numeric(mult_density(data[,i],mu=mu_int[,j],sigma=v_int))
      }
    }

    for(i in 1:n){
      d<-brob(mult_pdf[i,])*pi_int
      z_int[i,]<-as.numeric(d/sum(d))
    }

    pi<-pi_int
    v<-v_int
    mu<-mu_int
    z<-z_int

    #========= M step ==========#
    iter=1
    log_lik<-1 # initialize
    log_lik_up<-0 #initialize
    while(!is.na(abs(log_lik_up-log_lik)) & abs(log_lik_up-log_lik)>10^(-7) & iter <=max_iter){
      # log likelihood before updating
      log_lik<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf[x,])))-lambda*sum(abs(mu))
      #update pi
      pi_up<-apply(z,2,sum)/n
      #pi_up<-ifelse(pi_up==0,1*10^(-16),pi_up)
      #update v
      update_v<-function(x){
        sig<-sum(sapply(1:K,function(y) sum(z[,y]*(data[x,]-mu[x,y])^2)))/n
        return(sig)
      }
      v_up<-sapply(1:p,function(x) update_v(x))

      #update mu
      mu_tu<-matrix(,nrow=p,ncol=K)
      mu_up<-matrix(,nrow=p,ncol=K)

      for(i in 1:K){
        temp<-sapply(1:n,function(x) z[x,i]*data[,x])
        mu_tu[,i]<-apply(temp,1,sum)/sum(z[,i])
        mu_up[,i]<-ifelse(lambda<=abs(apply(temp,1,sum)/v),sign(mu_tu[,i])*(abs(mu_tu[,i])-(lambda/sum(z[,i]))*v),0)
      }


      #update z
      z_up<-matrix(,nrow=n,ncol=K)

      for(i in 1:n){
        for(j in 1:K){
          mult_pdf[i,j]<-as.numeric(mult_density(data[,i],mu=mu_up[,j],sigma=v_up))
        }
      }

      for(i in 1:n){
        d<-brob(mult_pdf[i,])*pi_up
        z_up[i,]<-as.numeric(d/sum(d))
      }

      # update parameters values
      pi<-pi_up
      v<-v_up
      mu<-mu_up
      z<-z_up
      #log likelihood after updating
      log_lik_up<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf[x,])))-lambda*sum(abs(mu))
      unpen_lik_up<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf[x,])))
      iter=iter+1
    }

    end<-Sys.time()
    time[init]<-end-start
    result_list[[init]] <- list('mu'=mu,'sigma'=v,'log_lik'=log_lik,'z'=z,'pi'=pi,'no_init'=init)
    log_lik_vector[init]<-log_lik_up
    unpen_lik_vector[init]<-unpen_lik_up
    if(is.na(log_lik_up)) print("NA ll found!!")
    #print(init)
  }
  max_lik<-unpen_lik_vector[which.max(log_lik_vector)]

  # the optimal initials is the one with the maximum log likelihood
  optimal_result<-result_list[[which.max(log_lik_vector)]]
  max_mu<-optimal_result$mu
  s<-apply(max_mu,1,function(x) sum(x==0))
  #d<-sum(s[which(s!=1)])
  q<-sum(s)
  #BIC
  if(lambda!=0) {
    BIC<--2*max_lik+log(n)*(K+p+K*p-1-q)
  } else if(lambda==0) {
    BIC<--2*max_lik+log(n)*(K+p+K*p-1)
  }
  selected.index=which(s!=K)
  list('optimal_result'=optimal_result,'result_list'=result_list,'time'=time,'BIC'=BIC,"selected.index"=selected.index)
}

predict.ogclust.test.select<-function(mod,K,D.test,D.train,X1,p=p,O.test,s_G,w){
  mu_est<-mod$result_list$mu
  pi_est<-mod$result_list$pi
  sigma_est<-mod$result_list$sigma
  sigma_int<-mod$result_list$int_sigma_coef
  coef_est<-mod$result_list$int_coef
  if(ncol(D.test)>1){
    D.test<-t(D.test)
  }
  #p<-ncol(D.test)
  #D.test<-t(D.test)
  #D.test<-as.matrix(D.test)
  mult_pdf<-matrix(,nrow=ncol(D.test),ncol=K)
  #------------------------------each sample in each cluster (sum of log likelihood)
  for(j1 in 1:K){
    temp<-mult_density1(as.numeric(D.test),mu=rep(mu_est[,j1],times=ncol(D.test)),sigma=rep(sigma_est,times=ncol(D.test)))
    temp_matrix<-matrix(temp,nrow=p,ncol=ncol(D.test),byrow=FALSE)
    mult_pdf[,j1]<-t(temp_matrix)%*%rep(1,p)
  }

  z_outcome<-matrix(,nrow=ncol(D.test),ncol=K)
  for(j in 1:K){
    z_outcome[,j]<-dnorm(O.test,mean=coef_est[j]+X1%*%coef_est[(K+1):length(coef_est)],sd=sigma_int,log=T)
    #z_outcome[,j]<-dnorm(y,mean=int_coef[j]+int_coef[4]*x.origin[,1]+int_coef[5]*x.origin[,2],sd=int_sigma_coef,log=T)
  }

  #z_outcome<-z_outcome/s_G
  mult_pdf.add<-(1-w)*mult_pdf+w*z_outcome

  max_pdf<-apply(mult_pdf.add,1,max)
  max_pdf_matrix<-matrix(rep(max_pdf,times=K),ncol=K)
  mult_pdf1<-mult_pdf.add-max_pdf_matrix
  mult_pdf2<-exp(mult_pdf1)*matrix(rep(pi_est,times=ncol(D.test)),byrow=T,ncol=K)
  sum_mult_pdf2<-mult_pdf2%*%rep(1,K)
  z_est<-mult_pdf2/matrix(rep(sum_mult_pdf2,times=K),ncol=K)
  #---------------------get the predicted probability for each sample
  # max_pdf<-apply(mult_pdf,1,max)
  # max_pdf_matrix<-matrix(rep(max_pdf,times=K),ncol=K)
  # mult_pdf1<-mult_pdf-max_pdf_matrix
  # mult_pdf2<-exp(mult_pdf1)*matrix(rep(pi_est,times=ncol(D.test)),byrow=T,ncol=K)
  # sum_mult_pdf2<-mult_pdf2%*%rep(1,K)
  # z_est<-mult_pdf2/matrix(rep(sum_mult_pdf2,times=K),ncol=K)
  beta_vector<-coef_est[1:K]

  if(nrow(z_est)==1){
    y.pred<-t(as.matrix(X1))%*%as.matrix(coef_est[(K+1):length(coef_est)])+ z_est%*%beta_vector
  }else{
    y.pred<-as.matrix(X1)%*%as.matrix(coef_est[(K+1):length(coef_est)])+ z_est%*%beta_vector
  }
  cluster<-apply(z_est,1,which.max)
  intercept<-rep(NA,nrow(X1))
  for(i in 1:K){
    intercept[which(cluster==i)]<-beta_vector[i]
  }
  y.pred.hard<-as.matrix(X1)%*%as.matrix(coef_est[(K+1):length(coef_est)])+intercept



  # cluster<-apply(mod$result_list$z,1,which.max)
  # G.mean<-matrix(NA,nrow=ncol(D.train),ncol=K)
  # for(i in 1:K){
  #   G.mean[,i]<-as.numeric(apply(D.train[which(cluster==i),],2,mean))
  # }
  #
  index<-which(apply(mod$result_list$mu,1,function(x){length(unique(x))})!=1)
  # G.pred<-as.matrix(G.mean)%*%as.matrix(t(z_est))
  if(nrow(z_est)==1){
    clus<-which.max(z_est)
  }else{
    clus<-apply(z_est,1,which.max)
  }
  #clus<-which.max(z_est)
  y.pred.hard
  #res<-list(Y=y.pred,G=G.pred,clus=clus,index=index)
  res<-list(Y=y.pred,clus=clus,index=index,Y.hard=y.pred.hard,z_est=z_est)
  return(res)
}

ogClust_Surv<-function(x,G,y,y.ind,c_center=NULL,lambda,v_int=NULL,pi_int=NULL,K=3,max_iter=200,w_outcome=0.5,w_G=0.5,z_int=NULL){
  x.origin<-x
  mult_density1<-function(x,mu,sigma){

    lik<-dnorm(x,mean=mu,sd=sqrt(sigma),log=TRUE)
    return(lik)
  }

  #Column refers to samples and row refers to genes

  p<-dim(G)[1] # number of variables
  n<-dim(G)[2] # number of subjects

  #========== E-step: initialization===========#
  mu_int=c_center


  if(length(v_int)!=p){
    v_int<-rep(1,p)
  }
  if(length(pi_int)!=p){
    pi_int<-rep(1/K,K)
  }
  pi<-pi_int
  v<-v_int
  mu<-mu_int


  if(is.null(z_int)){
    #-----------------------------------Set initial clustering assignment
    z_int<-matrix(,nrow=n,ncol=K) # initial value of prob in cluster k for each subject
    mult_pdf<-matrix(,nrow=n,ncol=K)
    #------------------------------each sample in each cluster (sum of log likelihood)
    for(j in 1:K){
      temp<-mult_density1(as.numeric(G),mu=rep(mu_int[,j],times=n),sigma=rep(v_int,times=n))
      temp_matrix<-matrix(temp,nrow=p,ncol=n,byrow=FALSE)
      mult_pdf[,j]<-t(temp_matrix)%*%rep(1,p)
    }

    #---------------------The initial clustering assignment is by Gene expression only
    max_pdf<-apply(mult_pdf,1,max)
    max_pdf_matrix<-matrix(rep(max_pdf,times=K),ncol=K)
    mult_pdf1<-mult_pdf-max_pdf_matrix
    mult_pdf2<-exp(mult_pdf1)*matrix(rep(pi_int,times=n),byrow=T,ncol=K)
    sum_mult_pdf2<-mult_pdf2%*%rep(1,K)
    z_int<-mult_pdf2/matrix(rep(sum_mult_pdf2,times=K),ncol=K)

    #initialize coefficient in outcome association
    #x<-cbind(1,x)
    #----------------------------------------
    x<-matrix(NA,nrow=nrow(x.origin)*K,ncol=ncol(x.origin))
    intercept<-matrix(0,nrow=n*K,ncol=K)
    for(i in 1:ncol(x.origin)){
      x[,i]<-rep(x.origin[,i],K)
    }
    for(i in 1:K){
      index<-((i-1)*n+1):(i*n)
      intercept[index,i]<-1
    }
    colnames(x)<-colnames(x.origin)
    colnames(intercept)<-paste("intercept",1:K,sep="")
    x<-cbind(intercept,x)
    W<-as.numeric(z_int)
    W[which(W == 0)] = 10^(-3)
    y1<-rep(y,K)
    y1.ind<-rep(y.ind,K)
    data1<-cbind(y1,y1.ind,x)
    data1<-as.data.frame(data1)
    mod<-survival::survreg(survival::Surv(y1, y1.ind) ~ -1+., data1,weights = W,robust=T,
                           dist="loglogistic")
    int_coef<-mod$coefficients
    int_sigma_coef<-mod$scale
    x<-x.origin
    #----------------------------------------

    #update the z matrix
    z_outcome<-matrix(,nrow=n,ncol=K)
    for(j in 1:K){
      Z<-log(y)
      #w<-(Z-int_coef[4]*x.origin[,1]-int_coef[5]*x.origin[,2]-int_coef[j])/int_sigma_coef
      w<-(Z-as.matrix(x.origin)%*%as.matrix(int_coef[(K+1):length(int_coef)])-int_coef[j])/int_sigma_coef
      lik_obs<-(1/int_sigma_coef)*(exp(w)/(1 + exp(w))^2)
      lik_unobs<-1/(1+exp(w))
      z_outcome[,j]<-(lik_obs^(y.ind))*(lik_unobs^(1-y.ind))

      #z_outcome[,j]<-dnorm(y,mean=int_coef[j]+int_coef[4]*x.origin[,1]+int_coef[5]*x.origin[,2],sd=int_sigma_coef,log=T)
    }
    z_outcome<-log(z_outcome)
    #mult_pdf<-mult_pdf*s_G
    #z_outcome<-z_outcome/s_G
    mult_pdf.add<-w_G*mult_pdf+w_outcome*z_outcome

    max_pdf<-apply(mult_pdf.add,1,max)
    max_pdf_matrix<-matrix(rep(max_pdf,times=K),ncol=K)
    mult_pdf1<-mult_pdf.add-max_pdf_matrix
    mult_pdf2<-exp(mult_pdf1)*matrix(rep(pi_int,times=n),byrow=T,ncol=K)
    sum_mult_pdf2<-mult_pdf2%*%rep(1,K)
    z<-mult_pdf2/matrix(rep(sum_mult_pdf2,times=K),ncol=K)

  }else{
    z<-z_int
  }

  pi_up<-apply(z,2,sum)/n
  v_up_matrix<-matrix(NA,ncol=K,nrow=p)
  for(j in 1:K){
    temp1<-(G-matrix(rep(mu[,j],times=n),ncol=n))^2
    v_up_matrix[,j]<-temp1%*%z[,j]/n
  }
  v_up<-v_up_matrix%*%rep(1,K)
  #update mu
  mu_tu<-matrix(,nrow=p,ncol=K)
  mu_up<-matrix(,nrow=p,ncol=K)
  for(i in 1:K){
    #temp<-sapply(1:n,function(x) z[x,i]*data[,x])
    G1<-G*matrix(rep(z[,i],times=p),ncol=n,byrow=T)
    mu_tu[,i]<-G1%*%rep(1,n)/sum(z[,i])
    #mu_tu[,i]<-apply(temp,1,sum)/sum(z[,i])
    #mu_up[,i]<-ifelse(lambda<=abs((G1%*%rep(1,n))/v),sign(mu_tu[,i])*(abs(mu_tu[,i])-(lambda/sum(z[,i]))*v),0)
    mu_up[,i]<-ifelse(lambda<=abs((G1%*%rep(1,n))/v),sign(mu_tu[,i])*(abs(mu_tu[,i])-(lambda/sum(z[,i]))*v),0)
  }

  mult_pdf<-matrix(,nrow=n,ncol=K)
  for(j in 1:K){
    temp<-mult_density1(as.numeric(G),mu=rep(mu_up[,j],times=n),sigma=rep(v_up,times=n))
    temp_matrix<-matrix(temp,nrow=p,ncol=n,byrow=FALSE)
    mult_pdf[,j]<-t(temp_matrix)%*%rep(1,p)
  }

  x<-matrix(NA,nrow=nrow(x.origin)*K,ncol=ncol(x.origin))
  intercept<-matrix(0,nrow=n*K,ncol=K)
  for(i in 1:ncol(x)){
    x[,i]<-rep(x.origin[,i],K)
  }
  for(i in 1:K){
    index<-((i-1)*n+1):(i*n)
    intercept[index,i]<-1
  }
  colnames(x)<-colnames(x.origin)
  colnames(intercept)<-paste("intercept",1:K,sep="")

  x<-cbind(intercept,x)
  W<-as.numeric(z)
  W[which(W == 0)] = 10^(-3)
  y1<-rep(y,K)
  y1.ind<-rep(y.ind,K)
  data1<-cbind(y1,y1.ind,x)
  # fit = survival::survreg(eval(parse(text = paste("Surv(Y,delta)~-1", paste(colnames(x), collapse = " + "), paste(paste0("mu", 1:K), collapse = " + "),
  #                                                 sep = " + "))), weights = weights, data = dt2, dist = dist, robust = TRUE)
  data1<-as.data.frame(data1)
  mod<-survival::survreg(survival::Surv(y1, y1.ind) ~ -1+., data1,weights = W,robust=T,
                         dist="loglogistic")
  int_coef<-mod$coefficients
  int_sigma_coef<-mod$scale
  x<-x.origin

  #update the z matrix
  z_outcome<-matrix(,nrow=n,ncol=K)
  for(j in 1:K){
    Z<-log(y)
    #w<-(Z-int_coef[4]*x.origin[,1]-int_coef[5]*x.origin[,2]-int_coef[j])/int_sigma_coef
    w<-(Z-as.matrix(x.origin)%*%as.matrix(int_coef[(K+1):length(int_coef)])-int_coef[j])/int_sigma_coef
    lik_obs<-(1/int_sigma_coef) * (exp(w)/(1 + exp(w))^2)
    lik_unobs<-1/(1+exp(w))
    z_outcome[,j]<-(lik_obs^(y.ind))*(lik_unobs^(1-y.ind))

    #z_outcome[,j]<-dnorm(y,mean=int_coef[j]+int_coef[4]*x.origin[,1]+int_coef[5]*x.origin[,2],sd=int_sigma_coef,log=T)
  }
  z_outcome<-log(z_outcome)
  #mult_pdf<-mult_pdf*s_G
  #z_outcome<-z_outcome/s_G
  mult_pdf.add<-w_G*mult_pdf+w_outcome*z_outcome

  max_pdf<-apply(mult_pdf.add,1,max)
  max_pdf_matrix<-matrix(rep(max_pdf,times=K),ncol=K)
  mult_pdf1<-mult_pdf.add-max_pdf_matrix
  mult_pdf2<-exp(mult_pdf1)*matrix(rep(pi_int,times=n),byrow=T,ncol=K)
  sum_mult_pdf2<-mult_pdf2%*%rep(1,K)
  z_up<-mult_pdf2/matrix(rep(sum_mult_pdf2,times=K),ncol=K)

  # update parameters values
  pi<-pi_up
  v<-v_up
  mu<-mu_up
  z<-z_up




  #table(apply(z,1,which.max))
  #========= M step ==========#
  iter=1
  log_lik<-1 # initialize
  log_lik_up<-0 #initialize

  while(abs(log_lik_up-log_lik)>10^(-7) & iter <=max_iter){
    #print(iter)
    log_lik<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf.add[x,])))-lambda*sum(abs(mu))
    pi_up<-apply(z,2,sum)/n
    v_up_matrix<-matrix(NA,ncol=K,nrow=p)
    for(j in 1:K){
      temp1<-(G-matrix(rep(mu[,j],times=n),ncol=n))^2
      v_up_matrix[,j]<-temp1%*%z[,j]/n
    }
    v_up<-v_up_matrix%*%rep(1,K)
    #update mu
    mu_tu<-matrix(,nrow=p,ncol=K)
    mu_up<-matrix(,nrow=p,ncol=K)
    for(i in 1:K){
      #temp<-sapply(1:n,function(x) z[x,i]*data[,x])
      G1<-G*matrix(rep(z[,i],times=p),ncol=n,byrow=T)
      mu_tu[,i]<-G1%*%rep(1,n)/sum(z[,i])
      #mu_tu[,i]<-apply(temp,1,sum)/sum(z[,i])
      #mu_up[,i]<-ifelse(lambda<=abs((G1%*%rep(1,n))/v),sign(mu_tu[,i])*(abs(mu_tu[,i])-(lambda/sum(z[,i]))*v),0)
      mu_up[,i]<-ifelse(lambda<=abs((G1%*%rep(1,n))/v),sign(mu_tu[,i])*(abs(mu_tu[,i])-(lambda/sum(z[,i]))*v),0)
    }

    mult_pdf<-matrix(,nrow=n,ncol=K)
    for(j in 1:K){
      temp<-mult_density1(as.numeric(G),mu=rep(mu_up[,j],times=n),sigma=rep(v_up,times=n))
      temp_matrix<-matrix(temp,nrow=p,ncol=n,byrow=FALSE)
      mult_pdf[,j]<-t(temp_matrix)%*%rep(1,p)
    }

    x<-matrix(NA,nrow=nrow(x.origin)*K,ncol=ncol(x.origin))
    intercept<-matrix(0,nrow=n*K,ncol=K)
    for(i in 1:ncol(x)){
      x[,i]<-rep(x.origin[,i],K)
    }
    for(i in 1:K){
      index<-((i-1)*n+1):(i*n)
      intercept[index,i]<-1
    }
    colnames(x)<-colnames(x.origin)
    colnames(intercept)<-paste("intercept",1:K,sep="")

    x<-cbind(intercept,x)
    W<-as.numeric(z)
    W[which(W == 0)] = 10^(-3)
    y1<-rep(y,K)
    y1.ind<-rep(y.ind,K)
    data1<-cbind(y1,y1.ind,x)
    # fit = survival::survreg(eval(parse(text = paste("Surv(Y,delta)~-1", paste(colnames(x), collapse = " + "), paste(paste0("mu", 1:K), collapse = " + "),
    #                                                 sep = " + "))), weights = weights, data = dt2, dist = dist, robust = TRUE)
    data1<-as.data.frame(data1)
    mod<-survival::survreg(survival::Surv(y1, y1.ind) ~ -1+., data1,weights = W,robust=T,
                           dist="loglogistic")
    int_coef<-mod$coefficients
    int_sigma_coef<-mod$scale
    x<-x.origin
    # intercept<-rep(NA,K)
    # for(i in 1:K){
    #   W<-diag(z[,i])
    #   intercept[i]<-(solve(t(x)%*%W%*%x)%*%t(x)%*%W%*%y)[1]
    #
    # }
    #update the z matrix
    z_outcome<-matrix(,nrow=n,ncol=K)
    for(j in 1:K){
      Z<-log(y)
      #w<-(Z-int_coef[4]*x.origin[,1]-int_coef[5]*x.origin[,2]-int_coef[j])/int_sigma_coef
      w<-(Z-as.matrix(x.origin)%*%as.matrix(int_coef[(K+1):length(int_coef)])-int_coef[j])/int_sigma_coef
      lik_obs<-(1/int_sigma_coef) * (exp(w)/(1 + exp(w))^2)
      lik_unobs<-1/(1+exp(w))
      z_outcome[,j]<-(lik_obs^(y.ind))*(lik_unobs^(1-y.ind))

      #z_outcome[,j]<-dnorm(y,mean=int_coef[j]+int_coef[4]*x.origin[,1]+int_coef[5]*x.origin[,2],sd=int_sigma_coef,log=T)
    }

    z_outcome<-log(z_outcome)
    #mult_pdf<-mult_pdf*s_G
    #z_outcome<-z_outcome/s_G
    mult_pdf.add<-w_G*mult_pdf+w_outcome*z_outcome

    max_pdf<-apply(mult_pdf.add,1,max)
    max_pdf_matrix<-matrix(rep(max_pdf,times=K),ncol=K)
    mult_pdf1<-mult_pdf.add-max_pdf_matrix
    mult_pdf2<-exp(mult_pdf1)*matrix(rep(pi_int,times=n),byrow=T,ncol=K)
    sum_mult_pdf2<-mult_pdf2%*%rep(1,K)
    z_up<-mult_pdf2/matrix(rep(sum_mult_pdf2,times=K),ncol=K)

    # update parameters values
    pi<-pi_up
    v<-v_up
    mu<-mu_up
    z<-z_up

    # if(sum(pi==0)!=0){
    #   pi[which(pi==0)]<-10^(-100)
    # }
    log_lik_up<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf.add[x,])))-lambda*sum(abs(mu))
    unpen_lik_up<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf.add[x,])))
    log_lik_uncomplete<-0
    max_mult<-apply(mult_pdf.add,1,max)
    for(i in 1:nrow(mult_pdf)){
      log_lik_uncomplete<-log_lik_uncomplete+sum(pi*exp(mult_pdf.add[i,]-max_mult[i]))+max_mult[i]
      #log_lik<-log_lik+log(sum(pi*exp(mult_pdf[i,]-max_mult[i])))+max_mult[i]
    }
    iter=iter+1
  }
  result_list <- list('mu'=mu,'sigma'=v,'z'=z,'pi'=pi,
                      'int_coef'=int_coef,'int_sigma_coef'=int_sigma_coef)
  log_lik<-log_lik_up
  unpen_lik<-unpen_lik_up
  max_lik<-log_lik_uncomplete
  max_mu<-result_list$mu
  s<-apply(max_mu,1,function(x) sum(x==0))
  q<-sum(s)
  #BIC
  if(lambda!=0) {
    BIC<--2*max_lik+log(n)*(K+p+K*p-q)
  } else if(lambda==0) {
    BIC<--2*max_lik+log(n)*(K+p+K*p)
  }

  res<-list('result_list'=result_list,'BIC'=BIC,'lik'=max_lik)
  return(res)


}

region.lambda.surv<-function (lambda1 =18,lambda2=0,iteration = 10, Y, G,X,center,w,s_G=s_G,K=K,delta.train){
  lambda.vector <- c(lambda1, lambda2)
  num.vector <- rep(-1, length(lambda.vector))
  for (i in 1:length(num.vector)) {
    mod1 <- ogClust_Surv(x=X,G=t(G),y=Y,y.ind = delta.train,c_center=center,lambda=lambda.vector[i],v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
    num.vector[i] <- sum(apply(mod1$result_list$mu,1,function(x){length(unique(x))})!=1)
  }
  num.vector_old <- num.vector
  lambda <- (lambda.vector[1]+lambda.vector[2])/2
  iter<-2
  while (iter <= iteration) {
    #print(iter)
    lambda.vector <- sort(c(lambda.vector, lambda),decreasing = T)
    num.vector <- rep(-1, length(lambda.vector))
    num.vector[-which(lambda.vector == lambda)] <- num.vector_old
    mod_new<-ogClust_Surv(x=X,G=t(G),y=Y,y.ind = delta.train,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)

    num.vector[which(lambda.vector == lambda)] <- sum(apply(mod_new$result_list$mu,1,function(x){length(unique(x))})!=1)
    d <- rep(-1, (length(lambda.vector) - 1))
    for (i in 1:length(d)) {
      if (num.vector[i] == num.vector[i + 1]) {
        d[i] <- 0
      }else {
        d[i] <- abs(log2(num.vector[i + 1]/num.vector[i]))
      }
    }
    index <- which.max(d)
    lambda <- (lambda.vector[index]+lambda.vector[index + 1])/2
    iter <- iter + 1
    num.vector_old <- num.vector
  }
  lambda.vector <- sort(c(lambda.vector, lambda),decreasing =T)
  num.vector <- rep(-1, length(lambda.vector))
  num.vector[-which(lambda.vector == lambda)] <- num.vector_old
  mod_new<-ogClust_Surv(x=X,G=t(G),y=Y,y.ind = delta.train,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
  num.vector[which(lambda.vector == lambda)] <- sum(apply(mod_new$result_list$mu,1,function(x){length(unique(x))})!=1)
  res<-list(lambda=lambda.vector,num=num.vector)
  return(res)
}


###GM method with penalized GLM#####
fit.ogClust <- function(n, K, np, NG, lambda, alpha, G, Y, X, theta_int, robust = "none", tau = 1.345) {
  stopifnot(robust %in% c("none", "huber", "median", "hubertf"))
  if (class(G) != "matrix")
    G = as.matrix(G)
  if (all(class(X) != "matrix"))
    X = as.matrix(X)
  theta_est = EM(theta_int, lambda = lambda, n = n, G = G, Y = Y, X = X, np = np, K = K, NG = NG, alpha = alpha, robust = robust, tau = tau)

  # estimated parameters
  beta_est = theta_est[1:np]
  gamma_est = theta_est[(np + 1):((K - 1) * (NG + 1) + np)]
  beta0_est = theta_est[((K - 1) * (NG + 1) + np + 1):((K - 1) * (NG + 1) + np + K)]
  sigma2_est = theta_est[((K - 1) * (NG + 1) + np + K + 1)]

  gamma_est_matrix = matrix(gamma_est, ncol = K - 1, byrow = T)
  gamma_est_matrix = cbind(gamma_est_matrix, 0)

  par_est<-list(beta0=beta0_est, beta=beta_est, sigma2=sigma2_est, gamma=gamma_est_matrix)
  G = cbind(1, G)
  pai_est = sapply(1:K, function(k) exp(G %*% gamma_est_matrix[, k, drop = F])/rowSums(exp(G %*% gamma_est_matrix)))
  f_est <- sapply(1:K, function(x) (1/sqrt(2 * pi * sigma2_est)) * exp(-(Y - beta0_est[x] - X %*% beta_est)^2/(2 * sigma2_est)))
  (ll = sum(log(diag(pai_est %*% t(f_est)))))

  # Calculate AIC BIC
  AIC = 2 * sum(theta_est != 0) - 2 * ll
  BIC = log(n) * sum(theta_est != 0) - 2 * ll
  # EBIC = BIC + 2 * (1 - 1/(2 * log(length(theta_est), base = n))) * log(choose(length(theta_est), sum(theta_est != 0)))

  # prosterior prob
  w_est = sapply(1:K, function(k) (pai_est[, k] * f_est[, k])/diag(pai_est %*% t(f_est)))
  cl.assign <- apply(w_est, 1, which.max)

  # calculate the expected value of Y and R2
  #Y_prd = apply(sapply(1:K, function(x) pai_est[, x] * (beta0_est[x] + X %*% beta_est)), 1, sum) #soft assignment
  Y_prd = beta0_est[cl.assign] + X %*% beta_est #hard assignment
  #R2 = 1 - sum((Y - Y_prd)^2)/sum((Y - mean(Y))^2)

  final.res <- list(par=par_est, ll = ll, AIC = AIC, BIC = BIC, lambda = lambda,
                    Y_prd=Y_prd, grp_assign=cl.assign)
  attr(final.res, "class") <- "ogClust"
  return(final.res)
}

EM <- function(theta, lambda, n, G, Y, X, np, K, NG, robust, alpha, tau = 1.345) {
  X = as.matrix(X)
  G = cbind(1, as.matrix(G))
  Y = Y

  l = 1
  repeat {
    # print(l)
    beta_old = theta[1:np]
    gamma_old = theta[(1 + np):((K - 1) * (NG + 1) + np)]
    miu_old = theta[((K - 1) * (NG + 1) + np + 1):((K - 1) * (NG + 1) + np + K)]
    sigma2_old = theta[((K - 1) * (NG + 1) + np + K + 1)]

    # ==E-STEP==#
    gamma_old_matrix = matrix(gamma_old, ncol = K - 1, byrow = T)
    gamma_old_matrix = cbind(gamma_old_matrix, 0)
    pai_old = sapply(1:K, function(k) exp(G %*% gamma_old_matrix[, k, drop = F])/rowSums(exp(G %*% gamma_old_matrix)))
    f_old <- sapply(1:K, function(x) (1/sqrt(2 * pi * sigma2_old)) * exp(-(Y - miu_old[x] - X %*% beta_old)^2/(2 * sigma2_old)))


    # calculate the expected value of Z
    w_old = sapply(1:K, function(k) (pai_old[, k] * f_old[, k])/diag(pai_old %*% t(f_old)))
    # ==M-STEP==#
    gamma_new_matrix = tryCatch({
      fit <- glmnet::glmnet(x = G[, -1], y = w_old, lambda = lambda, family = "multinomial", alpha = alpha, type.multinomial = "grouped")
      gamma_new_matrix = rbind(t(fit$a0), sapply(1:K, function(x) as.numeric(fit$beta[[x]])))
      gamma_new_matrix = sapply(1:K, function(x) gamma_new_matrix[, x] - gamma_new_matrix[, K])
    }, error = function(e) {
      return(gamma_old_matrix)
    })

    #---- non robust ----#
    if (robust == "none") {
      # update miu, alpha
      miu_new = sapply(1:K, function(k) {
        sum(w_old[, k] * (Y - X %*% beta_old), na.rm = T)/sum(w_old[, k], na.rm = T)
      })
      # update beta
      beta_new <- beta_old
      for (i in 1:np) {
        beta_new[i] <- sum(sapply(1:K, function(k) w_old[, k] * X[, i] * (Y - miu_new[k] - X[, -i, drop = F] %*% beta_new[-i])), na.rm = T)/sum(w_old *
                                                                                                                                                  X[, i]^2, na.rm = T)
      }
      # update sigma2
      sigma2_new = sum(sapply(1:K, function(k) w_old[, k] * (Y - miu_new[k] - X %*% beta_new)^2), na.rm = T)/n
    }

    #--- median truncated ---#
    if (robust == "median") {
      e = sapply(1:K, function(k) Y - miu_old[k] - X %*% beta_old)
      # update miu, alpha
      miu_new = sapply(1:K, function(k) sum((w_old[, k] * (Y - X %*% beta_old))[abs(e[, k]) <= median(abs(e[, k]))], na.rm = T)/sum(w_old[abs(e[, k]) <=
                                                                                                                                            median(abs(e[, k])), k], na.rm = T))
      # update beta1 beta2
      beta_new <- beta_old
      for (i in 1:np) {
        beta_new[i] = sum(unlist(sapply(1:K, function(k) (w_old[, k] * X[, i] * (Y - miu_new[k] - X[, -i, drop = F] %*% beta_new[-i]))[abs(e[, k]) <=
                                                                                                                                         median(abs(e[, k]))])), na.rm = T)/sum(unlist(sapply(1:K, function(k) (w_old[, k] * X[, i]^2)[abs(e[, k]) <= median(abs(e[, k]))])), na.rm = T)
      }
      # update sigma2
      sigma2_new = sum(unlist(sapply(1:K, function(k) (w_old[, k] * (Y - miu_new[k] - X %*% beta_new)^2)[abs(e[, k]) <= median(abs(e[, k]))])), na.rm = T)/sum(unlist(sapply(1:K,
                                                                                                                                                                             function(k) w_old[abs(e[, k]) <= median(abs(e[, k])), k])), na.rm = T)

    }

    #--- huber ---#
    if (robust == "huber") {
      e = sapply(1:K, function(k) Y - miu_old[k] - X %*% beta_old)
      # update miu, alpha
      miu_new = sapply(1:K, function(k) sum(c((w_old[, k] * (Y - X %*% beta_old))[abs(e[, k]) <= tau], (w_old[, k] * tau * sign(e[, k]))[abs(e[, k]) >
                                                                                                                                           tau]), na.rm = T)/sum(w_old[abs(e[, k]) <= tau, k], na.rm = T))
      # update beta
      beta_new = beta_old
      for (i in 1:np) {
        beta_new[i] = sum(sapply(1:K, function(k) c((w_old[, k] * X[, i] * (Y - miu_new[k] - X[, -i, drop = F] %*% beta_new[-i]))[abs(e[, k]) <= tau],
                                                    (w_old[, k] * tau * sign(e[, k]))[abs(e[, k]) > tau])), na.rm = T)/sum(sapply(1:K, function(k) sum((w_old[, k] * X[, i]^2)[abs(e[, k]) <= tau],
                                                                                                                                                       na.rm = T)), na.rm = T)
      }
      # update sigma2
      sigma2_new = sum(sapply(1:K, function(k) c((w_old[, k] * (Y - miu_new[k] - X %*% beta_new)^2)[abs(e[, k]) <= tau], (w_old[, k] * (2 * tau * abs(Y -
                                                                                                                                                        miu_new[k] - X %*% beta_new) - tau^2))[abs(e[, k]) > tau])), na.rm = T)/sum(sapply(1:K, function(k) sum(w_old[abs(e[, k]) <= tau, ], na.rm = T)),
                                                                                                                                                                                                                                    na.rm = T)
    }

    #--- adaptive huber ---#
    if (robust == "hubertf") {
      beta_new = beta_old
      miu_new = miu_old
      grp.id <- t(1 * apply(w_old, 1, function(x) x == max(x)))[, -1]
      X_tf <- cbind(grp.id, X)

      listHuber = adaHuber::huberReg(X_tf, Y)
      mm <- diag(K)
      mm[1, ] <- 1
      miu_new <- listHuber$theta[1:K] %*% mm
      beta_new <- listHuber$theta[-c(1:K)]
      tau = listHuber$tauCoef
      e = sapply(1:K, function(k) Y - miu_new[k] - X %*% beta_new)
      sigma2_new = sum(sapply(1:K, function(k) c((w_old[, k] * (Y - miu_new[k] - X %*% beta_new)^2)[abs(e[, k]) <= tau], (w_old[, k] * (2 * tau * abs(Y -
                                                                                                                                                        miu_new[k] - X %*% beta_new) - tau^2))[abs(e[, k]) > tau])), na.rm = T)/sum(sapply(1:K, function(k) sum(w_old[abs(e[, k]) <= tau, ], na.rm = T)),
                                                                                                                                                                                                                                    na.rm = T)
    }

    theta_new = c(beta_new, as.numeric(t(gamma_new_matrix[, -K])), miu_new, sigma2_new)
    dis = sqrt(sum((theta_new - theta)^2, na.rm = T))
    theta = theta_new
    l = l + 1
    if (dis < 1 * 10^(-6) | l > 200) {
      break
    }
  }
  return(theta)
}

fit.ogClust.surv <- function(n, K, np, NG, lambda, alpha, G, Y, X, delta, theta_int, dist = "loglogistic") {
  G = cbind(1, as.matrix(G))
  X = as.matrix(X)
  theta_est = EM.surv(theta_int, lambda = lambda, n = n, G = G, Y = Y, X = X, delta = delta, np = np, K = K, NG = NG, alpha = alpha, dist = dist)$theta

  # estimated parameters
  beta_est = theta_est[1:np]
  gamma_est = theta_est[(np + 1):((K - 1) * (NG + 1) + np)]
  beta0_est = theta_est[((K - 1) * (NG + 1) + np + 1):((K - 1) * (NG + 1) + np + K)]
  sigma2_est = theta_est[((K - 1) * (NG + 1) + np + K + 1)]

  gamma_est_matrix = matrix(gamma_est, ncol = K - 1, byrow = T)
  gamma_est_matrix = cbind(gamma_est_matrix, 0)

  par_est<-list(beta0=beta0_est, beta=beta_est, sigma2=sigma2_est,gamma=gamma_est_matrix)
  pai_est = sapply(1:K, function(k) exp(G %*% gamma_est_matrix[, k, drop = F])/rowSums(exp(G %*% gamma_est_matrix)))
  f_est = sapply(1:K, function(x) f_calc(Y1 = Y, X1 = X, beta = beta_est, mu = beta0_est[x], sigma2 = sigma2_est, delta = delta))
  f_est = t(apply(f_est, 1, function(x) x/sum(x)))
  idx = which.max(beta0_est)
  test = apply(f_est, 1, sum)

  f_est[which(test < 10^-3 | is.na(test)), idx] = 1
  f_est[is.na(f_est)] = 0

  (ll = sum(log(diag(pai_est %*% t(f_est)))))

  # Calculate AIC BIC
  AIC = 2 * sum(theta_est != 0) - 2 * ll
  BIC = log(n) * sum(theta_est != 0) - 2 * ll

  # prosterior prob
  w_est = sapply(1:K, function(k) (pai_est[, k] * f_est[, k])/diag(pai_est %*% t(f_est)))
  cl.assign <- apply(w_est, 1, which.max)

  # calculate the expected value of Y and R2
  #Y_prd = apply(sapply(1:K, function(x) pai_est[, x] * (beta0_est[x] + X %*% beta_est)), 1, sum) #soft assignment
  Y_prd = beta0_est[cl.assign] + X %*% beta_est #hard assignment
  #R2 = 1 - sum((Y - Y_prd)^2)/sum((Y - mean(Y))^2)

  final.res <- list(par=par_est, ll = ll, AIC = AIC, BIC = BIC, lambda = lambda,
                    Y_prd=Y_prd, grp_assign=cl.assign)
  attr(final.res, "class") <- "ogClust"
  return(final.res)
}


EM.surv <- function(theta, lambda, n, G, Y, X, delta, np, K, NG, alpha = 0.5, dist) {
  #-----------------------------------------------#
  l = 1
  repeat {
    beta_old = theta[1:np]
    gamma_old = theta[(1 + np):((K - 1) * (NG + 1) + np)]

    miu_old = theta[((K - 1) * (NG + 1) + np + 1):((K - 1) * (NG + 1) + np + K)]

    sigma2_old = theta[((K - 1) * (NG + 1) + np + K + 1):length(theta)]


    # ==E-STEP==#
    gamma_old_matrix = matrix(gamma_old, ncol = K - 1, byrow = T)
    gamma_old_matrix = cbind(gamma_old_matrix, 0)
    #pai_old = sapply(1:K, function(k) exp(G %*% gamma_old_matrix[, k, drop = F])/rowSums(exp(G %*% gamma_old_matrix)))
    pai_old = sapply(1:K, function(k) 1/rowSums(exp(sweep(G %*% gamma_old_matrix, 1, G %*% gamma_old_matrix[, k, drop = F]))))

    f_old = sapply(1:K, function(x) f_calc(Y1 = Y, X1 = X, beta = beta_old, mu = miu_old[x], sigma2 = sigma2_old, delta = delta))
    f_old = t(apply(f_old, 1, function(x) x/sum(x)))
    idx = which.max(miu_old)
    test = apply(f_old, 1, sum)

    f_old[which(test < 10^-3 | is.na(test)), idx] = 1
    f_old[is.na(f_old)] = 0
    # calculate the expected value of Z
    w_old = sapply(1:K, function(k) (pai_old[, k] * f_old[, k])/diag(pai_old %*% t(f_old)))

    # ==M-STEP==#
    #gamma_new_matrix = tryCatch({

    if(any(apply(w_old,2, function(x) sum(abs(x)))<1e-05)){
      w_old[,apply(w_old,2, function(x) sum(abs(x)))<1e-05]<-1e-05/n
    }
    fit <- glmnet::glmnet(x = G[, -1], y = w_old, lambda = lambda, family = "multinomial", alpha = alpha, type.multinomial = "grouped")
    gamma_new_matrix = rbind(t(fit$a0), sapply(1:K, function(x) as.numeric(fit$beta[[x]])))
    gamma_new_matrix = sapply(1:K, function(x) gamma_new_matrix[, x] - gamma_new_matrix[, K])
    #}, error = function(e) {
    #    return(gamma_old_matrix)
    #})

    if (is.null(colnames(X))) {
      colnames(X) = paste0("X", 1:np)
    }
    dt0 = data.frame(Y, delta, X)
    dt1 = dt0
    for (k in 1:(K - 1)) dt1 = rbind.data.frame(dt1, dt0)
    dt2 = dt1
    for (k in 1:K) dt2 = cbind.data.frame(dt2, rep(diag(K)[k, ], each = dim(dt0)[1]))
    colnames(dt2)[(ncol(dt2) - K + 1):ncol(dt2)] <- paste0("mu", 1:K)
    weights = vector()
    for (k in 1:K) {
      weights = c(weights, w_old[, k])
    }
    weights[which(weights == 0)] = 10^-3

    fit = survival::survreg(eval(parse(text = paste("Surv(Y,delta)~-1", paste(colnames(X), collapse = " + "), paste(paste0("mu", 1:K), collapse = " + "),
                                                    sep = " + "))), weights = weights, data = dt2, dist = dist, robust = FALSE)
    miu_new = fit$coefficients[(np + 1):(np + K)]
    sigma2_new = fit$scale
    beta_new = fit$coefficients[1:np]

    theta_new = c(beta_new, as.numeric(t(gamma_new_matrix[, -K])), miu_new, sigma2_new)
    dis = sqrt(sum((theta_new - theta)^2, na.rm = T))
    theta = theta_new
    l = l + 1
    if (dis < 1 * 10^(-6) | l > 500) {
      break
    }
  }
  return(list(theta = theta, w = w_old, l = l))
}

f_calc = function(Y1, X1, beta, mu, sigma2, delta, K) {
  Z = log(Y1)
  X0 = X1 %*% beta
  mu0 = mu
  W = (Z - X0 - mu0)/sigma2
  pdf_calc = function(delta, W, ind) {
    return((1/sigma2) * (exp(W[ind])/(1 + exp(W[ind]))^2)^(delta[ind]))
  }
  res = sapply(1:length(Y1), function(ind) (1/(1 + exp(W[ind])))^(1 - delta[ind]) * pdf_calc(delta, W, ind))
  return(res)
}

###GM method with penalized LDA#####
fit.ogClust.LDA <- function(K, lambda,G, Y, X, max_iter=200, tau = 1.345) {
  X = as.matrix(X)
  G = as.matrix(G)
  n=nrow(G)
  NG=ncol(G)
  np=ncol(X)
  EMout = EM.LDA(lambda = lambda, G = G, Y = Y, X = X, K = K, max_iter=max_iter,tau = tau)
  theta_est=EMout$theta
  # estimated parameters
  beta_est = theta_est[1:np]
  beta0_est = theta_est[(np + 1):(np + K)]
  sigma2_est = theta_est[(np + K + 1)]

  cluster=EMout$label
  out <- PenalizedLDA(G,cluster,xte=G,lambda=lambda,K=(K-1))
  PLDAlist=Classify(out$xproj,out$xteproj,cluster)
  pai_est=PLDAlist$pi_est
  pai_est[which(is.na(pai_est))]=1
  f_est <- sapply(1:K, function(x) (1/sqrt(2 * pi * sigma2_est)) * exp(-(Y - beta0_est[x] - X %*% beta_est)^2/(2 * sigma2_est)))

  # prosterior prob
  w_est = sapply(1:K, function(k) (pai_est[, k] * f_est[, k])/diag(pai_est %*% t(f_est)))
  cl.assign <- apply(w_est, 1, which.max)

  Y_prd = beta0_est[cl.assign] + X %*% beta_est

  final.res <- list(par=theta_est, label=cluster, lambda = lambda,
                    Y_prd=Y_prd, grp_assign=cl.assign)
  attr(final.res, "class") <- "ogClust"
  return(final.res)
}

EM.LDA <- function(lambda, G, Y, X, K, max_iter=200, tau = 1.345) {
  np=ncol(X)
  n=nrow(G)
  mod.kmeans<-kmeans(G,centers = K,nstart = 50)
  cluster<-mod.kmeans$cluster

  beta0_int=c()
  for (k in 1:K) {
    index1<-which(cluster==k)
    data1<-data.frame(y=Y[index1],X[index1,])
    mod1<-lm(y~.,data=data1)
    mod1<-summary(mod1)
    beta0_int=c(beta0_int,mod1$coefficients[1,1])
  }

  data<-data.frame(y=Y,X)
  mod<-lm(y~.,data=data)
  mod<-summary(mod)
  beta_int<-mod$coefficients[2:(np+1),1]

  sigma2_int<-1

  class_label=cluster
  out <- PenalizedLDA(G,class_label,xte=G,lambda=lambda,K=(K-1))
  PLDAlist=Classify(out$xproj,out$xteproj,class_label)

  theta= c(beta_int, beta0_int, sigma2_int,as.numeric(t(PLDAlist$mus)))
  l = 1
  repeat {
    #print(l)
    beta_old = theta[1:np]
    miu_old = theta[(np + 1):(np + K)]
    sigma2_old = theta[(np + K + 1)]
    # ==E-STEP==#
    f_old <- sapply(1:K, function(x) (1/sqrt(2 * pi * sigma2_old)) * exp(-(Y - miu_old[x] - X %*% beta_old)^2/(2 * sigma2_old)))
    pai_old=PLDAlist$pi_est #n by K
    pai_old[which(is.na(pai_old))]=1
    # calculate the expected value of Z
    w_old = sapply(1:K, function(k) (pai_old[, k] * f_old[, k])/diag(pai_old %*% t(f_old)))
    # ==M-STEP==#
    label_new=apply(w_old,1,which.max)
    if (is.list(label_new)) {
      empty.c=which(sapply(label_new,length)==0)
      label_new[empty.c]=sample(1:K,length(empty.c),replace =T)
      label_new=unlist(label_new)
    }
    #print(table(label_new))
    if (length(table(label_new))==K) {
      class_label=label_new
      out <- PenalizedLDA(G,class_label,xte=G,lambda=lambda,K=(K-1))
      PLDAlist=Classify(out$xproj,out$xteproj,class_label)
    }

    miu_new = sapply(1:K, function(k) {
      sum(w_old[, k] * (Y - X %*% beta_old), na.rm = T)/sum(w_old[, k], na.rm = T)
    })
    # update beta
    beta_new <- beta_old
    for (i in 1:np) {
      beta_new[i] <- sum(sapply(1:K, function(k) w_old[, k] * X[, i] * (Y - miu_new[k] - X[, -i, drop = F] %*% beta_new[-i])), na.rm = T)/sum(w_old *X[, i]^2, na.rm = T)
    }
    # update sigma2
    sigma2_new = sum(sapply(1:K, function(k) w_old[, k] * (Y - miu_new[k] - X %*% beta_new)^2), na.rm = T)/n


    theta_new = c(beta_new, miu_new, sigma2_new,as.numeric(t(PLDAlist$mus)))
    dis = sqrt(sum((theta_new - theta)^2, na.rm = T))
    theta = theta_new
    l = l + 1
    if (dis < 1 * 10^(-6) | l > max_iter) {
      break
    }
  }
  return(list(theta=theta,label=class_label))
}

# test=fit.ogClust(K=3,lambda=0.0333,G=X.GSE47460_1,Y=Y1,X=Z.GSE47460_1)
# cluster.train=test$grp_assign
# out.test <- PenalizedLDA(X.GSE47460_1,test$label,xte=X.GSE47460_2,lambda=test$lambda,K=(K-1)) #Time difference of 30.10021 secs
# cluster.test=out.test$ypred[,(K-1)]
# select.feature=apply(out.test$discrim,1,function(x) any(x!=0)) #about 472 with lambda=0.0333

#PLDAlist=Classify(out$xproj,out$xteproj,cluster)


fit.ogClust.LAD.surv <- function(K, lambda, G, Y, X, delta, dist = "loglogistic",max_iter=200) {
  X = as.matrix(X)
  G = as.matrix(G)
  n=nrow(G)
  NG=ncol(G)
  np=ncol(X)
  EMout = EM.LDA.surv(lambda = lambda, G = G, Y = Y, X = X, delta = delta, K = K, dist = dist,max_iter=max_iter)
  theta_est=EMout$theta
  # estimated parameters
  beta_est = theta_est[1:np]
  beta0_est = theta_est[(np + 1):(np + K)]
  sigma2_est = theta_est[(np + K + 1)]

  cluster=EMout$label
  out <- PenalizedLDA(G,cluster,xte=G,lambda=lambda,K=(K-1))
  PLDAlist=Classify(out$xproj,out$xteproj,cluster)
  pai_est=PLDAlist$pi_est
  pai_est[which(is.na(pai_est))]=1

  f_est = sapply(1:K, function(x) f_calc(Y1 = Y, X1 = X, beta = beta_est, mu = beta0_est[x], sigma2 = sigma2_est, delta = delta))
  f_est = t(apply(f_est, 1, function(x) x/sum(x)))
  idx = which.max(beta0_est)
  test = apply(f_est, 1, sum)

  f_est[which(test < 10^-3 | is.na(test)), idx] = 1
  f_est[is.na(f_est)] = 0

  # prosterior prob
  w_est = sapply(1:K, function(k) (pai_est[, k] * f_est[, k])/diag(pai_est %*% t(f_est)))
  cl.assign <- apply(w_est, 1, which.max)

  Y_prd = beta0_est[cl.assign] + X %*% beta_est #hard assignment

  final.res <- list(par=theta_est, label=cluster, lambda = lambda,
                    Y_prd=Y_prd, grp_assign=cl.assign)

  attr(final.res, "class") <- "ogClust"
  return(final.res)
}


EM.LDA.surv <- function(lambda, G, Y, X, delta, K, dist,max_iter) {

  np=ncol(X)
  n=nrow(G)
  mod.kmeans<-kmeans(G,centers = K,nstart = 50)
  cluster<-mod.kmeans$cluster

  beta0_int=c()
  for (k in 1:K) {
    index1<-which(cluster==k)
    data1<-data.frame(y=Y[index1],X[index1,])
    mod1<-lm(y~.,data=data1)
    mod1<-summary(mod1)
    beta0_int=c(beta0_int,mod1$coefficients[1,1])
  }

  data<-data.frame(y=Y,X)
  mod<-lm(y~.,data=data)
  mod<-summary(mod)
  beta_int<-mod$coefficients[2:(np+1),1]

  sigma2_int<-1

  class_label=cluster
  out <- PenalizedLDA(G,class_label,xte=G,lambda=lambda,K=(K-1))
  PLDAlist=Classify(out$xproj,out$xteproj,class_label)

  theta= c(beta_int, beta0_int, sigma2_int,as.numeric(t(PLDAlist$mus)))
  #-----------------------------------------------#
  l = 1
  repeat {
    beta_old = theta[1:np]
    miu_old = theta[(np + 1):(np + K)]
    sigma2_old = theta[(np + K + 1)]

    # ==E-STEP==#
    pai_old=PLDAlist$pi_est #n by K
    pai_old[which(is.na(pai_old))]=1

    f_old = sapply(1:K, function(x) f_calc(Y1 = Y, X1 = X, beta = beta_old, mu = miu_old[x], sigma2 = sigma2_old, delta = delta))
    f_old = t(apply(f_old, 1, function(x) x/sum(x)))
    idx = which.max(miu_old)
    test = apply(f_old, 1, sum)

    f_old[which(test < 10^-3 | is.na(test)), idx] = 1
    f_old[is.na(f_old)] = 0
    # calculate the expected value of Z
    w_old = sapply(1:K, function(k) (pai_old[, k] * f_old[, k])/diag(pai_old %*% t(f_old)))

    # ==M-STEP==#
    #gamma_new_matrix = tryCatch({

    if(any(apply(w_old,2, function(x) sum(abs(x)))<1e-05)){
      w_old[,apply(w_old,2, function(x) sum(abs(x)))<1e-05]<-1e-05/n
    }

    label_new=apply(w_old,1,which.max)
    if (is.list(label_new)) {
      empty.c=which(sapply(label_new,length)==0)
      label_new[empty.c]=sample(1:K,length(empty.c))
      label_new=unlist(label_new)
    }

    #print(table(label_new))
    if (length(table(label_new))==K) {
      class_label=label_new
      out <- PenalizedLDA(G,class_label,xte=G,lambda=lambda,K=(K-1))
      PLDAlist=Classify(out$xproj,out$xteproj,class_label)
    }

    if (is.null(colnames(X))) {
      colnames(X) = paste0("X", 1:np)
    }
    dt0 = data.frame(Y, delta, X)
    dt1 = dt0
    for (k in 1:(K - 1)) dt1 = rbind.data.frame(dt1, dt0)
    dt2 = dt1
    for (k in 1:K) dt2 = cbind.data.frame(dt2, rep(diag(K)[k, ], each = dim(dt0)[1]))
    colnames(dt2)[(ncol(dt2) - K + 1):ncol(dt2)] <- paste0("mu", 1:K)
    weights = vector()
    for (k in 1:K) {
      weights = c(weights, w_old[, k])
    }
    weights[which(weights == 0)] = 10^-3

    fit = survival::survreg(eval(parse(text = paste("Surv(Y,delta)~-1", paste(colnames(X), collapse = " + "), paste(paste0("mu", 1:K), collapse = " + "),
                                                    sep = " + "))), weights = weights, data = dt2, dist = dist, robust = FALSE)
    miu_new = fit$coefficients[(np + 1):(np + K)]
    sigma2_new = fit$scale
    beta_new = fit$coefficients[1:np]

    theta_new = c(beta_new, miu_new, sigma2_new,as.numeric(t(PLDAlist$mus)))
    dis = sqrt(sum((theta_new - theta)^2, na.rm = T))
    theta = theta_new
    l = l + 1
    if (dis < 1 * 10^(-6) | l > max_iter) {
      break
    }
  }
  return(list(theta=theta,label=class_label))
}

Classify <- function(xtr,xte,ytr,equalpriors=TRUE){ # I introduced unequal priors on 02/22/2010
  prior <- rep(1/length(unique(ytr)), length(unique(ytr)))
  if(!equalpriors){
    for(k in 1:length(unique(ytr))) prior[k] <- mean(ytr==k)
  }
  # classify test obs to nearest training centroid.
  if(is.matrix(xtr) && ncol(xtr)>1){
    mus <- matrix(0, nrow=ncol(xtr), ncol=length(unique(ytr)))
    for(k in 1:length(unique(ytr))){
      mus[,k] <- apply(xtr[ytr==k,], 2, mean)
    }
  } else {
    mus <- matrix(NA, nrow=1, ncol=length(unique(ytr)))
    for(k in 1:length(unique(ytr))) mus[1,k] <- mean(xtr[ytr==k])
  }
  negdists <- diag.disc(t(xte), mus, prior)
  class_label=apply(negdists,1,which.max)
  pi_est=exp(negdists)/apply(exp(negdists),1,sum)
  return(list(class_label=class_label,pi_est=pi_est,mus=mus))
}

wcsd <- function(vec, y){
  K <- length(unique(y))
  n <- length(vec)
  tots <- 0
  for(k in unique(y)){
    tots <- tots + sum((vec[y==k]-mean(vec[y==k]))^2)
  }
  return(sqrt(tots/n))
}

diag.disc <-function(x, centroids, prior) {
  dd <- t(x) %*% centroids
  dd0 <- (rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
  scale(dd, as.numeric(dd0), FALSE) # this is -.5*||x_i - mu_k||^2+log(pi_k)
}


#############Generate simulate data##########
Sim3<-function(n,beta1,beta2,q,q1,q2,q3,c1,var_g,mu,mu1,sigma_y){
  #simulate covariates from N(0,1)
  x1<-rnorm(n,0,1)
  x2<-rnorm(n,0,1)
  x<-cbind(x1,x2)
  #simulate y
  y<-rep(NA,n)
  y[1:(n/3)]<-rnorm(n/3,0+beta1*x1[1:(n/3)]+beta2*x2[1:(n/3)],sigma_y)
  y[(n/3+1):(2*n/3)]<-rnorm(n/3,c1+beta1*x1[(n/3+1):(2*n/3)]+beta2*x2[(n/3+1):(2*n/3)],sigma_y)
  y[(2*n/3+1):n]<-rnorm(n/3,2*c1+beta1*x1[(2*n/3+1):n]+beta2*x2[(2*n/3+1):n],sigma_y)
  #simulate x
  X1 = rbind(matrix(rnorm((n/3) * q), ncol = q),
             matrix(rnorm((n/3) * q), ncol = q),
             matrix(rnorm((n/3) * q), ncol = q))
  X2 = rbind(matrix(rnorm((n/3) * q1), ncol = q1),
             matrix(rnorm((n/3) * q1), ncol = q1),
             matrix(rnorm((n/3) * q1), ncol = q1))
  X3 = rbind(matrix(rnorm((n/3) * q2), ncol = q2),
             matrix(rnorm((n/3) * q2), ncol = q2),
             matrix(rnorm((n/3) * q2), ncol = q2))
  X.noise = rbind(matrix(rnorm((n/3) * q3), ncol = q3),
                  matrix(rnorm((n/3) * q3), ncol = q3),
                  matrix(rnorm((n/3) * q3), ncol = q3))
  X1[1:(n/3), 1:q] <- X1[1:(n/3), 1:q] - mu
  X1[(n/3+1):(2*n/3), 1:q] <- X1[(n/3+1):(2*n/3), 1:q]
  X1[(2*n/3+1):n, 1:q] <- X1[(2*n/3+1):n, 1:q] + mu
  index<-sample(1:n)
  X2[index[1:(n/3)], ] <- X2[index[1:(n/3)],]-mu1
  X2[index[(n/3+1):(2*n/3)], ] <- X2[index[(n/3+1):(2*n/3)], ]
  X2[index[(2*n/3+1):n], ] <- X2[index[(2*n/3+1):n], ]+mu1
  for(i in 1:ncol(X3)){
    X3[,i]<-X3[,i]+rnorm(1,0,var_g)*x1
  }
  X<-cbind(X1,X2,X3,X.noise)
  data<-list(x=x,y=y,G=X)
  return(data)
}

Sim_surv<-function(n,beta1,beta2,q,q1,q2,q3,c1,var_g,mu,mu1,sigma_y,censor){
  #simulatex
  x1<-rnorm(n,1,0.5)
  x2<-rnorm(n,1,0.5)
  x<-cbind(x1,x2)
  #simulate y
  y<-rep(NA,n)
  y[1:(n/3)]<-0+beta1*x1[1:(n/3)]+beta2*x2[1:(n/3)]+rlogis(n/3, location = 0, scale = 1)*sigma_y
  y[(n/3+1):(2*n/3)]<-c1+beta1*x1[(n/3+1):(2*n/3)]+beta2*x2[(n/3+1):(2*n/3)]+rlogis(n/3, location = 0, scale = 1)*sigma_y
  y[(2*n/3+1):n]<-2*c1+beta1*x1[(2*n/3+1):(n)]+beta2*x2[(2*n/3+1):(n)]+rlogis(n/3, location = 0, scale = 1)*sigma_y
  y<-exp(y)
  y.ind<-ifelse(y>censor,0,1)
  y<-ifelse(y>censor,censor,y)
  X1 = rbind(matrix(rnorm((n/3) * q), ncol = q),
             matrix(rnorm((n/3) * q), ncol = q),
             matrix(rnorm((n/3) * q), ncol = q))
  X2 = rbind(matrix(rnorm((n/3) * q1), ncol = q1),
             matrix(rnorm((n/3) * q1), ncol = q1),
             matrix(rnorm((n/3) * q1), ncol = q1))
  X3 = rbind(matrix(rnorm((n/3) * q2), ncol = q2),
             matrix(rnorm((n/3) * q2), ncol = q2),
             matrix(rnorm((n/3) * q2), ncol = q2))
  X.noise = rbind(matrix(rnorm((n/3) * q3), ncol = q3),
                  matrix(rnorm((n/3) * q3), ncol = q3),
                  matrix(rnorm((n/3) * q3), ncol = q3))
  X1[1:(n/3), 1:q] <- X1[1:(n/3), 1:q] - mu
  X1[(n/3+1):(2*n/3), 1:q] <- X1[(n/3+1):(2*n/3), 1:q]
  X1[(2*n/3+1):n, 1:q] <- X1[(2*n/3+1):n, 1:q] + mu
  index<-sample(1:n)
  X2[index[1:(n/3)], ] <- X2[index[1:(n/3)],]-mu1
  X2[index[(n/3+1):(2*n/3)], ] <- X2[index[(n/3+1):(2*n/3)], ]
  X2[index[(2*n/3+1):n], ] <- X2[index[(2*n/3+1):n], ]+mu1
  for(i in 1:ncol(X3)){
    X3[,i]<-X3[,i]+rnorm(1,0,var_g)*x1
  }
  X<-cbind(X1,X2,X3,X.noise)
  data<-list(x=x,y=y,y.ind=y.ind,G=X)
  return(data)
}

Jaccard.index = function(TrueFeature, SelectedFeature1){
  u = length(union(which(SelectedFeature1 == 1), which(TrueFeature == 1)))
  is = length(intersect(which(SelectedFeature1 == 1), which(TrueFeature == 1)))
  J1 = is/u
  return(J1)
}

SimUnbalance<-function(n,beta1,beta2,q,q1,q2,q3,c1,var_g,mu,mu1,sigma_y,frac3=0.5){
  #simulate covariates from N(0,1)
  x1<-rnorm(n,0,1)
  x2<-rnorm(n,0,1)
  x<-cbind(x1,x2)
  #simulate y
  y<-rep(NA,n)
  n1=floor(n*((1-frac3)/2))
  n3=n-2*n1
  y[1:n1]<-rnorm(n1,0+beta1*x1[1:n1]+beta2*x2[1:n1],sigma_y)
  y[(n1+1):(2*n1)]<-rnorm(n1,c1+beta1*x1[(n1+1):(2*n1)]+beta2*x2[(n1+1):(2*n1)],sigma_y)
  y[(2*n1+1):n]<-rnorm(n3,2*c1+beta1*x1[(2*n1+1):n]+beta2*x2[(2*n1+1):n],sigma_y)
  #simulate x
  X1 = rbind(matrix(rnorm(n1 * q), ncol = q),
             matrix(rnorm(n1 * q), ncol = q),
             matrix(rnorm(n3 * q), ncol = q))
  X2 = rbind(matrix(rnorm(n1 * q1), ncol = q1),
             matrix(rnorm(n1 * q1), ncol = q1),
             matrix(rnorm(n3 * q1), ncol = q1))
  X3 = rbind(matrix(rnorm(n1 * q2), ncol = q2),
             matrix(rnorm(n1 * q2), ncol = q2),
             matrix(rnorm(n3 * q2), ncol = q2))
  X.noise = rbind(matrix(rnorm(n1 * q3), ncol = q3),
                  matrix(rnorm(n1 * q3), ncol = q3),
                  matrix(rnorm(n3 * q3), ncol = q3))
  X1[1:n1, 1:q] <- X1[1:n1, 1:q] - mu
  X1[(n1+1):(2*n1), 1:q] <- X1[(n1+1):(2*n1), 1:q]
  X1[(2*n1+1):n, 1:q] <- X1[(2*n1+1):n, 1:q] + mu
  index<-sample(1:n)
  X2[index[1:n1], ] <- X2[index[1:(n1)],]-mu1
  X2[index[(n1+1):(2*n1)], ] <- X2[index[(n1+1):(2*n1)], ]
  X2[index[(2*n1+1):n], ] <- X2[index[(2*n1+1):n], ]+mu1
  for(i in 1:ncol(X3)){
    X3[,i]<-X3[,i]+rnorm(1,0,var_g)*x1
  }
  X<-cbind(X1,X2,X3,X.noise)
  data<-list(x=x,y=y,G=X)
  return(data)
}

