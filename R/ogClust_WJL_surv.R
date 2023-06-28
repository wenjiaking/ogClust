#' @title ogClust_WJL mehtod to cluster with survival outcome
#' @param x covariates in the outcome association model (e.g., age, gender, BMI),  n*q matrix where n is sample size and q is number of covariates
#' @param G Gene expression data (p*n matrix), where each row is a gene and each column is a sample
#' @param y outcome, a vector of length n
#' @param y.ind the censor indicator for the outcome. 1 means events and 0 means censor.
#' @param c_center a n*K matrix which determines the initial clusters.
#' @param lambda a tuning parameter determines gene selection.
#' @param v_int initial variance for each gene, optional
#' @param pi_int initial probability for each cluster, optional
#' @param K Number of Clusters
#' @param max_iter maximum number of iterations
#' @param w_outcome weight for outcome
#' @param w_G weight for gene expression
#' @param z_int a initial cluster assignment matrix. A n*K matrix where each row denotes the probability
#' for a sample belongs to each cluster. If z_int is given, c_center will be ignored.
#'
#' @return
#' result_list: A list with final estimators.
#' \itemize{
#' \item{mu: }{p*K matrix where each row represents the cluster centers for each feature}
#' \item{sigma: }{a vector of length p, the standard deviation for each feature}
#' \item{z: }{n*K matrix where each row represents the probability the sample belongs to each cluster}
#' \item{pi: }{a vector of length K, mixed probability for each cluster}
#' \item{int_coef: }{a K+q vector, the first K elements are the intercept estimated for each cluster. The next q values are the common coefficient for the covariates.}
#' \item{int_sigma_coef: }{standard deviation estimation in the outcome association model.}
#' }
#' Other three outputs:
#' \itemize{
#' \item{BIC: }{BIC from the final model}
#' \item{lik: }{optimal log likelihood achieved}
#' \item{iter: }{Number of iterations before convergence}
#' }
#' @export
#' @import survival
#' @examples
#' \dontrun{
#'
#'   data('Data.Metabric') #load lung dataset
#'
#'   # extract gene expression G, covariate X, outcome Y
#'   G=Data.Metabric$Expression
#'   X=as.matrix(Data.Metabric$covariate[,1])
#'   Y=Data.Metabric$OS
#'   delta<-Data.Metabric$OS.event
#'   Index.Mean<-order(apply(G,2,mean),decreasing = T)
#'   G <- G[,Index.Mean[1:(ncol(G)/2)]]
#'   Index.Sd<-order(apply(G,2,sd),decreasing = T)
#'   G<-G[,Index.Sd[1:(ncol(G)/2)]]
#'   G<-scale(G)
#'   #initialize cluster centers:
#'   mod.kmeans<-kmeans(G,centers = 2,nstart = 50)
#'   center<-t(mod.kmeans$centers)
#'   s_G<-300
#'   w<-0.8
#'   w<-(s_G*w)/(s_G*w+1-w)
#'   lambda=0.007
#'   #implement ogClust_WJL.surv:
#'   fit.res=ogClust_WJL.Surv(x=X,G=t(G),y=Y,y.ind=delta,c_center=center,lambda=lambda,
#'   v_int=NULL,pi_int=NULL,K=2,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
#'}

ogClust_WJL.Surv<-function(x,G,y,y.ind,c_center=NULL,lambda,v_int=NULL,pi_int=NULL,K=2,max_iter=200,w_outcome=0.5,w_G=0.5,z_int=NULL){
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

  res<-list('result_list'=result_list,'BIC'=BIC,'lik'=max_lik,'iter'=iter)
  return(res)


}

