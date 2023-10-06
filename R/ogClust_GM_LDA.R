#' Title Fit ogClust_GM mixture model embedding the LDA for gene disease subtyping and with continuous outcome
#'
#' @param K an integer defines the number of subgroups
#' @param lambda the regularization tuning parameter for sparsity
#' @param G the matrix for omics data. The rows are samples and columns are features
#' @param Y the vector of a single outcome
#' @param X the vector of covariates
#' @param max_iter the maximum number of iterations of EM algorithm. Default value is 200
#' @details The ogClust_GM.LDA is a unified latent generative model to perform clustering constructed
#' from omics data \code{G} with the guidance of outcome \code{Y}, and with covariate \code{X} to account for
#' the variability that is not related to subgrouping. A modified EM algorithm is applied for
#' numeric computation such that the liklihood is maximized. A posterior probability is obtain
#' for each subject belonging to each cluster.
#'
#' Similar to ogClust_GM, ogClust_GM.LDA method performs feature selection, latent subtype characterization and outcome prediction simultaneously.
#' The difference is that ogClust_GM.LDA embeds the sparse linear discriminant analysis (LDA) to model gene disease subtyping.
#'
#' @return An object with class \code{"ogClust"}
#' \itemize{
#'  \item{\code{par}}{ a list of parameter estimates}
#'  \item{\code{lambda}}{lambda}
#'  \item{\code{Y_prd}}{ predicted outcome}
#'  \item{\code{grp_assign}}{ prediced group assignement}
#' }
#' @export
#' @import penalizedLDA
#' @import adaHuber
#' @importFrom stats median
#'
#' @examples
#' \dontrun{
#'   data('GSE47460_GPL14550') #load lung dataset
#'
#'   # extract gene expression G, covariate X, outcome Y
#'   G=GSE47460_GPL14550$Expression
#'   X=GSE47460_GPL14550$Covariates
#'   Y=GSE47460_GPL14550$outcome
#'   g.mean<-apply(G,1,mean)
#'   cut.mean=quantile(g.mean,probs=0.5)
#'   G=G[g.mean>cut.mean,] # remove 50% lowest mean expression genes
#'   g.sd=apply(G,1, sd)
#'   cut.sd=quantile(g.sd,probs=0.5)
#'   G=G[g.sd>=cut.sd,] # further remove 50% lowest variance genes
#'   G<-t(G)
#'   G<-scale(G)
#'   # number of clusters
#'   K=3
#'   # tuning parameter
#'   lambda=0.001
#'
#'   # fit ogClust_GM.LDA
#'   fit.res<-ogClust_GM.LDA(K=K,lambda=lambda, G=G, Y=Y, X=X)
#'}

ogClust_GM.LDA <- function(K, lambda,G, Y, X, max_iter=200) {
  X = as.matrix(X)
  G = as.matrix(G)
  n=nrow(G)
  NG=ncol(G)
  np=ncol(X)
  EMout = EM.LDA(lambda = lambda, G = G, Y = Y, X = X, K = K, max_iter=max_iter)
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

EM.LDA <- function(lambda, G, Y, X, K, max_iter=200) {
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

diag.disc <-function(x, centroids, prior) {
  dd <- t(x) %*% centroids
  dd0 <- (rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
  scale(dd, as.numeric(dd0), FALSE) # this is -.5*||x_i - mu_k||^2+log(pi_k)
}

