#' Title Fit ogClust_GM mixture model embedding the LDA for gene disease subtyping and with survival outcome
#'
#' @param K an integer defines the number of subgroups
#' @param lambda the regularization tuning parameter for sparsity
#' @param G the matrix for omics data. The rows are samples and columns are features.
#' @param X the vector of covariates
#' @param Y the vector of time variable
#' @param max_iter the maximum number of iterations of EM algorithm. Default value is 200
#' @param delta a binary indicator of censoring for time-to-event outcome. 0 means censored and 1 means event observed.
#' @param dist distribution of the survival time, defualt is \code{"loglogistic"}. The choices are \code{"weibull"},
#' \code{"exponential"}, \code{"gaussian"}, \code{"logistic"},\code{"lognormal"} and \code{"loglogistic"}.
#'
#' @details The ogClust_GM.LDA.Surv is a unified latent generative model to perform clustering constructed
#' from omics data \code{G} with the guidance of outcome \code{Y}, and with covariate \code{X} to account for
#' the variability that is not related to subgrouping. A modified EM algorithm is applied for
#' numeric computation such that the liklihood is maximized. A posterior probability is obtain
#' for each subject belonging to each cluster.
#'
#' Similar to ogClust_GM.Surv, ogClust_GM.LDA.Surv method performs feature selection, latent subtype characterization and outcome prediction simultaneously.
#' The difference is that ogClust_GM.LDA.Surv embeds the sparse linear discriminant analysis (LDA) to model gene disease subtyping.
#'
#' Time variable \code{Y} and censoring indicator \code{delta} together defines a single time-to-event outcome.
#' We use accelerated-failure time(AFT) model to model time-to-event outcome, \code{dist} defines
#' the distribution of survival time.
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
#' @import survival
#' @examples
#' \dontrun{
#'   data('Data.Metabric') #load dataset
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
#'   # number of clusters
#'   K=2
#'   # tuning parameter
#'   lambda=0.007
#'
#'   # fit ogClust
#'   fit.res<-ogClust_GM.LDA.Surv(K=K,lambda=lambda,delta = delta, G=G, Y=Y, X=X)
#'}

ogClust_GM.LDA.Surv <- function(K, lambda, G, Y, X, delta, dist = "loglogistic",max_iter=200) {
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
      label_new[empty.c]=sample(1:K,length(empty.c),replace =T)
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

diag.disc <-function(x, centroids, prior) {
  dd <- t(x) %*% centroids
  dd0 <- (rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
  scale(dd, as.numeric(dd0), FALSE) # this is -.5*||x_i - mu_k||^2+log(pi_k)
}
