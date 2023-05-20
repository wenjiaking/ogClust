#' @title efficient lambda grids for a given K, for survival outcome
#' @param lambda1 a large lambda which gives close to 0 informative genes
#' @param lambda2 a small lambda which includes almost all the genes (default is lambda2=0)
#' @param iteration Number of iterations
#' @param Y outcome
#' @param G gene expression matrix (row is samples, column is genes)
#' @param X covariates matrix
#' @param center initial center for ogClust
#' @param w weight of outcome
#' @param K Number of clusters
#' @param delta censor indicator of survival outcome
#' @return a list with two components
#' \itemize{
#' \item{lambda: }{generated lambda grids}
#' \item{num: }{number of informative genes corresponding to each lambda}
#' }
#' @export
#'


region_lambda_surv<-function (lambda1 =18,lambda2=0,iteration = 10, Y, G,X,center,w,K=K,delta){
  lambda.vector <- c(lambda1, lambda2)
  num.vector <- rep(-1, length(lambda.vector))
  for (i in 1:length(num.vector)) {
    mod1 <- ogClust_WJL.Surv(x=X,G=t(G),y=Y,y.ind = delta,c_center=center,lambda=lambda.vector[i],v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
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
    mod_new<-ogClust_WJL.Surv(x=X,G=t(G),y=Y,y.ind = delta,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)

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
  mod_new<-ogClust_WJL.Surv(x=X,G=t(G),y=Y,y.ind = delta,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
  num.vector[which(lambda.vector == lambda)] <- sum(apply(mod_new$result_list$mu,1,function(x){length(unique(x))})!=1)
  res<-list(lambda=lambda.vector,num=num.vector)
  return(res)
}
