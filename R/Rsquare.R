#' @title Get the Rsquare value for gene expression and outcome (continuous)
#' @param cluster cluster assignment
#' @param Y outcome
#' @param G Gene expression data (p*n matrix), where each row is a gene and each column is a sample
#' @param X covariates matrix

#'
#' @return a list with two components
#' \itemize{
#' \item{R2.outcome: }{R2 of outcome}
#' \item{R2.gene: }{R2 of p genes}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' #mod is the model from ogClust
#' cluster<-apply(mod$result_list$z,1,which.max)
#' #y is continuous outcome, G is gene expression matrix, x is covariates matrix.
#' Rsquare(cluster,Y = y,G = G,X = x)
#' }

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
