# ogClust R package

The R package "ogClust" mainly implements two outcome-guided clustering methods for disease subtyping. The source code is under "master" branch.

## Install the pacakge

`library(devtools)` \
`install_github("wenjiaking/ogClust")` \
`library(ogClust)`

## ogClust with continuous outcome

### load lung disease data
GSE47460 contains gene expression profiling of chronic lung disease for the Lung Genomics Research Consortium. GPL14550 is one of the two platforms in this GEO dataset, and only COPD and ILD patients in this dataset are kept.

```
data('GSE47460_GPL14550')
G=GSE47460_GPL14550$Expression
X=GSE47460_GPL14550$Covariates
Y=GSE47460_GPL14550$outcome
```


### preprocessing
```
g.mean<-apply(G,1,mean)
cut.mean=quantile(g.mean,probs=0.5)
G=G[g.mean>cut.mean,]
g.sd=apply(G,1, sd)
cut.sd=quantile(g.sd,probs=0.5)
G=G[g.sd>=cut.sd,]
G<-t(G)
G<-scale(G)
```
### implement ogClust_GM method

```
n=nrow(G)
NG=ncol(G)
np=ncol(X)
K=3
lambda=0.001
beta_int = runif(np, 0, 3)
gamma_int = runif((K - 1) * (NG + 1), 0, 1)
beta0_int = runif(K, 0, 3)
sigma2_int = runif(1, 1, 3)
theta_int = c(beta_int, gamma_int, beta0_int, sigma2_int)
fit.res<-ogClust_GM(n=n, K=3, np=np, NG=NG, lambda=lambda, alpha=1, G=G, Y=Y, X=X,theta_int=theta_int)
```

### implement ogClust_WJL method
```
mod.kmeans<-kmeans(G,centers = 3,nstart = 50)
center<-t(mod.kmeans$centers)
s_G<-300
w<-0.7
w<-(s_G*w)/(s_G*w+1-w)
lambda<-32
fit.res = ogClust_WJL(x=as.matrix(X),G=t(G),y=Y,c_center=center,lambda=lambda,v_int=NULL,pi_int=NULL,K=3,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
```
## ogClust with survival outcome
This is gene expression profiling of METARBIC breast cancer dataset, containing 275 triple negative breast cancer samples.

### load breast cancer data

```
data('Data.Metabric')
G=Data.Metabric$Expression
X=as.matrix(Data.Metabric$covariate[,1])
Y=Data.Metabric$OS
delta<-Data.Metabric$OS.event
```
### preprocessing
```
Index.Mean<-order(apply(G,2,mean),decreasing = T)
G <- G[,Index.Mean[1:(ncol(G)/2)]]
Index.Sd<-order(apply(G,2,sd),decreasing = T)
G<-G[,Index.Sd[1:(ncol(G)/2)]]
G<-scale(G)
```
### implement ogClust_GM method

```
n=nrow(G)
NG=ncol(G)
np=ncol(X)
K=2
lambda=0.007
beta_int = runif(np, 0, 3)
gamma_int = runif((K - 1) * (NG + 1), 0, 1)
beta0_int = runif(K, 0, 3)
sigma2_int = runif(1, 1, 3)
theta_int = c(beta_int, gamma_int, beta0_int, sigma2_int)
fit.res<-ogClust_GM.Surv(n=n, K=K, np=np, NG=NG, lambda=lambda,delta = delta, alpha=1, G=G, Y=Y, X=X, theta_int=theta_int)
```
### implement ogClust_WJL method

```
mod.kmeans<-kmeans(G,centers = 2,nstart = 50)
center<-t(mod.kmeans$centers)
s_G<-300
w<-0.9
w<-(s_G*w)/(s_G*w+1-w)
lambda<-33
fit.res=ogClust_WJL.Surv(x=X,G=t(G),y=Y,y.ind=delta,c_center=center,lambda=lambda, v_int=NULL,pi_int=NULL,K=2,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
```

