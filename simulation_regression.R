library(MASS)
library(glmnet)
library(iterators)
library(foreach)
library(doParallel)
library(doRNG)
########################################################################
#`result' provides the analysis result used to construct Figure S1
#each row of `result' represent a parameter setting.
#column 1 of `result` is the power of regular lasso
#column 2 of `result` is the false positive rate of regular lasso
#column 3 of `result` is the false discovery rate (fdr) of regular lasso
#column 4 of `result` is the FDPEx of regular lasso

#column 5 of `result` is the power of debiased lasso
#column 6 of `result` is the false positive rate of debiased lasso
#column 7 of `result` is the false discovery rate (fdr) of debiased lasso
#column 8 of `result` is the FDPEx of debiased lasso

#The simulation results are further organized in `simulation-result.xlsx` to produce Figure S1.

#########################################################################
ARcov = function(p, rho){
  cov0 = matrix(0, p, p)
  for(i in 1 : p){
    for (j in 1 : p){
      cov0[i, j] = rho^(abs(i - j))
    }
  }
  return(cov0)
}


#function for data sampling; block diagonal structure
#meaning of parameters can be seen in the article
#Raw0 is the parameter for correlation structure, raw is strength of time dependence
#n is the number of subjects, p the dimension of covariance matrix, m the length of time series

Sample_1<-function(p,Raw0,q,Beta,gamma,n,m,s0,raw){
  Sigma0 = ARcov(p, Raw0)
  Sigmaw = ARcov(q, 0.7) * 0.5
  W = mvrnorm(n, rep(0, q), Sigmaw)
  Prob = as.vector(exp(W %*% Beta) / (1 + exp(W %*% Beta)))
  
  #--- Data simulation---
  D = c()
  for(i in 1 : n){
    D[i] = rbinom(1, 1, Prob[i])
  }
  
  epsilon = rnorm(n, 0, 0.1)
  delta = as.vector(exp(abs(W %*% gamma)) / (1 + exp(abs(W %*% gamma)))) - 0.5 + epsilon 
  # need to take absolute value here otherwise the mean of delta is close to zero, no average treatment effect.
  X = array(0, c(m, p, n))
  E = array(0, c(m, p, n))
  
  #block diagonal signal structure
  k=p/s0
  signa.vector=array(0,c(k,p))
  for(i in 1:k){
    signa.vector[i,(((i-1)*s0+1):(i*s0))]<-rep(1,s0)
  }
  signal<-t(signa.vector)%*%signa.vector
  signal=signal-diag(rep(1,p))
  
  
  for (i in 1 : n){
    Y0 = Sigma0
    Y1 = solve( solve(Sigma0) + signal* abs(delta[i])*Raw0 )
    Y = D[i] * Y1 + (1 - D[i]) * Y0
    
    # for time dependent data
    E[, , i] = mvrnorm(m, rep(0,p), Y)
    for(j in 1:m){
      if(j==1){X[j, , i]=E[j, , i]}
      if(j>=2){X[j, , i]=raw*X[j-1, , i]+sqrt(1-raw^2)*E[j, , i]}
    }
  }
  outcome<-list(W,D,X)
  # W the covariance matrix, D the treatment assignment, X the time series data.
  
  return(outcome)}



##function of the proposed procedure on regression parameters

node_regression_estimate_diag<-function(W,X,D,p,n,m,q,B,aph,c,s0){
  
  #lasso regresion penalty
  lambda = sqrt(2 * log(p)/n)
  
  #causal effect estimation - lasso
  Y.estimate = array(0, c(p, p, n))
  for (i in 1:n){
    for(j in 1:p){
      X_ = X[ , ,i]
      X_ = scale(X_, scale = FALSE)
      
      lasso = glmnet(x = X_[ ,-j], y = X_[, j], intercept = FALSE, 
                     standardize = FALSE)
      Coef = coef.glmnet(lasso, s = lambda)
      
      Y.estimate[j, -j, i] =  Coef[-1]
      Y.estimate[j, j, i] = -1
    }
  }
  
  
  #causal effect estimation - debiased lasso
  
  Y.estimate_d = array(0, c(p, p, n))
  for (i in 1:n){
    for(j1 in 1:p){
      for(j2 in (1:p)[-j1]){
        X_ = X[ , ,i]
        X_ = scale(X_, scale = FALSE)
        
        lasso_debias = glmnet(x = X_[ ,-c(j1,j2)], y = X_[, j2], intercept = FALSE, 
                              standardize = FALSE)
        Coef = coef.glmnet(lasso_debias, s = lambda)
        gamma.debias<-Coef[-1]
        A<-mean(X_[ ,j2]*(X_[ ,j2]-X_[,-c(j1,j2)]%*%gamma.debias))
        Y.estimate_d[j1,j2,i] = Y.estimate[j1,j2,i] - 1/A*mean((X_[,]%*%Y.estimate[j1, ,i])*(X_[,j2]-X_[,-c(j1,j2)]%*%gamma.debias))
      }
    }
  }
  
  #the rest is the same as the correlation simulation
  
  logistReg = glm(D ~ W + 0, family = binomial)
  beta.estimate = logistReg$coefficient
  Prob.estimate = as.vector(exp(W %*% beta.estimate) / (1 + exp(W %*% beta.estimate)))
  
  #tau
  tau.estimate = matrix(0, p, p)
  weight = D / Prob.estimate - (1 - D) / (1 - Prob.estimate)
  for (i in 1 : n){
    tau.estimate = tau.estimate + Y.estimate[, , i] * weight[i]
  }
  tau.estimate = tau.estimate / n
  
  #tau_d
  tau.estimate_d = matrix(0, p, p)
  weight = D / Prob.estimate - (1 - D) / (1 - Prob.estimate)
  for (i in 1 : n){
    tau.estimate_d = tau.estimate_d + Y.estimate_d[, , i] * weight[i]
  }
  tau.estimate_d = tau.estimate_d / n
  
  
  #estimation for eta, the influence function and variance
  Fisher.Information = 0
  for (i in 1 : n){
    Fisher.Information = Fisher.Information + Prob.estimate[i] * (1 - Prob.estimate[i]) * W[i, ] %*% t(W[i, ])
  }
  Fisher.Information = Fisher.Information / n
  Theta = solve(Fisher.Information)
  weight.H = D * (1 - Prob.estimate) / Prob.estimate + (1 - D) * Prob.estimate / (1 - Prob.estimate)
  H = array(0, c(p, p, q))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      H[j1, j2, ] = colMeans(weight.H * Y.estimate[j1, j2, ] * W)
    }
  }
  
  eta = array(0, c(p, p, n))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      for (i in 1 : n){
        eta[j1, j2, i] = Y.estimate[j1, j2, i] * weight[i] - H[j1, j2, ] %*% Theta %*% W[i, ] * (D[i] - Prob.estimate[i])
      }
    }
  }
  
  theta.var = matrix(0, p, p)
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      theta.var[j1, j2] = mean((eta[j1, j2, ] - tau.estimate[j1, j2])^2)
    }
  }
  
  
  #estimation for eta_d, the influence function and variance
  Fisher.Information_d = 0
  for (i in 1 : n){
    Fisher.Information_d = Fisher.Information_d + Prob.estimate[i] * (1 - Prob.estimate[i]) * W[i, ] %*% t(W[i, ])
  }
  Fisher.Information_d = Fisher.Information_d / n
  Theta_d = solve(Fisher.Information_d)
  weight.H = D * (1 - Prob.estimate) / Prob.estimate + (1 - D) * Prob.estimate / (1 - Prob.estimate)
  H_d = array(0, c(p, p, q))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      H_d[j1, j2, ] = colMeans(weight.H * Y.estimate_d[j1, j2, ] * W)
    }
  }
  
  eta_d = array(0, c(p, p, n))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      for (i in 1 : n){
        eta_d[j1, j2, i] = Y.estimate_d[j1, j2, i] * weight[i] - H_d[j1, j2, ] %*% Theta_d %*% W[i, ] * (D[i] - Prob.estimate[i])
      }
    }
  }
  
  theta.var_d = matrix(0, p, p)
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      theta.var_d[j1, j2] = mean((eta_d[j1, j2, ] - tau.estimate_d[j1, j2])^2)
    }
  }
  
  
  
  #standardized treatment effect
  Tstat = sqrt(n) * abs(tau.estimate) * theta.var^(-0.5)
  #set variance term, the diagonal to be zero
  Index = matrix(1, p, p)
  for (j in 1 : p){
    Index[j, j] = 0
  }
  Tstat0 = Tstat * Index
  Tstat0[which(abs(tau.estimate)<0.01&theta.var<0.01,arr.ind = TRUE)]<-0
  
  ##using c=0.01 to exclude tau<c, theta<c
  id_exclude = which(abs(tau.estimate)<0.01&theta.var<0.01,arr.ind = TRUE)
  
  #Proposed multiple testing procedure
  z = array(0, c(p, p, B))
  for(b in 1 : B){
    g = rnorm(n)
    temp = 0
    for (i in 1 : n){
      temp = temp + g[i] * (eta[, , i] - tau.estimate)
    }
    z[, , b] = theta.var^(-0.5) * temp / sqrt(n) * Index
    
    id_exclude_ = cbind(id_exclude,rep(b,nrow(id_exclude)))
    z[id_exclude_] = 0
  }
  
  result<-matrix(0,p,p)
  z.initial = z
  Tstat.initial = Tstat0
  repeat{
    Tstat.max = max(abs(Tstat.initial))
    index.temp = which(abs(Tstat0) == Tstat.max, arr.ind = TRUE)
    z.max = apply(abs(z.initial), 3, max)
    #z.max[which(abs(tau.estimate)<0.01&theta.var<0.01,arr.ind = TRUE)]<-0
    z.max.quan = quantile(z.max, 1-aph)
    if (z.max.quan == 0) break
    if (Tstat.max < z.max.quan) break
    for (i in 1 : dim(index.temp)[1]){
      Tstat.initial[index.temp[i, 1], index.temp[i, 2]] = 0
      z.initial[index.temp[i, 1], index.temp[i, 2], ] = 0
      result[index.temp[i,1],index.temp[i,2]]=1
    }
  }
  
  #Augment
  size<-sum(result)
  num_add<-floor(c*size/(1-c))
  if(num_add>=1){
    test_replace<-sort(Tstat.initial,decreasing=TRUE)
    for(i in 1:(num_add)){
      index.temp = which(abs(Tstat0) == test_replace[i], arr.ind = TRUE)
      for (i in 1 : dim(index.temp)[1]){
        result[index.temp[i,1],index.temp[i,2]]=1
      }
    }}
  
  
  
  #standardized treatment effect------debiased
  Tstat_d = sqrt(n) * abs(tau.estimate_d) * theta.var_d^(-0.5)
  #set variance term, the diagonal to be zero
  Index = matrix(1, p, p)
  for (j in 1 : p){
    Index[j, j] = 0
  }
  Tstat0_d = Tstat_d * Index
  Tstat0_d[which(abs(tau.estimate_d)<0.01&theta.var_d<0.01,arr.ind = TRUE)]<-0
  
  ##using c=0.01 to exclude tau<c, theta<c
  id_exclude_d = which(abs(tau.estimate_d)<0.01&theta.var_d<0.01,arr.ind = TRUE)
  
  #Proposed multiple testing procedure
  z_d = array(0, c(p, p, B))
  for(b in 1 : B){
    g = rnorm(n)
    temp = 0
    for (i in 1 : n){
      temp = temp + g[i] * (eta_d[, , i] - tau.estimate_d)
    }
    z_d[, , b] = theta.var_d^(-0.5) * temp / sqrt(n) * Index
    
    id_exclude_ = cbind(id_exclude_d,rep(b,nrow(id_exclude_d)))
    z_d[id_exclude_] = 0
  }
  
  result_d<-matrix(0,p,p)
  z.initial_d = z_d
  Tstat.initial_d = Tstat0_d
  repeat{
    Tstat.max_d = max(abs(Tstat.initial_d))
    index.temp_d = which(abs(Tstat0_d) == Tstat.max_d, arr.ind = TRUE)
    z.max_d = apply(abs(z.initial_d), 3, max)
    #z.max_d[which(abs(tau.estimate_d)<0.01&theta.var_d<0.01,arr.ind = TRUE)]<-0
    z.max.quan_d = quantile(z.max_d, 1-aph)
    if (z.max.quan_d == 0) break
    if (Tstat.max_d < z.max.quan_d) break
    for (i in 1 : dim(index.temp_d)[1]){
      Tstat.initial_d[index.temp_d[i, 1], index.temp_d[i, 2]] = 0
      z.initial_d[index.temp_d[i, 1], index.temp_d[i, 2], ] = 0
      result_d[index.temp_d[i,1],index.temp_d[i,2]]=1
    }
  }
  
  #Augment
  size_d<-sum(result_d)
  num_add_d<-floor(c*size_d/(1-c))
  if(num_add_d>=1){
    test_replace_d<-sort(Tstat.initial_d,decreasing=TRUE)
    for(i in 1:(num_add_d)){
      index.temp_d = which(abs(Tstat0_d) == test_replace_d[i], arr.ind = TRUE)
      for (i in 1 : dim(index.temp_d)[1]){
        result_d[index.temp_d[i,1],index.temp_d[i,2]]=1
      }
    }}
  
  
  
  
  
  #This part for calculation of simulation results
  p_true<-0
  p_false<-0
  fdr<-0
  
  
  k=p/s0
  signa.vector=array(0,c(k,p))
  for(i in 1:k){
    signa.vector[i,(((i-1)*s0+1):(i*s0))]<-rep(1,s0)
  }
  signal<-t(signa.vector)%*%signa.vector
  
  for(i in 1:p){
    for(j in (1:p)[-i]){
      if(signal[i,j]==1&&result[i,j]==1){p_true=p_true+1}
      if(signal[i,j]==0&&result[i,j]==1){p_false=p_false+1}
    }}
  
  if((p_true+p_false)>0){fdr<-p_false/(p_false+p_true)}else{
    fdr<-0}
  p_true_rate<-p_true/((sum(signal)-p))
  p_false_rate<-p_false/((p*p-sum(signal)))
  
  #This part for calculation of simulation results - debiased
  p_true_d<-0
  p_false_d<-0
  fdr_d<-0
  
  
  k=p/s0
  signa.vector=array(0,c(k,p))
  for(i in 1:k){
    signa.vector[i,(((i-1)*s0+1):(i*s0))]<-rep(1,s0)
  }
  signal<-t(signa.vector)%*%signa.vector
  
  for(i in 1:p){
    for(j in (1:p)[-i]){
      if(signal[i,j]==1&&result_d[i,j]==1){p_true_d=p_true_d+1}
      if(signal[i,j]==0&&result_d[i,j]==1){p_false_d=p_false_d+1
      print(c(i,j))}
    }}
  
  if((p_true_d+p_false_d)>0){fdr_d<-p_false_d/(p_false_d+p_true_d)}else{
    fdr_d<-0}
  p_true_rate_d<-p_true_d/(sum(signal)-p)
  p_false_rate_d<-p_false_d/((p*p-sum(signal)))
  
  
  #print results
  #print(Tstat)
  #print(result)
  #print(p_true_rate)
  #print(p_false_rate)
  #print(fdr)
  re<-c(p_true_rate,p_false_rate,fdr,p_true_rate_d,p_false_rate_d,fdr_d)
  return(re)}







result<-matrix(0,nrow=48,ncol=8)
colnames(result)<-c("power-debiased","fpr-debiased","fdr-debiased","FDPex-debiased","power","fpr","fdr","FDPex")

##Simulation1:lasso+block diagonal+time dependencce raw=0

raw=0

#parameters
m=300
q=4
gamma=rep(0.4,q)
Beta=rep(0.5,q)
c=0.1
B=1000
aph=0.05
s0=10


#1-1
Raw0=0.3
n=200
p=50


registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re1 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[1,c(1,2,3,5,6,7)]<-apply(l_re1,2,mean)
result[1,4]<-mean(l_re1[,3]>c)
result[1,8]<-mean(l_re1[,6]>c)


#1-2

p=50
n=100
registerDoParallel(30)

time_start <- Sys.time()

l_re2 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[2,c(1,2,3,5,6,7)]<-apply(l_re2,2,mean)
result[2,4]<-mean(l_re2[,3]>c)
result[2,8]<-mean(l_re2[,6]>c)

#1-3

n<-200
p<-100
registerDoParallel(30)

time_start <- Sys.time()

l_re3 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[3,c(1,2,3,5,6,7)]<-apply(l_re3,2,mean)
result[3,4]<-mean(l_re3[,3]>c)
result[3,8]<-mean(l_re3[,6]>c)

#1-4
n=100
p=100
registerDoParallel(30)

time_start <- Sys.time()

l_re4 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


stopImplicitCluster() 

result[4,c(1,2,3,5,6,7)]<-apply(l_re4,2,mean)
result[4,4]<-mean(l_re4[,3]>c)
result[4,8]<-mean(l_re4[,6]>c)




#1-5

n=200
p=50
Raw0=0.25


registerDoParallel(30)

time_start <- Sys.time()

l_re5 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


stopImplicitCluster() 

result[5,c(1,2,3,5,6,7)]<-apply(l_re5,2,mean)
result[5,4]<-mean(l_re5[,3]>c)
result[5,8]<-mean(l_re5[,6]>c)


#1-6

n=100
p=50

registerDoParallel(30)

time_start <- Sys.time()

l_re6 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


stopImplicitCluster() 

result[6,c(1,2,3,5,6,7)]<-apply(l_re6,2,mean)
result[6,4]<-mean(l_re6[,3]>c)
result[6,8]<-mean(l_re6[,6]>c)




#1-7

n=200
p=100

registerDoParallel(30)

time_start <- Sys.time()

l_re7 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


stopImplicitCluster() 

result[7,c(1,2,3,5,6,7)]<-apply(l_re7,2,mean)
result[7,4]<-mean(l_re7[,3]>c)
result[7,8]<-mean(l_re7[,6]>c)




#1-8

n=100
p=100

registerDoParallel(30)

time_start <- Sys.time()

l_re8 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


stopImplicitCluster() 

result[8,c(1,2,3,5,6,7)]<-apply(l_re8,2,mean)
result[8,4]<-mean(l_re8[,3]>c)
result[8,8]<-mean(l_re8[,6]>c)




#1-9

Raw0=0.2
n=200
p=50


registerDoParallel(30)

time_start <- Sys.time()

l_re9 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


stopImplicitCluster() 

result[9,c(1,2,3,5,6,7)]<-apply(l_re9,2,mean)
result[9,4]<-mean(l_re9[,3]>c)
result[9,8]<-mean(l_re9[,6]>c)




#1-10

n=100
p=50


registerDoParallel(30)

time_start <- Sys.time()

l_re10 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


stopImplicitCluster() 

result[10,c(1,2,3,5,6,7)]<-apply(l_re10,2,mean)
result[10,4]<-mean(l_re10[,3]>c)
result[10,8]<-mean(l_re10[,6]>c)




#1-11  
n=200
p=100

registerDoParallel(30)

time_start <- Sys.time()

l_re11 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


stopImplicitCluster() 

result[11,c(1,2,3,5,6,7)]<-apply(l_re11,2,mean)
result[11,4]<-mean(l_re11[,3]>c)
result[11,8]<-mean(l_re11[,6]>c)




#1-12

n=100
p=100

registerDoParallel(30)

time_start <- Sys.time()

l_re12 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


stopImplicitCluster() 

result[12,c(1,2,3,5,6,7)]<-apply(l_re12,2,mean)
result[12,4]<-mean(l_re12[,3]>c)
result[12,8]<-mean(l_re12[,6]>c)






##Simulation2:proposed method+block diagonal+time dependencce raw=0.3

raw=0.3

#parameters
m=300
q=4
gamma=rep(0.4,q)
Beta=rep(0.5,q)
c=0.1
B=1000
aph=0.05
s0=10


#2-1
Raw0=0.3
n=200
p=50



registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re13 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[13,c(1,2,3,5,6,7)]<-apply(l_re13,2,mean)
result[13,4]<-mean(l_re13[,3]>c)
result[13,8]<-mean(l_re13[,6]>c)




#2-2

p=50
n=100

registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re14 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[14,c(1,2,3,5,6,7)]<-apply(l_re14,2,mean)
result[14,4]<-mean(l_re14[,3]>c)
result[14,8]<-mean(l_re14[,6]>c)

#2-3

n<-200
p<-100

registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re15 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[15,c(1,2,3,5,6,7)]<-apply(l_re15,2,mean)
result[15,4]<-mean(l_re15[,3]>c)
result[15,8]<-mean(l_re15[,6]>c)


#2-4
n=100
p=100

registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re16 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[16,c(1,2,3,5,6,7)]<-apply(l_re16,2,mean)
result[16,4]<-mean(l_re16[,3]>c)
result[16,8]<-mean(l_re16[,6]>c)




#2-5

n=200
p=50
Raw0=0.25

registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re17 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[17,c(1,2,3,5,6,7)]<-apply(l_re17,2,mean)
result[17,4]<-mean(l_re17[,3]>c)
result[17,8]<-mean(l_re17[,6]>c)




#2-6

n=100
p=50


registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re18 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[18,c(1,2,3,5,6,7)]<-apply(l_re18,2,mean)
result[18,4]<-mean(l_re18[,3]>c)
result[18,8]<-mean(l_re18[,6]>c)





#2-7

n=200
p=100


registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re19 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[19,c(1,2,3,5,6,7)]<-apply(l_re19,2,mean)
result[19,4]<-mean(l_re19[,3]>c)
result[19,8]<-mean(l_re19[,6]>c)



#2-8

n=100
p=100


registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re20 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[20,c(1,2,3,5,6,7)]<-apply(l_re20,2,mean)
result[20,4]<-mean(l_re20[,3]>c)
result[20,8]<-mean(l_re20[,6]>c)




#2-9

Raw0=0.2
n=200
p=50


registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re21 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[21,c(1,2,3,5,6,7)]<-apply(l_re21,2,mean)
result[21,4]<-mean(l_re21[,3]>c)
result[21,8]<-mean(l_re21[,6]>c)


#2-10

n=100
p=50

registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re22 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[22,c(1,2,3,5,6,7)]<-apply(l_re22,2,mean)
result[22,4]<-mean(l_re22[,3]>c)
result[22,8]<-mean(l_re22[,6]>c)





#2-11  
n=200
p=100


registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re23 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[23,c(1,2,3,5,6,7)]<-apply(l_re23,2,mean)
result[23,4]<-mean(l_re23[,3]>c)
result[23,8]<-mean(l_re23[,6]>c)





#2-12

n=100
p=100

registerDoParallel(30)

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

l_re24 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-node_regression_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

stopImplicitCluster() 

result[24,c(1,2,3,5,6,7)]<-apply(l_re24,2,mean)
result[24,4]<-mean(l_re24[,3]>c)
result[24,8]<-mean(l_re24[,6]>c)
