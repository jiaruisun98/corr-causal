library(openxlsx)
source("Functions.R")

##This code is for the real case study. 

#Connection_diff1 provides the Nerwork in `Table1' and `Figure S2'
#Tau_sig provides the `Estimated Effect` in `Table1`
#CI_low and CI_high provides the 95% CI in `Table1`
#Ave-trt and Ace-cl provide Ave-trt and Ace-cl in `Table1`

data <- data.frame()
data <- array(data=0,dim=c(175,116,79))
#79 subjects,116ROIs,175 time points
#read time series data.

for(i in 50952:50962){
  path <- paste0("NYU_00",i,"_rois_aal",".1D")
  data[,,i-50951] <-as.matrix(read.table(file = path, header = TRUE))[1:175,]
}

for(i in 50964:51003){
  path <- paste0("NYU_00",i,"_rois_aal",".1D")
  data[,,i-50951-1] <-as.matrix(read.table(file = path, header = TRUE))[1:175,]
}


for(i in 51006:51021){
  path <- paste0("NYU_00",i,"_rois_aal",".1D")
  data[,,i-50951-3] <-as.matrix(read.table(file = path, header = TRUE))[1:175,]
}

for(i in 51023:51030){
  path <- paste0("NYU_00",i,"_rois_aal",".1D")
  data[,,i-50951-4] <-as.matrix(read.table(file = path, header = TRUE))[1:175,]
}


for(i in 51032:51035){
  path <- paste0("NYU_00",i,"_rois_aal",".1D")
  data[,,i-50951-5] <-as.matrix(read.table(file = path, header = TRUE))[1:175,]
}

X<-data

#read covariance
Covariance<-read.csv("phenotypic_NYU.csv")
#Covariance<-read.xlsx("COV_NYU.xlsx")
Covariance<-Covariance[which(Covariance[,"DX_GROUP"]==1),]
D<-Covariance[,"CURRENT_MED_STATUS"]

#deleting subjects with missing data
W<-Covariance[,c(5,6,8,9)]
W<-as.matrix(W)
X<-X[,,which(W[,3]!=-9999)]
D<-as.vector(D)
D<-D[which(W[,3]!=-9999)]
W<-W[which(W[,3]!=-9999),]
X<-X[,,-c(66,76)]
W<-W[-c(66,76),]
D<-D[-c(66,76)]

#parameters
p<-116
n<-76
m<-175
q<-4
B<-1000
aph<-0.05
c<-0.1


set.seed(12345678)



##Estimating and inference procedure
#causal effect estimation 

Y.estimate = array(0, c(p, p, n))
for (i in 1 : n){
  Y.estimate[, , i] = cor(X[, , i])
}
logistReg = glm(D ~ W + 0, family = binomial)
beta.estimate = logistReg$coefficient
Prob.estimate = as.vector(exp(W %*% beta.estimate) / (1 + exp(W %*% beta.estimate)))

tau.estimate = matrix(0, p, p)
weight = D / Prob.estimate - (1 - D) / (1 - Prob.estimate)

ID1=which(D==1)
ID0=which(D==0)
tau.estimate1 = matrix(0, p, p)
tau.estimate0 = matrix(0, p, p)
 for (i in 1 : length(ID1)){
  tau.estimate1 = tau.estimate1 + Y.estimate[, , ID1[i]] * weight[ID1[i]]
}
tau.estimate1 = tau.estimate1 / n

for (i in 1 : length(ID0)){
  tau.estimate0 = tau.estimate0 + Y.estimate[, , ID0[i]] * weight[ID0[i]]
}
tau.estimate0 = tau.estimate0 / n

for (i in 1 : n){
  tau.estimate = tau.estimate + Y.estimate[, , i] * weight[i]
}
tau.estimate = tau.estimate / n

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

#standardized treatment effect

Tstat = sqrt(n) * abs(tau.estimate) * theta.var^(-0.5)

Index = matrix(1, p, p)
for (j in 1 : p){
  Index[j, j] = 0
}

Tstat0 = Tstat * Index

#Proposed multiple testing procedure

z = array(0, c(p, p, B))
for(b in 1 : B){
  g = rnorm(n)
  temp = 0
  for (i in 1 : n){
    temp = temp + g[i] * (eta[, , i] - tau.estimate)
  }
  z[, , b] = theta.var^(-0.5) * temp / sqrt(n) * Index
}

result<-matrix(0,p,p)
z.initial = z
Tstat.initial = Tstat0
repeat{
  Tstat.max = max(abs(Tstat.initial))
  index.temp = which(abs(Tstat0) == Tstat.max, arr.ind = TRUE)
  z.max = apply(abs(z.initial), 3, max)
  z.max.quan = quantile(z.max, 1-aph)
  if (Tstat.max < z.max.quan) break
  for (i in 1 : dim(index.temp)[1]){
    Tstat.initial[index.temp[i, 1], index.temp[i, 2]] = 0
    z.initial[index.temp[i, 1], index.temp[i, 2], ] = 0
    result[index.temp[i,1],index.temp[i,2]]=1
  }
}

size<-sum(result)/2
num_add<-floor(c*size/(1-c))
if(num_add>=1){
  test_replace<-sort(Tstat.initial,decreasing=TRUE)
  for(i in 1:(2*num_add)){
    index.temp = which(abs(Tstat0) == test_replace[i], arr.ind = TRUE)
    for (i in 1 : dim(index.temp)[1]){
      result[index.temp[i,1],index.temp[i,2]]=1
    }
  }}

#detected signals
Connection1<-result
Connection_diff1<-which(Connection1 ==1, arr.ind = TRUE)
print(Connection_diff1)

#CI for the causal effect for detected connections
Tau_sig=tau.estimate[Connection_diff1]
theta_sig=theta.var[Connection_diff1]
CI_low=Tau_sig-sqrt(theta_sig)/sqrt(n)*qnorm(0.975,0,1)
CI_high=Tau_sig+sqrt(theta_sig)/sqrt(n)*qnorm(0.975,0,1)

#Ave-trt+Ave-cl
Ave_trt<-tau.estimate1[Connection_diff1]
Ave_cl<-tau.estimate0[Connection_diff1]




##########Non-Causal Method############



p<-116
n<-76
m<-175
q<-4
B<-1000
aph<-0.05
c<-0.1



##Correlation estimation
Y.estimate = array(0, c(p, p, n))
for (i in 1 : n){
  Y.estimate[, , i] = cor(X[, , i])
}
Y.estimate.y = Y.estimate[,,which(D==1)]
Y.estimate.n = Y.estimate[,,which(D==0)]


m1 = sum(D==1)
m2 = sum(D==0)

##Difference between treatment and control
Difference = array(0,c(p,p))
for(i in 1:m1){
  Difference = Difference + 1/sqrt(m1)*Y.estimate.y[,,i]
}
for(i in 1:m2){
  Difference = Difference - sqrt(m1)/m2*Y.estimate.n[,,i]
}

Difference = abs(Difference)

Y.estimate.y.average = array(0,c(p,p))
Y.estimate.n.average = array(0,c(p,p))

for(i in 1:m1){
  Y.estimate.y.average = Y.estimate.y.average + Y.estimate.y[,,i]
}
Y.estimate.y.average = Y.estimate.y.average/m1

for(i in 1:m2){
  Y.estimate.n.average = Y.estimate.n.average + Y.estimate.n[,,i]
}
Y.estimate.n.average = Y.estimate.n.average/m2


Index = matrix(1, p, p)
for (j in 1 : p){
  Index[j, j] = 0
}

Difference0 = Difference * Index

##proposed procedure 
z = array(0, c(p, p, B))
for(b in 1 : B){
  g1 = rnorm(m1)
  temp = 0
  for (i in 1 : m1){
    temp = temp + g1[i] * (Y.estimate.y[, , i] - Y.estimate.y.average)/sqrt(m1)
  }
  g2 = rnorm(m2)
  for (i in 1 : m2){
    temp = temp + g2[i] * (Y.estimate.n[, , i] - Y.estimate.n.average)*sqrt(m1)/m2
  }
  
  
  
  z[, , b] = temp * Index
}

result<-matrix(0,p,p)

Tstat.initial = Difference0
z.initial = z

repeat{
  Tstat.max = max(abs(Tstat.initial))
  index.temp = which(abs(Difference0) == Tstat.max, arr.ind = TRUE)
  z.max = apply(abs(z.initial), 3, max)
  z.max.quan = quantile(z.max, 1-aph)
  if (Tstat.max < z.max.quan) break
  # reject = rbind(reject, index.temp)
  for (i in 1 : dim(index.temp)[1]){
    Tstat.initial[index.temp[i, 1], index.temp[i, 2]] = 0
    z.initial[index.temp[i, 1], index.temp[i, 2], ] = 0
    result[index.temp[i,1],index.temp[i,2]]=1
  }
  #  cat("Max Stat = ", Tstat.max, "quantile = ", z.max.quan, "\n")
}

size<-sum(result)/2
num_add<-floor(c*size/(1-c))
if(num_add>=1){
  test_replace<-sort(Tstat.initial,decreasing=TRUE)
  for(i in 1:(2*num_add)){
    index.temp = which(abs(Tstat0) == test_replace[i], arr.ind = TRUE)
    for (i in 1 : dim(index.temp)[1]){
      result[index.temp[i,1],index.temp[i,2]]=1
    }
  }}


result_diff2<-which(result ==1, arr.ind = TRUE)


####################Causal+BH Method###########################


Y.estimate = array(0, c(p, p, n))
for (i in 1 : n){
  Y.estimate[, , i] = cor(X[, , i])
}
logistReg = glm(D ~ W + 0, family = binomial)
beta.estimate = logistReg$coefficient
Prob.estimate = as.vector(exp(W %*% beta.estimate) / (1 + exp(W %*% beta.estimate)))

tau.estimate = matrix(0, p, p)
weight = D / Prob.estimate - (1 - D) / (1 - Prob.estimate)
for (i in 1 : n){
  tau.estimate = tau.estimate + Y.estimate[, , i] * weight[i]
}
tau.estimate = tau.estimate / n



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


##Use the BH procedure instead of the proposed procedure

Tstat = sqrt(n) * abs(tau.estimate) * theta.var^(-0.5)
P=2*pnorm(Tstat,0,1,lower.tail = FALSE)

P1=P
for(i in 1:p){
  for(j in 1:p){
    if(j<=i){P1[i,j]=1}
  }
}
test=sort(P1)
num=1
while(test[num]<num*aph/(p*(p-1)/2)){
  num=num+1
}
num=num-1



result<-matrix(0,p,p)
if(num>1){
  for(i in 1:num){
    index.temp = which(P == test[i], arr.ind = TRUE)
    for (j in 1 : dim(index.temp)[1]){
      result[index.temp[j,1],index.temp[j,2]]=1
    }
  }
}


#detected signals
Connection_bh<-result
Connection_diff_bh<-which(Connection_bh ==1, arr.ind = TRUE)
print(Connection_diff_bh)

#CI for the causal effect for detected connections
Tau_sig=tau.estimate[Connection_diff1]
theta_sig=theta.var[Connection_diff1]
CI_low=Tau_sig-sqrt(theta_sig)/sqrt(n)*qnorm(0.975,0,1)
CI_high=Tau_sig+sqrt(theta_sig)/sqrt(n)*qnorm(0.975,0,1)

#Ave-trt+Ave-cl
Ave_trt<-tau.estimate1[Connection_diff1]
Ave_cl<-tau.estimate0[Connection_diff1]












