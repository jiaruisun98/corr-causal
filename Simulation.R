source("Functions.R")

library(MASS)
library(iterators)
#library(foreach)
library(doParallel)

registerDoParallel(30)
#################################################################################
#Simulation_block_proposed and Simulation_block_BH are simulation results used to produce Figure 1.
#Simulation_off_proposed and Simulation_off_BH are simulation results used to produce Figure2.
#Each row of Simulation_block_proposed, Simulation_block_BH, Simulation_off_proposed and Simulation_off_BH represent a parameter setting.

#column 1 of `Simulation_block_proposed` is the power of proposed method under block diagonal setting.
#column 2 of `Simulation_block_proposed` is the false positive rate of proposed method under block diagonal setting.
#column 3 of `Simulation_block_proposed` is the false discovery rate (fdr) of proposed method under block diagonal setting.
#column 4 of `Simulation_block_proposed` is the FDPEx of proposed method under block diagonal setting.

#column 1 of `Simulation_block_BH` is the power of BH method under block diagonal setting.
#column 2 of `Simulation_block_BH` is the false positive rate of BH method under block diagonal setting.
#column 3 of `Simulation_block_BH` is the false discovery rate (fdr) of BH method under block diagonal setting.
#column 4 of `Simulation_block_BH` is the FDPEx of BH method under block diagonal setting.

#column 1 of `Simulation_off_proposed` is the power of proposed method under off diagonal setting.
#column 2 of `Simulation_off_proposed` is the false positive rate of proposed method under off diagonal setting.
#column 3 of `Simulation_off_proposed` is the false discovery rate (fdr) of proposed method under off diagonal setting.
#column 4 of `Simulation_off_proposed` is the FDPEx of proposed method under off diagonal setting.

#column 1 of `Simulation_off_BH` is the power of BH method under off diagonal setting.
#column 2 of `Simulation_off_BH` is the false positive rate of BH method under off diagonal setting.
#column 3 of `Simulation_off_BH` is the false discovery rate (fdr) of BH method under block diagonal setting.
#column 4 of `Simulation_off_BH` is the FDPEx of BH method under off diagonal setting.

#The simulation results are further organized in `simulation-result.xlsx` to produce Figure 1 and Figure 2.

#################################################################################
Simulation_block_proposed<-matrix(0,nrow=24,ncol=4)
Simulation_block_BH<-matrix(0,nrow=24,ncol=4)
Simulation_off_proposed<-matrix(0,nrow=24,ncol=4)
Simulation_off_BH<-matrix(0,nrow=24,ncol=4)
##############################
##Simulation1:proposed method+block diagonal+time dependencce raw=0

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

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

a_re1 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[1,c(1,2,3)]<-apply(a_re1,2,mean)
Simulation_block_proposed[1,4]<-mean(a_re1[,3]>c)

#1-2
 
p=50
n=100

time_start <- Sys.time()

a_re2 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[2,c(1,2,3)]<-apply(a_re2,2,mean)
Simulation_block_proposed[2,4]<-mean(a_re2[,3]>c)

#1-3

n<-200
p<-100

time_start <- Sys.time()

a_re3 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[3,c(1,2,3)]<-apply(a_re3,2,mean)
Simulation_block_proposed[3,4]<-mean(a_re3[,3]>c)

#1-4
n=100
p=100

time_start <- Sys.time()

a_re4 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[4,c(1,2,3)]<-apply(a_re4,2,mean)
Simulation_block_proposed[4,4]<-mean(a_re4[,3]>c)

#1-5
  
n=200
p=50
Raw0=0.25


time_start <- Sys.time()

a_re5 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[5,c(1,2,3)]<-apply(a_re5,2,mean)
Simulation_block_proposed[1,5]<-mean(a_re5[,3]>c)

#1-6
  
n=100
p=50

time_start <- Sys.time()

a_re6 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[6,c(1,2,3)]<-apply(a_re6,2,mean)
Simulation_block_proposed[6,4]<-mean(a_re6[,3]>c)

#1-7
  
n=200
p=100

time_start <- Sys.time()

a_re7 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[7,c(1,2,3)]<-apply(a_re7,2,mean)
Simulation_block_proposed[7,4]<-mean(a_re7[,3]>c)

#1-8
  
n=100
p=100

time_start <- Sys.time()

a_re8 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[8,c(1,2,3)]<-apply(a_re8,2,mean)
Simulation_block_proposed[8,4]<-mean(a_re8[,3]>c)

#1-9
  
Raw0=0.2
n=200
p=50


time_start <- Sys.time()

a_re9 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[9,c(1,2,3)]<-apply(a_re9,2,mean)
Simulation_block_proposed[9,4]<-mean(a_re9[,3]>c)

#1-10
  
n=100
p=50


time_start <- Sys.time()

a_re10 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[10,c(1,2,3)]<-apply(a_re10,2,mean)
Simulation_block_proposed[10,4]<-mean(a_re10[,3]>c)

#1-11  
n=200
p=100

time_start <- Sys.time()

a_re11 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[11,c(1,2,3)]<-apply(a_re11,2,mean)
Simulation_block_proposed[11,4]<-mean(a_re11[,3]>c)

#1-12
  
n=100
p=100

time_start <- Sys.time()

a_re12 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[12,c(1,2,3)]<-apply(a_re12,2,mean)
Simulation_block_proposed[12,4]<-mean(a_re12[,3]>c)

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

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

b_re1 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[13,c(1,2,3)]<-apply(b_re1,2,mean)
Simulation_block_proposed[13,4]<-mean(b_re1[,3]>c)

#2-2

p=50
n=100

time_start <- Sys.time()

b_re2 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[14,c(1,2,3)]<-apply(b_re2,2,mean)
Simulation_block_proposed[14,4]<-mean(b_re2[,3]>c)

#2-3

n<-200
p<-100

time_start <- Sys.time()

b_re3 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[15,c(1,2,3)]<-apply(b_re3,2,mean)
Simulation_block_proposed[15,4]<-mean(b_re3[,3]>c)
#2-4
n=100
p=100

time_start <- Sys.time()

b_re4 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[16,c(1,2,3)]<-apply(b_re4,2,mean)
Simulation_block_proposed[16,4]<-mean(b_re4[,3]>c)

#2-5

n=200
p=50
Raw0=0.25


time_start <- Sys.time()

b_re5 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[17,c(1,2,3)]<-apply(b_re5,2,mean)
Simulation_block_proposed[17,4]<-mean(b_re5[,3]>c)

#2-6

n=100
p=50

time_start <- Sys.time()

b_re6 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[18,c(1,2,3)]<-apply(b_re6,2,mean)
Simulation_block_proposed[18,4]<-mean(b_re6[,3]>c)


#2-7

n=200
p=100

time_start <- Sys.time()

b_re7 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[19,c(1,2,3)]<-apply(b_re7,2,mean)
Simulation_block_proposed[19,4]<-mean(b_re7[,3]>c)

#2-8

n=100
p=100

time_start <- Sys.time()

b_re8 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[20,c(1,2,3)]<-apply(b_re8,2,mean)
Simulation_block_proposed[20,4]<-mean(b_re8[,3]>c)

#2-9

Raw0=0.2
n=200
p=50


time_start <- Sys.time()

b_re9 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[21,c(1,2,3)]<-apply(b_re9,2,mean)
Simulation_block_proposed[21,4]<-mean(b_re9[,3]>c)

#2-10

n=100
p=50


time_start <- Sys.time()

b_re10 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[22,c(1,2,3)]<-apply(b_re10,2,mean)
Simulation_block_proposed[22,4]<-mean(b_re10[,3]>c)

#2-11  
n=200
p=100

time_start <- Sys.time()

b_re11 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[23,c(1,2,3)]<-apply(b_re11,2,mean)
Simulation_block_proposed[23,4]<-mean(b_re11[,3]>c)

#2-12

n=100
p=100

time_start <- Sys.time()

b_re12 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_proposed[24,c(1,2,3)]<-apply(b_re12,2,mean)
Simulation_block_proposed[24,4]<-mean(b_re12[,3]>c)

##Simulation3:BH method+block diagonal+time dependencce raw=0


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


#3-1
Raw0=0.3
n=200
p=50

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

c_re1 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[1,c(1,2,3)]<-apply(c_re1,2,mean)
Simulation_block_BH[1,4]<-mean(c_re1[,3]>c)

#3-2

p=50
n=100

time_start <- Sys.time()

c_re2 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[2,c(1,2,3)]<-apply(c_re2,2,mean)
Simulation_block_BH[2,4]<-mean(c_re2[,3]>c)

#3-3

n<-200
p<-100

time_start <- Sys.time()

c_re3 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[3,c(1,2,3)]<-apply(c_re3,2,mean)
Simulation_block_BH[3,4]<-mean(c_re3[,3]>c)

#3-4
n=100
p=100

time_start <- Sys.time()

c_re4 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[4,c(1,2,3)]<-apply(c_re4,2,mean)
Simulation_block_BH[4,4]<-mean(c_re4[,3]>c)

#3-5

n=200
p=50
Raw0=0.25


time_start <- Sys.time()

c_re5 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_block_BH[5,c(1,2,3)]<-apply(c_re5,2,mean)
Simulation_block_BH[5,4]<-mean(c_re5[,3]>c)
#3-6

n=100
p=50

time_start <- Sys.time()

c_re6 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[6,c(1,2,3)]<-apply(c_re6,2,mean)
Simulation_block_BH[6,4]<-mean(c_re6[,3]>c)


#3-7

n=200
p=100

time_start <- Sys.time()

c_re7 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_block_BH[7,c(1,2,3)]<-apply(c_re7,2,mean)
Simulation_block_BH[7,4]<-mean(c_re7[,3]>c)

#3-8

n=100
p=100

time_start <- Sys.time()

c_re8 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[8,c(1,2,3)]<-apply(c_re8,2,mean)
Simulation_block_BH[8,4]<-mean(c_re8[,3]>c)

#3-9

Raw0=0.2
n=200
p=50


time_start <- Sys.time()

c_re9 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[9,c(1,2,3)]<-apply(c_re9,2,mean)
Simulation_block_BH[9,4]<-mean(c_re9[,3]>c)


#3-10

n=100
p=50


time_start <- Sys.time()

c_re10 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[10,c(1,2,3)]<-apply(c_re10,2,mean)
Simulation_block_BH[10,4]<-mean(c_re10[,3]>c)


#3-11  
n=200
p=100

time_start <- Sys.time()

c_re11 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_block_BH[11,c(1,2,3)]<-apply(c_re11,2,mean)
Simulation_block_BH[11,4]<-mean(c_re11[,3]>c)

#3-12

n=100
p=100

time_start <- Sys.time()

c_re12 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[12,c(1,2,3)]<-apply(c_re12,2,mean)
Simulation_block_BH[12,4]<-mean(c_re12[,3]>c)

##Simulation4:BH method+block diagonal+time dependencce raw=0.3

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


#4-1
Raw0=0.3
n=200
p=50

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

d_re1 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[13,c(1,2,3)]<-apply(d_re1,2,mean)
Simulation_block_BH[13,4]<-mean(d_re1[,3]>c)

#4-2

p=50
n=100

time_start <- Sys.time()

d_re2 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[14,c(1,2,3)]<-apply(d_re2,2,mean)
Simulation_block_BH[14,4]<-mean(d_re2[,3]>c)
#4-3

n<-200
p<-100

time_start <- Sys.time()

d_re3 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[15,c(1,2,3)]<-apply(d_re3,2,mean)
Simulation_block_BH[15,4]<-mean(d_re3[,3]>c)

#4-4
n=100
p=100

time_start <- Sys.time()

d_re4 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_block_BH[16,c(1,2,3)]<-apply(d_re4,2,mean)
Simulation_block_BH[16,4]<-mean(d_re4[,3]>c)

#4-5

n=200
p=50
Raw0=0.25


time_start <- Sys.time()

d_re5 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[17,c(1,2,3)]<-apply(d_re5,2,mean)
Simulation_block_BH[17,4]<-mean(d_re5[,3]>c)

#4-6

n=100
p=50

time_start <- Sys.time()

d_re6 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_block_BH[18,c(1,2,3)]<-apply(d_re6,2,mean)
Simulation_block_BH[18,4]<-mean(d_re6[,3]>c)

#4-7

n=200
p=100

time_start <- Sys.time()

d_re7 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[19,c(1,2,3)]<-apply(d_re7,2,mean)
Simulation_block_BH[19,4]<-mean(d_re7[,3]>c)


#4-8

n=100
p=100

time_start <- Sys.time()

d_re8 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[20,c(1,2,3)]<-apply(d_re8,2,mean)
Simulation_block_BH[20,4]<-mean(d_re8[,3]>c)

#4-9

Raw0=0.2
n=200
p=50


time_start <- Sys.time()

d_re9 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_block_BH[21,c(1,2,3)]<-apply(d_re9,2,mean)
Simulation_block_BH[21,4]<-mean(d_re9[,3]>c)


#4-10

n=100
p=50


time_start <- Sys.time()

d_re10 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_block_BH[22,c(1,2,3)]<-apply(d_re10,2,mean)
Simulation_block_BH[22,4]<-mean(d_re10[,3]>c)


#4-11  
n=200
p=100

time_start <- Sys.time()

d_re11 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_block_BH[23,c(1,2,3)]<-apply(d_re11,2,mean)
Simulation_block_BH[23,4]<-mean(d_re11[,3]>c)

#4-12

n=100
p=100

time_start <- Sys.time()

d_re12 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_1(p,Raw0,q,Beta,gamma,n,m,s0,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh(W,X,D,p,n,m,q,B,aph,c,s0)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_block_BH[24,c(1,2,3)]<-apply(d_re12,2,mean)
Simulation_block_BH[24,4]<-mean(d_re12[,3]>c)












##Simulation5:proposed method+super diagonal+time dependencce raw=0

raw=0

#parameters
m=300
q=4
gamma=rep(0.4,q)
Beta=rep(0.5,q)
c=0.1
B=1000
aph=0.05
s0=6
s1=6

#5-1
Raw0=0.18
n=200
p=50

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

e_re1 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_proposed[1,c(1,2,3)]<-apply(e_re1,2,mean)
Simulation_off_proposed[1,4]<-mean(e_re1[,3]>c)

#5-2

p=50
n=100

time_start <- Sys.time()

e_re2 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[2,c(1,2,3)]<-apply(e_re2,2,mean)
Simulation_off_proposed[2,4]<-mean(e_re2[,3]>c)

#5-3

n<-200
p<-100

time_start <- Sys.time()

e_re3 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[3,c(1,2,3)]<-apply(e_re3,2,mean)
Simulation_off_proposed[3,4]<-mean(e_re3[,3]>c)

#5-4
n=100
p=100

time_start <- Sys.time()

e_re4 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[4,c(1,2,3)]<-apply(e_re4,2,mean)
Simulation_off_proposed[4,4]<-mean(e_re4[,3]>c)

#5-5

n=200
p=50
Raw0=0.20


time_start <- Sys.time()

e_re5 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[5,c(1,2,3)]<-apply(e_re5,2,mean)
Simulation_off_proposed[5,4]<-mean(e_re5[,3]>c)

#5-6

n=100
p=50

time_start <- Sys.time()

e_re6 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[6,c(1,2,3)]<-apply(e_re6,2,mean)
Simulation_off_proposed[6,4]<-mean(e_re6[,3]>c)


#5-7

n=200
p=100

time_start <- Sys.time()

e_re7 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_proposed[7,c(1,2,3)]<-apply(e_re7,2,mean)
Simulation_off_proposed[7,4]<-mean(e_re7[,3]>c)

#5-8

n=100
p=100

time_start <- Sys.time()

e_re8 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[8,c(1,2,3)]<-apply(e_re8,2,mean)
Simulation_off_proposed[8,4]<-mean(e_re8[,3]>c)

#5-9

Raw0=0.22
n=200
p=50


time_start <- Sys.time()

e_re9 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[9,c(1,2,3)]<-apply(e_re9,2,mean)
Simulation_off_proposed[9,4]<-mean(e_re9[,3]>c)


#5-10

n=100
p=50


time_start <- Sys.time()

e_re10 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[10,c(1,2,3)]<-apply(e_re10,2,mean)
Simulation_off_proposed[10,4]<-mean(e_re10[,3]>c)


#5-11  
n=200
p=100

time_start <- Sys.time()

e_re11 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_proposed[11,c(1,2,3)]<-apply(e_re11,2,mean)
Simulation_off_proposed[11,4]<-mean(e_re11[,3]>c)

#5-12

n=100
p=100

time_start <- Sys.time()

e_re12 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[12,c(1,2,3)]<-apply(e_re12,2,mean)
Simulation_off_proposed[12,4]<-mean(e_re12[,3]>c)

##Simulation6:proposed method+super diagonal+time dependencce raw=0.3


raw=0.3

#parameters
m=300
q=4
gamma=rep(0.4,q)
Beta=rep(0.5,q)
c=0.1
B=1000
aph=0.05
s0=6
s1=6

#6-1
Raw0=0.18
n=200
p=50

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

f_re1 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_proposed[13,c(1,2,3)]<-apply(f_re1,2,mean)
Simulation_off_proposed[13,4]<-mean(f_re1[,3]>c)

#6-2

p=50
n=100

time_start <- Sys.time()

f_re2 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_proposed[14,c(1,2,3)]<-apply(f_re2,2,mean)
Simulation_off_proposed[14,4]<-mean(f_re2[,3]>c)

#6-3

n<-200
p<-100

time_start <- Sys.time()

f_re3 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[15,c(1,2,3)]<-apply(f_re3,2,mean)
Simulation_off_proposed[15,4]<-mean(f_re3[,3]>c)

#6-4
n=100
p=100

time_start <- Sys.time()

f_re4 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[16,c(1,2,3)]<-apply(f_re4,2,mean)
Simulation_off_proposed[16,4]<-mean(f_re4[,3]>c)

#6-5

n=200
p=50
Raw0=0.20


time_start <- Sys.time()

f_re5 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[17,c(1,2,3)]<-apply(f_re5,2,mean)
Simulation_off_proposed[17,4]<-mean(f_re5[,3]>c)

#6-6

n=100
p=50

time_start <- Sys.time()

f_re6 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[18,c(1,2,3)]<-apply(f_re6,2,mean)
Simulation_off_proposed[18,4]<-mean(f_re6[,3]>c)


#6-7

n=200
p=100

time_start <- Sys.time()

f_re7 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[19,c(1,2,3)]<-apply(f_re7,2,mean)
Simulation_off_proposed[19,4]<-mean(f_re7[,3]>c)

#6-8

n=100
p=100

time_start <- Sys.time()

f_re8 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[20,c(1,2,3)]<-apply(f_re8,2,mean)
Simulation_off_proposed[20,4]<-mean(f_re8[,3]>c)

#6-9

Raw0=0.22
n=200
p=50


time_start <- Sys.time()

f_re9 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[21,c(1,2,3)]<-apply(f_re9,2,mean)
Simulation_off_proposed[21,4]<-mean(f_re9[,3]>c)

#6-10

n=100
p=50


time_start <- Sys.time()

f_re10 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[22,c(1,2,3)]<-apply(f_re10,2,mean)
Simulation_off_proposed[22,4]<-mean(f_re10[,3]>c)

#6-11  
n=200
p=100

time_start <- Sys.time()

f_re11 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[23,c(1,2,3)]<-apply(f_re11,2,mean)
Simulation_off_proposed[23,4]<-mean(f_re11[,3]>c)

#6-12

n=100
p=100

time_start <- Sys.time()

f_re12 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_diag_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_proposed[24,c(1,2,3)]<-apply(f_re12,2,mean)
Simulation_off_proposed[24,4]<-mean(f_re12[,3]>c)

##Simulation7:BH method+super diagonal+time dependencce raw=0


raw=0

#parameters
m=300
q=4
gamma=rep(0.4,q)
Beta=rep(0.5,q)
c=0.1
B=1000
aph=0.05
s0=6
s1=6

#7-1
Raw0=0.18
n=200
p=50

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

g_re1 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_BH[1,c(1,2,3)]<-apply(g_re1,2,mean)
Simulation_off_BH[1,4]<-mean(g_re1[,3]>c)

#7-2

p=50
n=100

time_start <- Sys.time()

g_re2 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_BH[2,c(1,2,3)]<-apply(g_re2,2,mean)
Simulation_off_BH[2,4]<-mean(g_re2[,3]>c)

#7-3

n<-200
p<-100

time_start <- Sys.time()

g_re3 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_BH[3,c(1,2,3)]<-apply(g_re3,2,mean)
Simulation_off_BH[3,4]<-mean(g_re3[,3]>c)


#7-4
n=100
p=100

time_start <- Sys.time()

g_re4 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[4,c(1,2,3)]<-apply(g_re4,2,mean)
Simulation_off_BH[4,4]<-mean(g_re4[,3]>c)

#7-5

n=200
p=50
Raw0=0.20


time_start <- Sys.time()

g_re5 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[5,c(1,2,3)]<-apply(g_re5,2,mean)
Simulation_off_BH[5,4]<-mean(g_re5[,3]>c)

#7-6

n=100
p=50

time_start <- Sys.time()

g_re6 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[6,c(1,2,3)]<-apply(g_re6,2,mean)
Simulation_off_BH[6,4]<-mean(g_re6[,3]>c)


#7-7

n=200
p=100

time_start <- Sys.time()

g_re7 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[7,c(1,2,3)]<-apply(g_re7,2,mean)
Simulation_off_BH[7,4]<-mean(g_re7[,3]>c)


#7-8

n=100
p=100

time_start <- Sys.time()

g_re8 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[8,c(1,2,3)]<-apply(g_re8,2,mean)
Simulation_off_BH[8,4]<-mean(g_re8[,3]>c)

#7-9

Raw0=0.22
n=200
p=50


time_start <- Sys.time()

g_re9 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[9,c(1,2,3)]<-apply(g_re9,2,mean)
Simulation_off_BH[9,4]<-mean(g_re9[,3]>c)


#7-10

n=100
p=50


time_start <- Sys.time()

g_re10 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[10,c(1,2,3)]<-apply(g_re10,2,mean)
Simulation_off_BH[10,4]<-mean(g_re10[,3]>c)

#7-11  
n=200
p=100

time_start <- Sys.time()

g_re11 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_BH[11,c(1,2,3)]<-apply(g_re11,2,mean)
Simulation_off_BH[11,4]<-mean(g_re11[,3]>c)



#7-12

n=100
p=100

time_start <- Sys.time()

g_re12 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_BH[12,c(1,2,3)]<-apply(g_re12,2,mean)
Simulation_off_BH[12,4]<-mean(g_re12[,3]>c)


##Simulation8:BH method+super diagonal+time dependencce raw=0.3


raw=0.3

#parameters
m=300
q=4
gamma=rep(0.4,q)
Beta=rep(0.5,q)
c=0.1
B=1000
aph=0.05
s0=6
s1=6

#8-1
Raw0=0.18
n=200
p=50

#timestart is used to calculate the time need for a simulation, can be omitted
time_start <- Sys.time()

h_re1 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_BH[13,c(1,2,3)]<-apply(h_re1,2,mean)
Simulation_off_BH[13,4]<-mean(h_re1[,3]>c)

#8-2

p=50
n=100

time_start <- Sys.time()

h_re2 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_BH[14,c(1,2,3)]<-apply(h_re2,2,mean)
Simulation_off_BH[14,4]<-mean(h_re2[,3]>c)

#8-3

n<-200
p<-100

time_start <- Sys.time()

h_re3 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[15,c(1,2,3)]<-apply(h_re3,2,mean)
Simulation_off_BH[15,4]<-mean(h_re3[,3]>c)

#8-4
n=100
p=100

time_start <- Sys.time()

h_re4 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[16,c(1,2,3)]<-apply(h_re4,2,mean)
Simulation_off_BH[16,4]<-mean(h_re4[,3]>c)

#8-5

n=200
p=50
Raw0=0.20


time_start <- Sys.time()

h_re5 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_BH[17,c(1,2,3)]<-apply(h_re5,2,mean)
Simulation_off_BH[17,4]<-mean(h_re5[,3]>c)


#8-6

n=100
p=50

time_start <- Sys.time()

h_re6 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[18,c(1,2,3)]<-apply(h_re6,2,mean)
Simulation_off_BH[18,4]<-mean(h_re6[,3]>c)


#8-7

n=200
p=100

time_start <- Sys.time()

h_re7 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)

Simulation_off_BH[19,c(1,2,3)]<-apply(h_re7,2,mean)
Simulation_off_BH[19,4]<-mean(h_re7[,3]>c)



#8-8

n=100
p=100

time_start <- Sys.time()

h_re8 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_2(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[20,c(1,2,3)]<-apply(h_re8,2,mean)
Simulation_off_BH[20,4]<-mean(h_re8[,3]>c)

#8-9

Raw0=0.22
n=200
p=50


time_start <- Sys.time()

h_re9 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[21,c(1,2,3)]<-apply(h_re9,2,mean)
Simulation_off_BH[21,4]<-mean(h_re9[,3]>c)


#8-10

n=100
p=50


time_start <- Sys.time()

h_re10 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[20,c(1,2,3)]<-apply(h_re10,2,mean)
Simulation_off_BH[20,4]<-mean(h_re10[,3]>c)

#8-11  
n=200
p=100

time_start <- Sys.time()

h_re11 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[23,c(1,2,3)]<-apply(h_re11,2,mean)
Simulation_off_BH[23,4]<-mean(h_re11[,3]>c)


#8-12

n=100
p=100

time_start <- Sys.time()

h_re12 <- foreach (
  k = 1:1000,
  .combine = rbind,
  .inorder = TRUE,
  .export = c("re")
) %dopar% {
  library(MASS)
  DATA<-Sample_3(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw)
  W<-DATA[[1]]
  D<-DATA[[2]]
  X<-DATA[[3]]
  re<-corr_estimate_bh_2(W,X,D,p,n,m,q,B,aph,c,s0,s1)
}

time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)


Simulation_off_BH[24,c(1,2,3)]<-apply(h_re12,2,mean)
Simulation_off_BH[24,4]<-mean(h_re12[,3]>c)










