rm(list=ls())
setwd("")
library(tensorflow)
use_python("")
source("EBtimecourse.R")
library(mvtnorm)

set.seed(1)

P=0.2
N=5000;

nDE=P*N; nEE=N-nDE;
cp_real_seq_index=1:nDE;
timePoint=8; replicate=3;
Ti=timePoint*replicate; 

changePointTable = data.frame(matrix(NA, nrow=(timePoint-1)+(timePoint-1)*(timePoint-2)/2, ncol=3), stringsAsFactors=F)
colnames(changePointTable) = c("n1", "n2", "n3")
changePointTable[1:(timePoint-1),"n1"]=1:(timePoint-1)
changePointTable[1:(timePoint-1),"n2"]=(timePoint-1):1
changePointTable[1:(timePoint-1),"n3"]=0
combT = as.data.frame(t(combn(timePoint-1, 2)))
combT$V2 = combT$V2 - combT$V1
changePointTable[timePoint:nrow(changePointTable),c("n1","n2")]=combT
changePointTable[timePoint:nrow(changePointTable),"n3"]=timePoint-rowSums(changePointTable[timePoint:nrow(changePointTable),c("n1", "n2")])
combNumber = nrow(changePointTable)

ss=sample(1:nrow(changePointTable), P*N, replace=T)
n1=changePointTable[ss,"n1"]
n2=changePointTable[ss,"n2"]
n3=changePointTable[ss,"n3"]

n1=c(n1, rep(timePoint,nEE)); n2=c(n2, rep(0, nEE)); n3=c(n3, rep(0, nEE));

nu=0; tau=3;
mu=matrix(rnorm(N*2,mean=nu,sd=tau),nc=2);
parameter=data.frame(mu, n1=n1, n2=n2, n3=n3);
colnames(parameter)=c("mu1", "mu2", "n1", "n2", "n3");

mean.de.list = apply(parameter[1:nDE,], 1, function(x) 
  {c(rep(x[1], x[3]*replicate), rep(x[2], x[4]*replicate), rep(x[1], x[5]*replicate))})
mean.de = matrix(unlist(mean.de.list), byrow=T, ncol=Ti)

mean.ee.list = apply(parameter[(nDE+1):N,], 1, function(x) {rep(x[1], x[3]*replicate)})
mean.ee = matrix(unlist(mean.ee.list), byrow=T, ncol=Ti)

mean.matrix <- rbind(mean.de, mean.ee)

precision.block <- matrix(0, nrow=30, ncol=30)
for(i in 1:30) {
  precision.block[i, i] <- 1
  if(i > 1) precision.block[i-1, i] <- 0.2
  if(i < 30) precision.block[i+1, i] <- 0.2
}
sigma.block <- solve(precision.block)
sigma.matrix <- matrix(0, nrow=N, ncol=N)
for(i in 1:floor(N/30)) {
  sigma.matrix[((i-1)*30+1):(i*30), ((i-1)*30+1):(i*30)] <- sigma.block
}
residual.index <- N-floor(N/30)*30
sigma.matrix[(floor(N/30)*30+1):N, (floor(N/30)*30+1):N] <- sigma.block[1:residual.index, 1:residual.index]

for(i in 1:Ti) {
  if(!exists("gene.exp")) {
    gene.exp <- matrix(rmvnorm(1, mean=mean.matrix[,i], sigma=sigma.matrix), nrow=N, ncol=1)
  } else {
    gene.exp <- cbind(gene.exp, matrix(rmvnorm(1, mean=mean.matrix[,i], sigma=sigma.matrix), nrow=N, ncol=1))
  }
}

ptm <- proc.time()
result = EBtimecourse(exp.dat = gene.exp, timepoint = timePoint, replicate = replicate, FDR=0.1, verbose=T)
print(proc.time() - ptm)

#### test the power and FDR of NNG model #####
cp_power = sum(result$cp.index %in% cp_real_seq_index)/nDE
cp_FDR = sum(!(result$cp.index %in% cp_real_seq_index))/length(result$cp.index)
cp_position_correct_numer = 0
cp_position_not_correct_numer = 0
for(i in 1:nrow(result$cp.position)) {
  if(parameter[result$cp.position$SeqID[i],"n1"]==result$cp.position[i, "n1"] & parameter[result$cp.position$SeqID[i],"n2"]==result$cp.position[i, "n2"])
    cp_position_correct_numer=cp_position_correct_numer+1
  else
    cp_position_not_correct_numer = cp_position_not_correct_numer+1
}
cp_position_power=cp_position_correct_numer/nDE
cp_position_FDR=cp_position_not_correct_numer/nrow(result$cp.position)

