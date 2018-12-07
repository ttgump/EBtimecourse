rm(list=ls())
source("EBtimecourse.R")

P=0.1
N=1000;

nDE=P*N; nEE=N-nDE;
cp_real_seq_index=1:nDE;
timePoint=10; replicate=3;
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
mu0=0; kappa0=0.01; alpha0=10; beta0=2;
lambda = rgamma(N*2, shape=alpha0, rate=beta0)
mu = rnorm(N*2,mean=mu0,sd=1/sqrt(kappa0*lambda))
lambda = matrix(lambda, ncol=2)
mu = matrix(mu, ncol=2)
sigma = 1/sqrt(lambda)
parameter=data.frame(mu, sigma, n1=n1, n2=n2, n3=n3)
colnames(parameter)=c("mu1", "mu2", "sigma1", "sigma2", "n1", "n2", "n3")

gene.de=t(apply(parameter[1:nDE,],1, function(x) 
  c(rnorm(x[5]*replicate, mean=x[1], sd=x[3]), rnorm(x[6]*replicate, mean=x[2], sd=x[4]),
    rnorm(x[7]*replicate, mean=x[1], sd=x[3]))))
gene.ee=t(apply(parameter[(nDE+1):N,],1, function(x) c(rnorm(x[5]*replicate, mean=x[1], sd=x[3]))))

gene.exp=rbind(gene.de, gene.ee)

result = EBtimecourse(exp.dat = gene.exp, timepoint = timePoint, replicate = replicate, model = "NNG", FDR = 0.1)

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
