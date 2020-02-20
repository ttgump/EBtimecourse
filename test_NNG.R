rm(list=ls())
setwd("~/Documents/Research/Microarray_time_course/manuscript/EBtimecourse_tf")
library(tensorflow)
use_python("/usr/local/bin/python3")
source("EBtimecourse.R")

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
mu0=0; kappa0=0.1; alpha0=1; beta0=10;
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


#### test baseline methods
FDR = 0.1
library(timecourse)
library(gprege)
library(limma)

# Compare with Yu Chuan Tai's method
assay = rep(c(LETTERS[1:replicate]), timePoint)
time.grp = rep(c(1:timePoint), each=replicate)
size = rep(replicate, N)
out = mb.long(gene.exp, times=timePoint, reps=size, rep.grp=assay, time.grp=time.grp)
cp_FDR_Tai=0;
cp_power_Tai=0;
for(i in 1:N) {
  cp_estimated_sq_index_Tai = out$pos.HotellingT2[1:i];
  cp_power_Tai = sum(cp_estimated_sq_index_Tai %in% cp_real_seq_index)/nDE;
  cp_FDR_Tai = sum(!(cp_estimated_sq_index_Tai %in% cp_real_seq_index))/length(cp_estimated_sq_index_Tai);
  if(cp_FDR_Tai >= cp_FDR) {break;}
}

# Compare with GP method
tTrue <- matrix(rep(1:timePoint, each = replicate), ncol=1)
gpregeOptions <- list()
gpregeOptions$indexRange <- 1:nrow(gene.exp)

# Explore individual profiles in interactive mode.
gpregeOptions$explore <- FALSE
# Exhaustive plot resolution of the LML function.
gpregeOptions$exhaustPlotRes <- 30
# Exhaustive plot contour levels.
gpregeOptions$exhaustPlotLevels <- 10
# Exhaustive plot maximum lengthscale.
gpregeOptions$exhaustPlotMaxWidth <- 100
# Noisy ground truth labels: which genes are in the top 786 ranks of the TSNI ranking.
gpregeOptions$labels <- rownames(gene.exp)
# SCG optimisation: maximum number of iterations.
gpregeOptions$iters <- 100
# SCG optimisation: no messages.
gpregeOptions$display <- FALSE
# Matrix of different hyperparameter configurations as rows:
# [inverse-lengthscale   percent-signal-variance   percent-noise-variance].
gpregeOptions$inithypers <- matrix( c(  1/1000, 1e-3, 0.999,
                                        1/30, 0.999, 1e-3,
                                        1/80, 2/3, 1/3), ncol=3, byrow=TRUE)


gpregeOutput<-gprege(data=gene.exp,inputs=tTrue,gpregeOptions=gpregeOptions)
gpregeOrder <- order(gpregeOutput$rankingScores, decreasing = T)
cp_FDR_GP=0;
cp_power_GP=0;
for(i in 1:N) {
  cp_estimated_sq_index_GP = gpregeOrder[1:i];
  cp_power_GP = sum(cp_estimated_sq_index_GP %in% cp_real_seq_index)/nDE;
  cp_FDR_GP = sum(!(cp_estimated_sq_index_GP %in% cp_real_seq_index))/length(cp_estimated_sq_index_GP);
  if(cp_FDR_GP >= cp_FDR) {break;}
}

#### limma
sample <- data.frame(Sample = 1:Ti,
                     Time = rep(1:timePoint, each = replicate))
sample$Time <- factor(sample$Time, levels = 1:timePoint)
design <- model.matrix(~sample$Time)
fit <- lmFit(gene.exp, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef=2:timePoint, adjust.method="BH", number = nrow(gene.exp))
s.tt <- tt[tt$adj.P.Val <= FDR, ]

limma_cp_estimated_sq_index <- as.numeric(rownames(s.tt))

cp_power_limma = sum(limma_cp_estimated_sq_index %in% cp_real_seq_index)/nDE;
cp_FDR_limma = sum(!(limma_cp_estimated_sq_index %in% cp_real_seq_index))/length(limma_cp_estimated_sq_index);

# Frequentist method to get change point position
cp_position_Tai = c()
for(i in 1:length(cp_estimated_sq_index_Tai)) {
  expSeq = gene.exp[cp_estimated_sq_index_Tai[i],]
  cp.fr.pvalue = c()
  for(j in 1:combNumber) {
    if((changePointTable[j,"n1"]+changePointTable[j,"n2"])<timePoint) {
      seq1 = expSeq[c(1:(changePointTable[j,"n1"]*replicate), 
                      ((changePointTable[j,"n1"]+changePointTable[j,"n2"])*replicate+1):Ti)]
      seq2 = expSeq[(changePointTable[j,"n1"]*replicate+1):((changePointTable[j,"n1"]+changePointTable[j,"n2"])*replicate)]
    } else {
      seq1 = expSeq[1:(changePointTable[j,"n1"]*replicate)]
      seq2 = expSeq[(changePointTable[j,"n1"]*replicate+1):Ti]
    }
    
    cp.fr.pvalue = c(cp.fr.pvalue, t.test(seq1, seq2)$p.value)
  }
  cp_position_Tai = c(cp_position_Tai, which.min(cp.fr.pvalue))
}
cp.position.Tai = data.frame(SeqID=cp_estimated_sq_index_Tai, changePointTable[cp_position_Tai,], stringsAsFactors=F)
cp_position_Tai_correct_numer = 0;
cp_position_Tai_not_correct_numer = 0;
for(i in 1:nrow(cp.position.Tai)) {
  if(parameter[cp.position.Tai$SeqID[i],"n1"]==cp.position.Tai[i, "n1"] & parameter[cp.position.Tai$SeqID[i],"n2"]==cp.position.Tai[i, "n2"])
    cp_position_Tai_correct_numer=cp_position_Tai_correct_numer+1
  else
    cp_position_Tai_not_correct_numer = cp_position_Tai_not_correct_numer+1
}
cp_position_power_Tai=cp_position_Tai_correct_numer/nDE;
cp_position_FDR_Tai=cp_position_Tai_not_correct_numer/nrow(cp.position.Tai);

# Frequentist method to get change point position
cp_position_GP = c()
for(i in 1:length(cp_estimated_sq_index_GP)) {
  expSeq = gene.exp[cp_estimated_sq_index_GP[i],]
  cp.fr.pvalue = c()
  for(j in 1:combNumber) {
    if((changePointTable[j,"n1"]+changePointTable[j,"n2"])<timePoint) {
      seq1 = expSeq[c(1:(changePointTable[j,"n1"]*replicate), 
                      ((changePointTable[j,"n1"]+changePointTable[j,"n2"])*replicate+1):Ti)]
      seq2 = expSeq[(changePointTable[j,"n1"]*replicate+1):((changePointTable[j,"n1"]+changePointTable[j,"n2"])*replicate)]
    } else {
      seq1 = expSeq[1:(changePointTable[j,"n1"]*replicate)]
      seq2 = expSeq[(changePointTable[j,"n1"]*replicate+1):Ti]
    }
    
    cp.fr.pvalue = c(cp.fr.pvalue, t.test(seq1, seq2)$p.value)
  }
  cp_position_GP = c(cp_position_GP, which.min(cp.fr.pvalue))
}
cp.position.GP = data.frame(SeqID=cp_estimated_sq_index_GP, changePointTable[cp_position_GP,], stringsAsFactors=F)
cp_position_GP_correct_numer = 0;
cp_position_GP_not_correct_numer = 0;
for(i in 1:nrow(cp.position.GP)) {
  if(parameter[cp.position.GP$SeqID[i],"n1"]==cp.position.GP[i, "n1"] & parameter[cp.position.GP$SeqID[i],"n2"]==cp.position.GP[i, "n2"])
    cp_position_GP_correct_numer=cp_position_GP_correct_numer+1
  else
    cp_position_GP_not_correct_numer = cp_position_GP_not_correct_numer+1
}
cp_position_power_GP=cp_position_GP_correct_numer/nDE;
cp_position_FDR_GP=cp_position_GP_not_correct_numer/nrow(cp.position.GP);

# Frequentist method to get change point position
cp_position_limma = c()
for(i in 1:length(limma_cp_estimated_sq_index)) {
  expSeq = gene.exp[limma_cp_estimated_sq_index[i],]
  cp.fr.pvalue = c()
  for(j in 1:combNumber) {
    if((changePointTable[j,"n1"]+changePointTable[j,"n2"])<timePoint) {
      seq1 = expSeq[c(1:(changePointTable[j,"n1"]*replicate), 
                      ((changePointTable[j,"n1"]+changePointTable[j,"n2"])*replicate+1):Ti)]
      seq2 = expSeq[(changePointTable[j,"n1"]*replicate+1):((changePointTable[j,"n1"]+changePointTable[j,"n2"])*replicate)]
    } else {
      seq1 = expSeq[1:(changePointTable[j,"n1"]*replicate)]
      seq2 = expSeq[(changePointTable[j,"n1"]*replicate+1):Ti]
    }
    
    cp.fr.pvalue = c(cp.fr.pvalue, t.test(seq1, seq2)$p.value)
  }
  cp_position_limma = c(cp_position_limma, which.min(cp.fr.pvalue))
}
cp.position.limma = data.frame(SeqID=limma_cp_estimated_sq_index, changePointTable[cp_position_limma,], stringsAsFactors=F)
cp_position_limma_correct_numer = 0;
cp_position_limma_not_correct_numer = 0;
for(i in 1:nrow(cp.position.limma)) {
  if(parameter[cp.position.limma$SeqID[i],"n1"]==cp.position.limma[i, "n1"] & parameter[cp.position.limma$SeqID[i],"n2"]==cp.position.limma[i, "n2"])
    cp_position_limma_correct_numer=cp_position_limma_correct_numer+1
  else
    cp_position_limma_not_correct_numer = cp_position_limma_not_correct_numer+1
}
cp_position_power_limma=cp_position_limma_correct_numer/nDE;
cp_position_FDR_limma=cp_position_limma_not_correct_numer/nrow(cp.position.limma);

# Output to file
FDR.list = list(Prob=nDE/N, FDR=FDR, cp.power=cp_power, cp.FDR=cp_FDR, 
                cp.power.Tai=cp_power_Tai, cp.FDR.Tai=cp_FDR_Tai,
                cp.power.GP=cp_power_GP, cp.FDR.GP=cp_FDR_GP,
                cp.power.limma=cp_power_limma, cp.FDR.limma=cp_FDR_limma,
                cp.position.power=cp_position_power, cp.position.FDR=cp_position_FDR,
                cp.position.Tai.power=cp_position_power_Tai, cp.position.Tai.FDR=cp_position_FDR_Tai,
                cp.position.GP.power=cp_position_power_GP, cp.position.GP.FDR=cp_position_FDR_GP,
                cp.position.limma.power=cp_position_power_limma, cp.position.limma.FDR=cp_position_FDR_limma,
                convergence=fit$convcode);

FDR.list.df = do.call("rbind", lapply(FDR.list, as.data.frame))
write.csv(FDR.list.df, paste(paste("P",P,sep=""), indx, "FDR.csv", sep="_"))
write.csv(gene.exp, paste(paste("P",P,sep=""), indx, "gene_exp.csv", sep="_"))
save(result, file=paste(paste("P",P,sep=""), indx, "fit.RData", sep="_"))
write.csv(cp_estimated_sq_index_Tai, paste(paste("P",P,sep=""), indx, "Tai_cp_index.csv", sep="_"))
write.csv(cp_estimated_sq_index_GP, paste(paste("P",P,sep=""), indx, "GP_cp_index.csv", sep="_"))
write.csv(limma_cp_estimated_sq_index, paste(paste("P",P,sep=""), indx, "limma_cp_index.csv", sep="_"))
write.csv(cp.position.Tai, paste(paste("P",P,sep=""), indx, "Tai_cp_position.csv", sep="_"))
write.csv(cp.position.GP, paste(paste("P",P,sep=""), indx, "GP_cp_position.csv", sep="_"))
write.csv(cp.position.limma, paste(paste("P",P,sep=""), indx, "limma_cp_position.csv", sep="_"))