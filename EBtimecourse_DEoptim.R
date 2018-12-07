require(DEoptim)
require(foreach)
library(doMC)
registerDoMC()

### auxiliary functions
NN_logL_ij <- function(P, sigma, nu, tau, i, cp, changePointTable, combNumber, Ti, left_seq_index_table, left_mean_row, right_mean_row, left_square_rowsum) {
  if(cp==0) {
    mean_xi = left_mean_row[i, combNumber+1]
    log_pr = log(1-P) + log(sigma) - Ti*(0.5*log(2*pi)+log(sigma)) - 0.5*log(Ti*(tau^2)+sigma^2) +
      (-left_square_rowsum[i, combNumber+1]/2/(sigma^2) - nu^2/2/(tau^2)) +
      (tau^2*Ti^2*mean_xi^2/(sigma^2) + sigma^2*nu^2/(tau^2)+2*Ti*mean_xi*nu)/(2*(Ti*tau^2+sigma^2));
  }else {
    n1 = length(left_seq_index_table[[cp]])
    n2 = Ti - length(left_seq_index_table[[cp]])
    
    mean_xi1 = left_mean_row[i, cp]
    mean_xi2 = right_mean_row[i, cp]

    log_pr = log(P) - log(combNumber) +
      2*log(sigma) - Ti*(0.5*log(2*pi)+log(sigma)) -0.5*log(n1*(tau^2)+sigma^2) - 0.5*log(n2*(tau^2)+sigma^2) +
      (-left_square_rowsum[i, combNumber+1]/2/(sigma^2) - nu^2/(tau^2)) +
      (tau^2*n1^2*mean_xi1^2/(sigma^2)+sigma^2*nu^2/(tau^2)+2*n1*mean_xi1*nu)/(2*(n1*tau^2+sigma^2)) +
      (tau^2*n2^2*mean_xi2^2/(sigma^2)+sigma^2*nu^2/(tau^2)+2*n2*mean_xi2*nu)/(2*(n2*tau^2+sigma^2));
  }
  return(log_pr);
}

NN_logL_i <- function(P, sigma, m, tau, changePointTable, combNumber, N, Ti, left_seq_index_table, left_mean_row, right_mean_row, left_square_rowsum) {
  r = 0;
  
  for(i in 1:N) {
    prob_i_log_v = c();
    for(j in 0:combNumber) {
      prob_i_log_v = c(prob_i_log_v, NN_logL_ij(P, sigma, nu, tau, i, j, changePointTable, combNumber, Ti, left_seq_index_table, left_mean_row, right_mean_row, left_square_rowsum));
    }
    prob_i_log = max(prob_i_log_v) + log1p(sum(exp(prob_i_log_v[-which.max(prob_i_log_v)] - max(prob_i_log_v))));
    r = r + prob_i_log;
  }
  return(r);
}

# object function for normal normal model
NN.obj = function(parameter, changePointTable, combNumber, N, Ti, left_seq_index_table, left_mean_row, right_mean_row, left_square_rowsum) {
  #P, sigma, mu, tau
  r = NN_logL_i(parameter[1], parameter[2], parameter[3], parameter[4],
                changePointTable, combNumber, N, Ti, left_seq_index_table, left_mean_row, right_mean_row, left_square_rowsum)
  r = -r
  return(r)
}

NNG_logL_ij <- function(P, mu, kappa, alpha, beta, i, cp, changePointTable, combNumber, Ti, left_seq_index_table, left_mean_row, right_mean_row, left_squareDeviation, right_squareDeviation) {
  if(cp==0) {
    mean_xi = left_mean_row[i, combNumber+1]
    square_deviation_xi = left_squareDeviation[i, combNumber+1]
    
    kappa_n = kappa + Ti
    alpha_n = alpha + Ti/2
    beta_n = beta + 0.5*square_deviation_xi + (kappa*Ti*(mean_xi-mu)^2)/2/(kappa+Ti)
    
    log_pr = log(1-P) + lgamma(alpha_n) - lgamma(alpha) + alpha*log(beta) - alpha_n*log(beta_n) +
      0.5*(log(kappa) - log(kappa_n)) - Ti/2*log(2*pi)
  }else {
    mean_xi_1 = left_mean_row[i, cp]
    square_deviation_xi_1 = left_squareDeviation[i, cp]
    mean_xi_2 = right_mean_row[i, cp]
    square_deviation_xi_2 = right_squareDeviation[i, cp]
    n1 = length(left_seq_index_table[[cp]])
    n2 = Ti - length(left_seq_index_table[[cp]])
    
    kappa_n_1 = kappa + n1
    alpha_n_1 = alpha + n1/2
    beta_n_1 = beta + 0.5*square_deviation_xi_1 + (kappa*n1*(mean_xi_1-mu)^2)/2/(kappa+n1)
    kappa_n_2 = kappa + n2
    alpha_n_2 = alpha + n2/2
    beta_n_2 = beta + 0.5*square_deviation_xi_2 + (kappa*n2*(mean_xi_2-mu)^2)/2/(kappa+n2)
    
    log_pr = log(P) - log(combNumber) +
      lgamma(alpha_n_1) - alpha_n_1*log(beta_n_1) - 0.5*log(kappa_n_1) +
      lgamma(alpha_n_2) - alpha_n_2*log(beta_n_2) - 0.5*log(kappa_n_2) - 
      Ti/2*log(2*pi) - 2*lgamma(alpha) + 2*alpha*log(beta) + log(kappa)
  }
  return(log_pr);
}

NNG_logL_i <- function(P, mu, kappa, alpha, beta, changePointTable, combNumber, N, Ti, left_seq_index_table, left_mean_row, right_mean_row, left_squareDeviation, right_squareDeviation) {
  r = 0;  
  for(i in 1:N) {
    prob_i_log_v = c();
    for(cp in 0:combNumber) {
      prob_i_log_v = c(prob_i_log_v, NNG_logL_ij(P, mu, kappa, alpha, beta, i, cp, changePointTable, combNumber, Ti, left_seq_index_table, left_mean_row, right_mean_row, left_squareDeviation, right_squareDeviation));
    }
    prob_i_log = max(prob_i_log_v) + log1p(sum(exp(prob_i_log_v[-which.max(prob_i_log_v)] - max(prob_i_log_v))));
    r = r + prob_i_log;
  }
  return(r);
}

# object function for normal normal gamma model
NNG.obj = function(parameter, changePointTable, combNumber, N, Ti, left_seq_index_table, left_mean_row, right_mean_row, left_squareDeviation, right_squareDeviation) {
  #prob, mu0, kappa0, alpha0, beta0
  r = NNG_logL_i(parameter[1], parameter[2], parameter[3], parameter[4], parameter[5],
                 changePointTable, combNumber, N, Ti, left_seq_index_table, left_mean_row, right_mean_row, left_squareDeviation, right_squareDeviation)
  r = -r
  return(r)
}

my_mean = function(x) apply(x, 1, mean)
my_sumsquare = function(x) apply(x, 1, function(x) sum(x^2))
my_squareDeviation = function(x) apply(x, 1, function(x) sum((x-mean(x))^2))

# Empirical Bayes model to identify differentially expressed genes in transcriptome timecourse data
# Tian Tian
# Created 23 Jun 2017
# Parameters:
##  exp.dat - timecourse gene expression data matrix, with rows are genes and columns are time points
##  timepoint - number of timepoints
##  replicate - number of replicated (columns of data matrix = timepoint*replicate)
##  model - one of "NN" or "NNG" for normal normal modle or normal normal gamma model
##  FDR - FDR level for detecting change points
#  output:
##  cp.index - gene indices that detected to have change points
##  cp.position - change point positions
##  estimated.parameter - hyperparameter estimation
EBtimecourse = function(exp.dat=NA, timepoint=NA, replicate=NA, model="NN", FDR=0.1) {
  stopifnot(model %in% c("NN", "NNG"))
  stopifnot(class(exp.dat) == "matrix")
  stopifnot(class(timepoint) == "numeric" & class(replicate) == "numeric" & class(FDR) == "numeric")
  stopifnot(ncol(exp.dat) == timepoint*replicate)
  
  # calculate global values
  N = nrow(exp.dat)
  changePointTable = data.frame(matrix(NA, nrow=(timepoint-1)+(timepoint-1)*(timepoint-2)/2, ncol=3), stringsAsFactors=F)
  colnames(changePointTable) = c("n1", "n2", "n3")
  changePointTable[1:(timepoint-1),"n1"] = 1:(timepoint-1)
  changePointTable[1:(timepoint-1),"n2"] = (timepoint-1):1
  changePointTable[1:(timepoint-1),"n3"] = 0
  combT = as.data.frame(t(combn(timepoint-1, 2)))
  combT$V2 = combT$V2 - combT$V1
  changePointTable[timepoint:nrow(changePointTable),c("n1","n2")] = combT
  changePointTable[timepoint:nrow(changePointTable),"n3"] = timepoint-rowSums(changePointTable[timepoint:nrow(changePointTable),c("n1", "n2")])
  combNumber = nrow(changePointTable)
  Ti = timepoint*replicate
  left_seq_index_table = foreach(cp = 1:(combNumber+1)) %do% {
    if(cp==combNumber+1) 
      return(1:Ti)
    else {
      rho1 = changePointTable[cp,"n1"]*replicate
      rho2 = (changePointTable[cp,"n1"]+changePointTable[cp,"n2"])*replicate
      
      if(rho2<Ti)
        return(c(1:rho1, (rho2+1):Ti))
      else
        return(c(1:rho1))
    }
  }
  
  left_submatrix = lapply(left_seq_index_table, function(i) exp.dat[, i, drop=F])
  right_submatrix = lapply(left_seq_index_table, function(i) exp.dat[, -i, drop=F])
  
  left_mean_row = sapply(left_submatrix, my_mean)
  right_mean_row = sapply(right_submatrix, my_mean)
  
  left_square_rowsum = sapply(left_submatrix, my_sumsquare)
  right_square_rowsum = sapply(right_submatrix, my_sumsquare)
  
  left_squareDeviation = sapply(left_submatrix, my_squareDeviation)
  right_squareDeviation = sapply(right_submatrix, my_squareDeviation)
  
  print("Optimizing hyperparameters...")
  if(model == "NN") {
    fit = DEoptim(lower=c(1e-9,1e-9,-100,1e-9), 
                 upper=c(0.999999999,100,100,100),
                 fn=NN.obj, changePointTable=changePointTable, combNumber=combNumber, 
                 N=N, Ti=Ti, left_seq_index_table=left_seq_index_table, 
                 left_mean_row=left_mean_row, right_mean_row=right_mean_row, 
                 left_square_rowsum=left_square_rowsum,
                 DEoptim.control(parallelType=2))
    
    P_hat = fit$optim$bestmem[1]; sigma_hat = fit$optim$bestmem[2]; 
    nu_hat = fit$optim$bestmem[3]; tau_hat = fit$optim$bestmem[4];
    
    parameter_hat = list(P_hat = P_hat, sigma_hat = sigma_hat, nu_hat = nu_hat, tau_hat = tau_hat)
    log_likelihood = matrix(data=rep(0, N*(combNumber+1)), nrow=N, ncol=(combNumber+1))
    for(i in 1:N) {
      for(j in 0:combNumber) {
        log_likelihood[i,j+1] = NN_logL_ij(P_hat, sigma_hat, nu_hat, tau_hat, i, j,
                                           changePointTable, combNumber, Ti, left_seq_index_table, 
                                           left_mean_row, right_mean_row, left_square_rowsum)
      }
    }
    likelihood = exp(log_likelihood)
    prob = likelihood/rowSums(likelihood)
  }
  else {
    fit = DEoptim(lower=c(1e-9,-100,1e-9,1e-9,1e-9), 
                 upper=c(0.999999999,100,100,100,100), 
                 fn=NNG.obj, changePointTable=changePointTable, combNumber=combNumber, 
                 N=N, Ti=Ti, left_seq_index_table=left_seq_index_table, 
                 left_mean_row=left_mean_row, right_mean_row=right_mean_row, 
                 left_squareDeviation=left_squareDeviation, right_squareDeviation=right_squareDeviation,
                 DEoptim.control(parallelType=2))
    
    P_hat = fit$optim$bestmem[1]; mu0_hat = fit$optim$bestmem[2]; 
    kappa0_hat = fit$optim$bestmem[3]; alpha0_hat = fit$optim$bestmem[4]; 
    beta0_hat = fit$optim$bestmem[5];
    
    parameter_hat = list(P_hat = P_hat, mu0_hat = mu0_hat, kappa0_hat = kappa0_hat, alpha0_hat = alpha0_hat, beta0_hat = beta0_hat)
    log_likelihood = matrix(data=rep(0, N*(combNumber+1)), nrow=N, ncol=(combNumber+1))
    for(i in 1:N) {
      for(j in 0:combNumber) {
        log_likelihood[i,j+1] = NNG_logL_ij(P_hat, mu0_hat, kappa0_hat, alpha0_hat, beta0_hat, i, j,
                                            changePointTable, combNumber, Ti, left_seq_index_table, 
                                            left_mean_row, right_mean_row, 
                                            left_squareDeviation, right_squareDeviation)
      }
    }
    likelihood = exp(log_likelihood)
    prob = likelihood/rowSums(likelihood)
  }
  # Get which genes have change points
  pi_0 = prob[,1];
  pi_0_order = sort(pi_0);
  pi_0_k=0;
  for(i in 1:N) {
    if(mean(pi_0_order[1:i]) <= FDR) {pi_0_k=pi_0_order[i];}
    else {break}
  }
  cp.index = which(pi_0 <= pi_0_k);
  
  # Get change points positions
  pi_star = apply(prob[,2:ncol(prob)], 1, max);
  pi_star_j = apply(prob[,2:ncol(prob)], 1, which.max);
  pi_star_1 = 1-pi_star
  pi_star_1_order = sort(pi_star_1);
  pi_star_1_k=0;
  for(i in 1:N) {
    if(mean(pi_star_1_order[1:i]) <= FDR) {pi_star_1_k=pi_star_1_order[i];}
    else {break}
  }
  pi_star_estimated = which(pi_star_1 <= pi_star_1_k);
  cp.position = data.frame(SeqID=pi_star_estimated, changePointTable[pi_star_j[pi_star_estimated],], stringsAsFactors=F)
  
  result = list(cp.index = cp.index, cp.position = cp.position, estimated_parameter = parameter_hat)
  return(result)
}
