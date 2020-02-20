library(tensorflow)
library(foreach)

### auxiliary functions for normal normal gamm model
NNG_logL_cp <- function(prior_P, mu0, kappa0, alpha0, beta0, changePointTable, combNumber, Ti, left_seq_n, 
                         left_mean_row, right_mean_row, left_squareDeviation, right_squareDeviation, tf) {
  n1 = left_seq_n
  n2 = tf$constant(Ti, dtype=tf$float64) - left_seq_n
  
  kappa_n_1 = kappa0 + n1
  alpha_n_1 = alpha0 + n1/tf$constant(2, dtype=tf$float64)
  beta_n_1 = beta0 + left_squareDeviation*tf$constant(0.5, dtype=tf$float64) + 
    (n1/tf$constant(2, dtype=tf$float64)*kappa0*(left_mean_row-mu0)**2)/(kappa0+n1)
  kappa_n_2 = kappa0 + n2
  alpha_n_2 = alpha0 + n2/tf$constant(2, dtype=tf$float64)
  beta_n_2 = beta0 + right_squareDeviation*tf$constant(0.5, dtype=tf$float64) + 
    (n2/tf$constant(2, dtype=tf$float64)*kappa0*(right_mean_row-mu0)**2)/(kappa0+n2)
  
  log_pr = tf$log(prior_P) +
    tf$lgamma(alpha_n_1) - alpha_n_1*tf$log(beta_n_1) - tf$constant(0.5, dtype=tf$float64)*tf$log(kappa_n_1) +
    tf$lgamma(alpha_n_2) - alpha_n_2*tf$log(beta_n_2) - tf$constant(0.5, dtype=tf$float64)*tf$log(kappa_n_2) - 
    Ti/2*log(2*pi) - tf$constant(2, dtype=tf$float64)*tf$lgamma(alpha0) + 
    tf$constant(2, dtype=tf$float64)*alpha0*tf$log(beta0) + tf$log(kappa0)
  return(log_pr);
}

# object function for normal normal gamma model
NNG_obj <- function(prior_P, mu0, kappa0, alpha0, beta0, changePointTable, combNumber, Ti, left_seq_n, 
                    left_mean_row, right_mean_row, left_squareDeviation, right_squareDeviation, tf) {
  loss_cp <- NNG_logL_cp(prior_P, mu0, kappa0, alpha0, beta0, changePointTable, combNumber, Ti, 
                             left_seq_n, left_mean_row, right_mean_row, left_squareDeviation, right_squareDeviation, tf)
  loss <- tf$reduce_logsumexp(loss_cp, 1L)
  loss <- - tf$reduce_sum(loss)
  return(loss);
}


my_mean <- function(x) apply(x, 1, mean)
my_sumsquare <- function(x) apply(x, 1, function(x) sum(x^2))
my_squareDeviation <- function(x) apply(x, 1, function(x) sum((x-mean(x))^2))


EBtimecourse = function(exp.dat=NA, timepoint=NA, replicate=NA, FDR=0.1, 
                        learning_rate=0.001, max_iter=1e5, rel_tol=1e-10, threads=0L, verbose=F) {
  stopifnot(class(exp.dat) == "matrix")
  stopifnot(class(timepoint) == "numeric" & class(replicate) == "numeric" & class(FDR) == "numeric")
  stopifnot(ncol(exp.dat) == timepoint*replicate)
  
  tf <- tf$compat$v1
  tf$disable_v2_behavior()
  
  
  tf$reset_default_graph()
  
  # calculate global values
  N <- nrow(exp.dat)
  changePointTable <- data.frame(matrix(NA, nrow=(timepoint-1)+(timepoint-1)*(timepoint-2)/2, ncol=3), stringsAsFactors=F)
  colnames(changePointTable) <- c("n1", "n2", "n3")
  changePointTable[1:(timepoint-1),"n1"] <- 1:(timepoint-1)
  changePointTable[1:(timepoint-1),"n2"] <- (timepoint-1):1
  changePointTable[1:(timepoint-1),"n3"] <- 0
  combT <- as.data.frame(t(combn(timepoint-1, 2)))
  combT$V2 <- combT$V2 - combT$V1
  changePointTable[timepoint:nrow(changePointTable),c("n1","n2")] <- combT
  changePointTable[timepoint:nrow(changePointTable),"n3"] <- timepoint-rowSums(changePointTable[timepoint:nrow(changePointTable),c("n1", "n2")])
  combNumber <- nrow(changePointTable)
  Ti <- timepoint*replicate
  left_seq_index_table = foreach(cp = 1:(combNumber+1)) %do% {
    if(cp==combNumber+1) 
      return(1:Ti)
    else {
      rho1 <- changePointTable[cp,"n1"]*replicate
      rho2 <- (changePointTable[cp,"n1"]+changePointTable[cp,"n2"])*replicate
      
      if(rho2<Ti)
        return(c(1:rho1, (rho2+1):Ti))
      else
        return(c(1:rho1))
    }
  }
  left_seq_n <- sapply(left_seq_index_table, length)
  left_seq_n_ <- tf$constant(left_seq_n, dtype=tf$float64)
  left_seq_n_ <- tf$expand_dims(left_seq_n_, 0L)
  
  left_submatrix <- lapply(left_seq_index_table, function(i) exp.dat[, i, drop=F])
  right_submatrix <- lapply(left_seq_index_table, function(i) exp.dat[, -i, drop=F])
  
  left_mean_row <- sapply(left_submatrix, my_mean)
  right_mean_row <- sapply(right_submatrix, my_mean)
  right_mean_row[is.nan(right_mean_row)] <- 0
    
  left_squareDeviation <- sapply(left_submatrix, my_squareDeviation)
  right_squareDeviation <- sapply(right_submatrix, my_squareDeviation)
  
  # Data placeholders
  exp.dat_inp <- tf$placeholder(tf$float64, shape = shape(NULL, Ti), name = "exp_dat_inp")
  
  prior_P0 <- tf$Variable(tf$random_uniform(shape = shape(1), minval=1e-9, maxval=1-1e-9, dtype = tf$float64), 
                          dtype = tf$float64, name = "prior_P0", 
                          constraint = function(x) {
                            tf$clip_by_value(x, tf$constant(1e-9, dtype = tf$float64), 
                                             tf$constant(1-1e-9, dtype = tf$float64))
                          })
  W <- tf$Variable(tf$random_uniform(shape = shape(1), minval=1, maxval=100, dtype = tf$float64), 
                   dtype = tf$float64, name = "W", 
                   constraint = function(x) {
                     tf$clip_by_value(x, tf$constant(1, dtype = tf$float64), 
                                      tf$constant(1e5, dtype = tf$float64))
                   })
  mu0 <- tf$Variable(tf$random_normal(shape = shape(1), dtype = tf$float64), 
                     dtype = tf$float64, name = "mu0",
                     constraint = function(x) {
                       tf$clip_by_value(x, tf$constant(-100, dtype = tf$float64), 
                                        tf$constant(100, dtype = tf$float64))
                     })
  kappa0 <- tf$Variable(tf$random_uniform(shape = shape(1), minval=0.1, maxval=10, dtype = tf$float64), 
                        dtype = tf$float64, name = "kappa0",
                        constraint = function(x) {
                          tf$clip_by_value(x, tf$constant(1e-9, dtype = tf$float64), 
                                           tf$constant(100, dtype = tf$float64))
                        })
  alpha0 <- tf$Variable(tf$random_uniform(shape = shape(1), minval=0.1, maxval=10, dtype = tf$float64),
                        dtype = tf$float64, name = "alpha0",
                        constraint = function(x) {
                          tf$clip_by_value(x, tf$constant(1e-9, dtype = tf$float64), 
                                           tf$constant(100, dtype = tf$float64))
                        })
  beta0 <- tf$Variable(tf$random_uniform(shape = shape(1), minval=0.1, maxval=10, dtype = tf$float64),
                       dtype = tf$float64, name = "beta0",
                       constraint = function(x) {
                         tf$clip_by_value(x, tf$constant(1e-9, dtype = tf$float64), 
                                          tf$constant(100, dtype = tf$float64))
                       })
  p0 <- tf$constant(c(rep(1, combNumber), 0), dtype=tf$float64) + 
    tf$constant(c(as.integer(changePointTable$n1 != 1 & changePointTable$n3 != 1), 0), dtype=tf$float64) * (W-tf$constant(1, dtype=tf$float64))
  p0 <- p0 / tf$reduce_sum(p0) * (tf$constant(1, dtype=tf$float64) - prior_P0)
  t1 <- tf$constant(c(rep(0, combNumber), 1), dtype = tf$float64)
  prior_P <- p0 + t1 * prior_P0
  
  left_mean_row_inp <- tf$placeholder(tf$float64, shape = shape(NULL, ncol(left_mean_row)), name = "left_mean_row_inp")
  right_mean_row_inp <- tf$placeholder(tf$float64, shape = shape(NULL, ncol(right_mean_row)), name = "right_mean_row_inp")
  left_squareDeviation_inp <- tf$placeholder(tf$float64, shape = shape(NULL, ncol(left_squareDeviation)), name = "left_squareDeviation_inp")
  right_squareDeviation_inp <- tf$placeholder(tf$float64, shape = shape(NULL, ncol(right_squareDeviation)), name = "right_squareDeviation_inp")
  
  log_loss <- NNG_obj(prior_P, mu0, kappa0, alpha0, beta0, changePointTable, combNumber, Ti, left_seq_n_, 
                 left_mean_row_inp, right_mean_row_inp, left_squareDeviation_inp, right_squareDeviation_inp, tf)
  
  optimizer = tf$train$AdamOptimizer(learning_rate=learning_rate)
  train = optimizer$minimize(log_loss)
  
  
  # Start the graph and inference
  session_conf <- tf$ConfigProto(intra_op_parallelism_threads = threads,
                                 inter_op_parallelism_threads = threads)
  sess <- tf$Session(config = session_conf)
  init <- tf$global_variables_initializer()
  sess$run(init)
  
  fd_full <- dict(left_mean_row_inp = left_mean_row, 
                  right_mean_row_inp = right_mean_row, 
                  left_squareDeviation_inp = left_squareDeviation, 
                  right_squareDeviation_inp = right_squareDeviation)
  ll <- ll_old <- sess$run(log_loss, feed_dict = fd_full)
  
  ll_diff <- c()
  for(i in seq_len(max_iter)) {
    op <- sess$run(train, feed_dict = fd_full)

    ll <- sess$run(log_loss, feed_dict = fd_full)
    
    if(length(ll_diff) == 20) {
      ll_diff <- ll_diff[-1]
    }
    ll_diff <- c(ll_diff, (ll_old - ll))
    if(verbose & (i-1) %% 20 == 0) {
      message(paste(i, ll))
      message(paste("Max delta ll:", max(ll_diff)))
    }
    if(max(ll_diff) < rel_tol) {
      message(paste(i, ll))
      message(paste("Max delta ll:", max(ll_diff)))
      break
    }
    ll_old <- ll
  }
  
  variable_list <- list(prior_P0, W, mu0, kappa0, alpha0, beta0)
  variable_names <- c("P", "W", "mu0", "kappa0", "alpha0", "beta0")
  mle_params <- sess$run(variable_list, feed_dict = fd_full)
  names(mle_params) <- variable_names
  
  print("Converge params:")
  print(mle_params)
  
  log_likelihood_ <- NNG_logL_cp(prior_P, mu0, kappa0, alpha0, beta0, changePointTable, combNumber, Ti, left_seq_n_,
                                 left_mean_row_inp, right_mean_row_inp, left_squareDeviation_inp, right_squareDeviation_inp, tf)
  log_likelihood <- sess$run(log_likelihood_, feed_dict = fd_full)
  sess$close()
  
  likelihood = exp(log_likelihood)
  prob = likelihood/rowSums(likelihood)
  
  # Get which genes have change points
  pi_0 = prob[,ncol(prob)];
  pi_0_order = sort(pi_0);
  pi_0_k=0;
  for(i in 1:N) {
    if(mean(pi_0_order[1:i]) <= FDR) {pi_0_k=pi_0_order[i];}
    else {break}
  }
  cp.index = which(pi_0 <= pi_0_k);
  
  # Get change points positions
  pi_star = apply(prob[,1:(ncol(prob)-1)], 1, max);
  pi_star_j = apply(prob[,1:(ncol(prob)-1)], 1, which.max);
  pi_star_1 = 1-pi_star
  pi_star_1_order = sort(pi_star_1);
  pi_star_1_k=0;
  for(i in 1:N) {
    if(mean(pi_star_1_order[1:i]) <= FDR) {pi_star_1_k=pi_star_1_order[i];}
    else {break}
  }
  pi_star_estimated = which(pi_star_1 <= pi_star_1_k);
  cp.position = data.frame(SeqID=pi_star_estimated, changePointTable[pi_star_j[pi_star_estimated],], stringsAsFactors=F)
  
  result = list(cp.index = cp.index, cp.position = cp.position, mle_parameter = mle_params, ll=likelihood)
  return(result)
}