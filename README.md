# EBtimecourse

An empirical Bayes change point model for gene expression timecourse data.

Parameters for EBtimecourse function:

exp.dat - matrix of data, rows stand for genes and columns stand for time points.

timepoints -  vector of time points. For example, if we have 5 time points and each time point has 2 replicates, then timepoints shold be c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)

FDR - expected false discovery rate (default 0.1).

learning_rate - learning rate for the Adam optimizer (default 0.001).

max_iter - max iterations (defualt 1e5).

rel_tol - relative tolerance threshold for optimizing termination (defualt 1e-10).

threads - number of threads.

