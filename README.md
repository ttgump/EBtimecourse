# EBtimecourse

An empirical Bayes change point model for gene expression timecourse data.

Parameters for EBtimecourse function:

  exp.dat - matrix of data, rows stand for genes and columns stand for time points.

  timepoint -  number of time points

  replicate - number of replicates for each time point.

  FDR - expected false discovery rate (default 0.1).

  learning_rate - learning rate for the Adam optimizer (default 0.001).

  max_iter - max iterations (defualt 1e5).

  rel_tol - relative tolerance threshold for optimizing termination (defualt 1e-10).

  threads - number of threads.


The folder "Simulate_sample" stores code for simulation experiments. The folder "data" stores the real data we used in the paper.
