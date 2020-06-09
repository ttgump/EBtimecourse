# EBtimecourse

An empirical Bayes change point model for gene expression timecourse data.

Required R packages: tensorflow, foreach

Use soure("EBtimecourse.R") to load the function in the R enviroment.

The function only accepts one conditional timecourse data. If you have two conditional data, say cases and controls, you can transform it to one conditional data: case - control for paired data, averaged cases - averaged controls for unpaired data.

Parameters for EBtimecourse function:
```
   exp.dat - matrix of data, rows stand for genes and columns stand for time points.

   timepoint -  number of time points

   replicate - number of replicates for each time point.

   FDR - expected false discovery rate (default 0.1).

   learning_rate - learning rate for the Adam optimizer (default 0.001).

   max_iter - max iterations (defualt 1e5).

   rel_tol - relative tolerance threshold for optimizing termination (defualt 1e-10).

   threads - number of threads.
```

If we have a data of 8 time points, and each time point has 3 replicates, then timepoint=8 and replicate=c(3,3,3,3,3,3,3,3).

The folder "Simulate_sample" stores code for simulation experiments. The folder "data" stores the real data we used in the paper.
