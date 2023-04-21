# EBtimecourse

An empirical Bayes change point model for gene expression timecourse data. Time-course experiments are commonly conducted to capture temporal changes. It is generally of interest to detect if any changes happen over time, which we deﬁne as a detection problem. If there is a change, it is informative to know when the change is, which we deﬁne as an identiﬁcation problem. It is often desired to control Type I error rate at a nominal level while applying a testing procedure to detect or identify these changes. The EBtimecourse model provides an unified multiple-testing framework built upon an empirical Bayes change-point model to solve these two problems. The detail is described in our paper "An Empirical Bayes Change Point Model for Transcriptome Time Course Data" published in *Annals of Applied Statistics*. https://doi.org/10.1214/20-AOAS1403

![alt text](https://github.com/ttgump/EBtimecourse/blob/master/Diagram.png?raw=True)<br/>
In the figure, Xij means the expression level of gene i at time point j.

Required R packages: tensorflow, foreach

Use ```soure("EBtimecourse.R")``` to load the function in the R enviroment.

The function only accepts one conditional timecourse data. The input data should be log scaled microarray data or normalized and log scaled RNA-seq data. If you have two conditional data, say cases and controls, you can transform it to one conditional data: case - control for paired data, averaged cases - averaged controls for unpaired data.

Parameters for EBtimecourse function:
```
   exp.dat - matrix of data, rows stand for genes and columns stand for time points.

   timepoint -  number of time points

   replicate - a vector of number of replicates for each time point.

   FDR - expected false discovery rate (default 0.1).

   learning_rate - learning rate for the Adam optimizer (default 0.001).

   max_iter - max iterations (defualt 1e5).

   rel_tol - relative tolerance threshold for optimizing termination (defualt 1e-10).

   threads - number of threads.
```

If you have a data of 8 time points, and each time point has 3 replicates, then timepoint=8 and replicate=c(3,3,3,3,3,3,3,3).

The folder "Simulate_sample" stores code for simulation experiments. The folder "data" stores the real data we used in the paper.

Here is version information of my R enviroment:

```
R version 3.6.0 (2019-04-26)
Platform: x86_64-apple-darwin18.6.0 (64-bit)
Running under: macOS  10.15.4

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /usr/local/Cellar/openblas/0.3.6_1/lib/libopenblasp-r0.3.6.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] foreach_1.4.7    tensorflow_2.0.0

loaded via a namespace (and not attached):
 [1] compiler_3.6.0   magrittr_1.5     Matrix_1.2-17    tools_3.6.0      whisker_0.4      base64enc_0.1-3  Rcpp_1.0.3      
 [8] reticulate_1.13  codetools_0.2-16 grid_3.6.0       iterators_1.0.12 jsonlite_1.6     tfruns_1.4       lattice_0.20-38 
```
