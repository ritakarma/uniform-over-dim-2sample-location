# uniform-over-dim-2sample-location

This repository contains the relevant R code from the paper:

Ritabrata Karmakar, Joydeep Chowdhury, Subhajit Dutta, Marc G. Genton. (2025) Uniform-over-dimension location tests for
high-dimensional data. [arxiv preprint](https://arxiv.org/abs/2403.16328)

---

# Repository structure

- `functions.R`: Contains the main functions for the two sample test
  - `perform_test()`: Performs the main two sample location test and returns p-value
  - `choi_marden_2sample()`: Performs the kernel based test described by Choi and Marden (1997), used in the simulation study
- `simulations.R`: Script for the simulation study mentioned in the paper
- `colon_tissue_data_analysis.R`: Script for the real data analysis using colon tissue data from Alon et al. (1999)

---

# Example usage

To run the two sample test on the data matrices `dataX` and `dataY` using the spatial kernel and 10000 simulations for cut-off estimation, run the code
```r
result <- perform_test(dataX = dataX, dataY = dataY, h = h_spatial, estimators = c(1,0), nsim = 10000)
print(result)
```
To use tapering estimator based test with parameter values 0.1, 0.25 and 0.4, run the script
```r
result <- perform_test(dataX = dataX, dataY = dataY, h = h_spatial, estimators = c(0,1), nsim = 10000, vec_beta = c(0.1, 0.25, 0.4))
print(result)
```
In this case, an array of 3 p-values is returned. To use the difference kernel and both plain and tapering estimators run,
```r
result <- perform_test(dataX = dataX, dataY = dataY, h = h_diff, estimators = c(1,1), nsim = 10000, vec_beta = c(0.1, 0.25, 0.4))
print(result)
```
In this case, an array of 4 p-values is returned with the first element corresponding to the plain estimator based test and the remaining elements corresponding to the tapering estimator based tests.
