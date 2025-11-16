h_diff <- function(dataX, dataY) {
  ##----------------------------------------------------------------------------
  ## Input: 
  ##   dataX: n1 * p matrix -rows are X observations  
  ##   dataY: n2 * p matrix -rows are Y observations
  ## Output: array of size n1 * n2 * p with (i, j, )-th entry X_i - Y_j
  ##----------------------------------------------------------------------------
  n1 = nrow(dataX)
  n2 = nrow(dataY)
  p = ncol(dataX)
  X = array(dataX, dim = c(n1, p, n2))  # n1 * p * n2 array where each slice ( , , k) is a copy of dataX
  X = aperm(X, c(1, 3, 2))              # swapped dimensions to make n1 * n2 * p array  
  
  Y = array(dataY, dim = c(n2, p, n1))  # n2 * p * n1 array where each slice ( , , k) is a copy of dataY
  Y = aperm(Y, c(3, 1, 2))              # swapped dimensions to make n1 * n2 * p array  
  return(X - Y)
}




h_spatial <- function(dataX, dataY, tol = 1e-8) {
  ##--------------------------------------------------------------------------------------
  ## Input: 
  ##   dataX: n1 * p matrix -rows are X observations  
  ##   dataY: n2 * p matrix -rows are Y observations
  ##   tol: tolerance to check if a number is zero
  ## Output: array of size n1 * n2 * p with (i, j, )-th entry (X_i - Y_j) / ||X_i - Y_j||
  ##         if ||X_i - Y_j|| > 0 and 0 otherwise
  ##--------------------------------------------------------------------------------------
  n1 = nrow(dataX)
  n2 = nrow(dataY)
  p = ncol(dataX)
  X = array(dataX, dim = c(n1, p, n2))       # n1 * p * n2 array where each slice ( , , k) is a copy of dataX
  X = aperm(X, c(1, 3, 2))                   # swapped dimensions to make n1 * n2 * p array  
  
  Y = array(dataY, dim = c(n2, p, n1))       # n2 * p * n1 array where each slice ( , , k) is a copy of dataY
  Y = aperm(Y, c(3, 1, 2))                   # swapped dimensions to make n1 * n2 * p array  
  D = X - Y                                  # n1 * n2 * p array with (i, j, )-th entry X_i - Y_j
  
  # swap columns and then sum over 1st dimension 
  weights = colSums(aperm((D * D), c(3, 1, 2)))  # n1 * n2 matrix with (i,j)-th entry ||D_(ij)||^2
  weights = sqrt(weights)
  id = (weights > tol)
  weights[id] = 1 / weights[id]                # (i,j)-th entry reciprocal of ||D_(ij)||
  weights[!id] = 0                             # (i,j)-th entry 0 if ||D_(ij)|| small
  
  return(D * array(weights, dim = c(n1, n2, p)))
}





compute_T_Sigma <- function(dataX, dataY, h) {
  ##----------------------------------------------------------------------------
  ## Input: 
  ##   dataX: n1 * p matrix -rows are X observations  
  ##   dataY: n2 * p matrix -rows are Y observations
  ##   h: kernel function - should take inputs dataX, dataY and give output an 
  ##      array of size n1 * n2 * p with (i, j, )-th entry h(X_i, Y_j)
  ## Output: 
  ##   A list containing -
  ##   T: the test statistic T_{n,p}
  ##   Sigma: estimator of Sigma_p  
  ##----------------------------------------------------------------------------
  
  n1 = nrow(dataX)
  n2 = nrow(dataY)
  p = ncol(dataX)
  
  #if (p != ncol(dataY)) {
  #  stop("Datasets have different dimensions")
  #}
  #if (n1 <= 1 || n2 <= 1 || p == 0) {
  #  stop("Must have rowsize of both datasets at least 2 and columnsize at least 1")
  #}
  
  n = n1 + n2
  H = h(dataX, dataY)     #dim: n1 * n2 * p
  
  if (!is.array(H) || !all(dim(H) == c(n1, n2, p))) {
    stop("Output of 'h' must be of dimension: n1 * n2 * p where n1 = nrow(dataX), n2 = nrow(dataY), p = ncol(dataX) = ncol(dataY)")
  }
  
  
  # sums over 1st dimension of H
  S2 = colSums(H)                    # S2: (n2 * p), j-th row is \sum_{i=1}^{n_2} h(X_i, Y_j)^t
  
  # swaps 1st and 2nd dimension and then sums over 1st dimension of H
  S1 = colSums(aperm(H, c(2,1,3)))   # S1: (n1 * p), i-th row is \sum_{j=1}^{n_2} h(X_i, Y_j)^t
  
  
  S = colSums(S1)             # vector (p * 1) - sum_i sum_j h(X_i, Y_j) 
  
  
  S = tcrossprod(S)           # p*p - (sum_i sum_j h(X_i, Y_j)) (sum_i sum_j h(X_i, Y_j))^t 
  
  S1 = crossprod(S1)          # p*p - sum_i (sum_j h(X_i, Y_j)) (sum_j h(X_i, Y_j))^t
  
  S2 = crossprod(S2)          # p*p - sum_j (sum_i h(X_i, Y_j)) (sum_i h(X_i, Y_j))^t
  
  
  sum_S_sq = sum(diag(S))     # scalar - ||sum_i sum_j h(X_i, Y_j)||^2
  
  sum_S1_sq = sum(diag(S1))   # scalar - sum_i ||sum_j h(X_i, Y_j)||^2
  
  sum_S2_sq = sum(diag(S2))   # scalar - sum_j ||sum_i h(X_i, Y_j)||^2
  
  sum_sq_h = sum(H * H)       # scalar - sum_i sum_j ||h(X_i, Y_j)||^2
  
  
  T_denom = n * n1 * n2
  
  T_val = ( sum_S_sq - sum_S1_sq - sum_S2_sq + sum_sq_h ) / T_denom
  
  Sigma = ( S1 + S2 ) / T_denom - S / (n1 * n1 * n2 * n2)
  
  return(list(T = T_val, Sigma = Sigma))
}






compute_Sigma_tap <- function(Sigma, n, beta) {
  ##--------------------------------------------------
  ## Input: 
  ##   Sigma: p * p matrix Sigma_p
  ##   n: (n1 + n2) - total number of observations
  ##   beta: tuning parameter beta > 0
  ## Output: 
  ##   tapered matrix Sigma_{p, tap} 
  ##--------------------------------------------------
  
  p = nrow(Sigma)
  
  pow = 1 / (2 * beta + 2)
  k = floor(n^pow)                                     
  if (k > p) {
    k = p
  }                                                    # k = min{p, floor(n^(1 / (2 * beta + 2)))}
  
  D = abs(outer(1:p, 1:p, "-"))                        # D_(ij) = |i-j|
  
  W = matrix(0, nrow = p, ncol = p)
  
  W[D <= (k/2)] = 1                                    # W_(ij) = 1 if |i-j| <= k/2
  
  middle_indices = (D > (k/2) & D < k)
  W[middle_indices] = 2 * (1 - D[middle_indices] / k)  # W_(ij) = 2 * (1 - D_(ij) / k) if |i-j| > k/2 and |i-j| < k
  
  return (Sigma * W)
}





perform_test <- function(dataX, dataY, h, nsim, estimators = c(1,1), vec_beta, batchsize = 10000) {
  ##----------------------------------------------------------------------------
  ## Input: 
  ##   dataX: n1 * p matrix -rows are X observations  
  ##   dataY: n2 * p matrix -rows are Y observations
  ##   h: kernel function - should take inputs dataX, dataY and give output an 
  ##      array of size n1 * n2 * p with (i, j, )-th entry h(X_i, Y_j)
  ##   nsim: number of simulations for cut-off calculation - default 1000
  ##   estimators: a vector of length 2 - if 1st entry is 1 plain estimator is 
  ##               used and if 2nd entry is 1 tapering estimator is used, both
  ##               estimators can be used simultaneously
  ##   vec_beta:  vector specifying beta parameter values for second estimator
  ##   batchsize: size of batches for simulation
  ## Output: 
  ##   a vector of p-values with names given by: 0 - plain estimator, 
  ##   i - tapering estimator corresponding to parameter vec_beta[i] for i > 0
  ##----------------------------------------------------------------------------
  
  if (all(estimators == 0)) {
    return(NULL)
  }
  
  res = compute_T_Sigma(dataX, dataY, h)
  T = res$T
  n = nrow(dataX) + nrow(dataY)
  
  # vector of indices corresponding to the estimators
  if (estimators[1] == 1 & estimators[2] == 1 & length(vec_beta) > 0) {
    indices = 0:length(vec_beta)
    vec_name = c("Plain", paste0("Tap-beta-", as.character(vec_beta))) 
  }
  else if (estimators[1] == 1){
    indices = 0
    vec_name = "Plain"
  }
  else {
    indices = 1:length(vec_beta)
    vec_name = paste0("Tap: beta ", as.character(vec_beta))
  }
  
  # list of p-values
  results = numeric(length(indices))
  names(results) = vec_name
  
  for (i in 1:length(indices)) {
    if (indices[i] == 0) {
      Sigma = res$Sigma  # plain estimator
    }
    else {
      Sigma = compute_Sigma_tap(res$Sigma, n, vec_beta[indices[i]])  #tapering estimator
    }
    
    # vector of eigenvalues
    lambda <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
    lambda <- pmax(lambda, 0)  # clip small negative eigenvalues if exists
    p = length(lambda)
    
    # vector of simulated random variables
    Q = numeric(nsim)
    
    # number of batches
    nbatches = ceiling(nsim / batchsize)
    
    # loop over batches
    for (b in 1:nbatches) {
      idx = ((b - 1) * batchsize + 1) : (min(b * batchsize, nsim))
      
      # generate standard normals
      Z = matrix(rnorm(length(idx) * p), nrow = length(idx), ncol = p)
      
      # compute quadratic forms for this batch
      Q[idx] = (Z^2 - 1) %*% lambda
    }
    
    # pvalue
    results[i] = mean(Q >= T) 
  }
  
  return(results)
}



choi_marden_2sample <- function(dataX, dataY, h = h_spatial) {
  ##----------------------------------------------------------------------------
  ## Input: 
  ##   dataX: n1 * p matrix -rows are X observations  
  ##   dataY: n2 * p matrix -rows are Y observations
  ##   h: kernel function (should take inputs dataX, dataY and give output an 
  ##      array of size n1 * n2 * p with (i, j, )-th entry h(X_i, Y_j))
  ##      --- h = h_diff: Hotelling's T^2 test with chi-square approximation
  ##      --- h = h_spatial (default): Spatial Rank test by Choi and Marden (1997) 
  ## Output: 
  ##   p-value of the test
  ##----------------------------------------------------------------------------
  
  n1 = nrow(dataX)
  n2 = nrow(dataY)
  n = n1 + n2
  p = ncol(dataX)
  
  #### Calculating Sigma_hat
  
  # h(dataX, dataX) is n1 * n1 * p array with (i,j, )-th entry being h(X_i, X_j)
  # S is n1 * p matrix with j-th row - \sum_{i=1}^{n_1} h(X_i, X_j)^t / n1 =  -R_N1(X^(j))^t as defined in Choi and Marden (1997)
  S = colSums(h(dataX, dataX)) / n1
  
  # S^t S = \sum_{j=1}^{n_1} R_N1(X^(j)) R_N1(X^(j))^t
  Sigma_hat = crossprod(S) 
  
  # Same calculation with the second group
  S = colSums(h(dataY, dataY)) / n2
  Sigma_hat = Sigma_hat + crossprod(S)
  Sigma_hat = Sigma_hat / (n-2)
  
  #### Calculating test statistic KW
  # n1*n*\bar{R}^(1) = \sum_{i=1}^n1 n * R(X^(i)) = \sum_{i=1}^n1 [\sum_i' h(X_i, X_i') + \sum_j h(X_i, Y_j)]
  #                 = \sum_{i=1}^n1 \sum_{i'=1}^n1 h(X_i, X_i') + \sum_{i=1}^n1 \sum_{j=1}^n1 h(X_i, Y_j)
  #                 = \sum_{i=1}^n1 \sum_{j=1}^n1 h(X_i, Y_j)   (the first term is zero since h(x,y) = -h(y,x))
  #                 = TS (say)
  # Similar result holds for group 2. Therefore, 
  # KW = n1 * \bar{R}^(1)^t Sigma_hat^{-1} \bar{R}^(1) + n2 * \bar{R}^(2)^t Sigma_hat^{-1} \bar{R}^(2)
  #    = (n1 / (n1^2 * n^2) + n2 / (n2^2 * n^2) ) TS^t Sigma_hat^{-1} TS
  #    = (1 / (n * n1 * n2)) TS^t Sigma_hat^{-1} TS
  
  TS = colSums(colSums(h(dataX, dataY)))
  KW = sum(TS * solve(Sigma_hat, TS)) / (n * n1 * n2)
  
  return(1-pchisq(KW, p))
}
