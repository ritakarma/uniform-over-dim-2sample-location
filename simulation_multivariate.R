source("functions.R")
 
## Load Packages
library("HDNRA")
library("mvtnorm")
library("highmean")
library("DescTools")
  
## specify the tests and local parameters below
## specify the global parameters below

N1 = 40
N2 = 50
P = 5
  
Mu = 6 #0:6
Mu = 2*Mu/10 #for p = 5
  
alpha = 0.05

# number of replications to calculate the proportion of rejections
Rep = 5000
  
## specify data generating models
model_indices = 1:9 #c(4,5,6,7,8,9)
  
Beta = 0.25 #c(0.1,0.2,0.3,0.4,0.5)  # tests corresponding to tapering estimators for different beta's
  
simulate_tests <- function() {
  ## Specify the parameters, models, tests etc below
  ## Output is a matrix with entries equal to the proportion of rejections
  ## and columns and rows corresponding to the models and tests respectively
    
  if (length(Beta) > 0) {
    tests = c("difplain", paste0("diftap_", Beta), "spaplain", paste0("spatap_", Beta))
  }
  else {
    tests = c("difplain", "spaplain")
  }
  
  ## Other tests
  #tests = c(tests, "ZGZC2020", "BS1996", "CLX2014", "CQ2010", "CLZ2014", "SD2008", "T^2")
  tests = c(tests, "ZGZC2020", "T^2", "CM1992")
  
  results = matrix(0, nrow = length(tests), ncol = length(model_indices), dimnames = list(tests, paste0("Ex ", as.character(model_indices))))
  
  for (id in 1:length(model_indices)) {
    vec_final = numeric(length(tests))  # vector to store total number of rejections
    print(model_indices[id])
    print(date())
    
    for (i  in 1:Rep) {
      res = model_sim(model_indices[id]) ## Generating dataX and dataY
      dataX = res$dataX
      dataY = res$dataY
      
      ## Performing our test with kernel h_diff, h_spatial and specified beta's
      ## 
      vec1 = perform_test(dataX = dataX, dataY = dataY, h = h_diff, nsim = 10000, vec_beta = Beta)
      vec2 = perform_test(dataX = dataX, dataY = dataY, h = h_spatial, nsim = 10000, vec_beta = Beta)
      vec = c(vec1, vec2)
      
      ## Performing other tests
      pval = HDNRA::ZGZC2020.TS.2cNRT(dataX, dataY)$p.value # Zhang et al 2020 Jasa (ZGZC2020)
      vec = c(vec, as.numeric(pval))   
      
      if(F){
      pval = highmean::apval_Bai1996(dataX, dataY)$pval # Bai and Saranadasa 1996 (BS1996)
      vec = c(vec, as.numeric(pval))
      
      pval = highmean::apval_Cai2014(dataX, dataY)$pval # Cai et al 2014 (CLX2014)
      vec = c(vec, as.numeric(pval))
      
      pval = highmean::apval_Chen2010(dataX, dataY)$pval # Chen and Qin 2010 (CQ2010)
      vec = c(vec, as.numeric(pval))
      
      pval = highmean::apval_Chen2014(dataX, dataY)$pval # Chen et al 2014 (CLZ2014)
      vec = c(vec, as.numeric(pval))
      
      pval = highmean::apval_Sri2008(dataX, dataY)$pval # Srivastava and Du 2008 (SD2008)
      vec = c(vec, as.numeric(pval))
      }
      
      pval = DescTools::HotellingsT2Test(x = dataX, y = dataY)$p.value # Hotelling's T^2 test
      vec = c(vec, as.numeric(pval))

      pval = choi_marden_2sample(dataX, dataY, h = h_spatial) # Choi and Marden (1992)
      vec = c(vec, as.numeric(pval))

      # reject null when p-value is less than alpha
      vec = as.numeric(vec <= alpha)
      
      vec_final = vec_final + vec
    }
    results[,id] = vec_final / Rep  # proportion of rejections
  }
  
  return(results)
}

model_sim <- function(index) {
  ## Input: index - model index
  ## Output: dataX, dataY
  
  n1 = N1
  n2 = N2
  p = P
  mu = Mu
  h_vec = (1:p)/sqrt(sum((1:p)^2))
  
  ## Models
  
  ##Model 1
  if (index == 1) {    
    ## X: N(0, Sigma), Y: N(0, Sigma), Sigma: equicorrelation matrix with non-diagonal
    ## elements 0.5
    
    Sigma = matrix(0.5, nrow = p, ncol = p)
    diag(Sigma) = rep(1, p)
    dataX = mvtnorm::rmvnorm(n1, mean = rep(0, p), sigma = Sigma)
    dataY = mvtnorm::rmvnorm(n2, mean = mu*h_vec, sigma = Sigma)
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ## Model 2
  if (index == 2) {
    ## X: t_4(0, Sigma), Y: t_4(0, Sigma), Sigma: equicorrelation matrix with non-diagonal
    ## elements 0.5
    
    Sigma = matrix(0.5, nrow = p, ncol = p)
    diag(Sigma) = rep(1, p)
    dataX = mvtnorm::rmvt(n1, sigma = Sigma, df = 4, delta = rep(0, p), type = "shifted")
    dataY = mvtnorm::rmvt(n2, sigma = Sigma, df = 4, delta = mu*h_vec, type = "shifted")
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ## Model 3
  if (index == 3) {
    ## X: t_1(0, Sigma), Y: t_1(0, Sigma), Sigma: equicorrelation matrix with non-diagonal
    ## elements 0.5
    
    Sigma = matrix(0.5, nrow = p, ncol = p)
    diag(Sigma) = rep(1, p)
    dataX = mvtnorm::rmvt(n1, sigma = Sigma, df = 1, delta = rep(0, p), type = "shifted")
    dataY = mvtnorm::rmvt(n2, sigma = Sigma, df = 1, delta = 10*mu*h_vec, type = "shifted")
    return(list(dataX = dataX, dataY = dataY))
  }

  ##Model 4
  if (index == 4) {
    ## X: N(0, Sigma), Y: N(0, Sigma), Sigma: identity
    
    Sigma = diag(1, nrow = p)
    dataX = mvtnorm::rmvnorm(n1, mean = rep(0, p), sigma = Sigma)
    dataY = mvtnorm::rmvnorm(n2, mean = mu*h_vec, sigma = Sigma)
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ## Model 5
  if (index == 5) {
    ## X: t_4(0, Sigma), Y: t_4(0, Sigma), Sigma: identity
    
    Sigma = diag(1, nrow = p)
    dataX = mvtnorm::rmvt(n1, sigma = Sigma, df = 4, delta = rep(0, p), type = "shifted")
    dataY = mvtnorm::rmvt(n2, sigma = Sigma, df = 4, delta = mu*h_vec, type = "shifted")
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ## Model 6
  if (index == 6) {
    ## X: t_1(0, Sigma), Y: t_1(0, Sigma), Sigma: identity
    
    Sigma = diag(1, nrow = p)
    dataX = mvtnorm::rmvt(n1, sigma = Sigma, df = 1, delta = rep(0, p), type = "shifted")
    dataY = mvtnorm::rmvt(n2, sigma = Sigma, df = 1, delta = 10*mu*h_vec, type = "shifted")
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ##Model 7
  if (index == 7) {
    ## X: N(0, Sigma), Y: N(0, Sigma), Sigma: (i,j)-th entry (0.7)^|i-j|
    
    DD = abs(outer(1:p, 1:p, "-"))
    Sigma = (0.75)^DD
    dataX = mvtnorm::rmvnorm(n1, mean = rep(0, p), sigma = Sigma)
    dataY = mvtnorm::rmvnorm(n2, mean = mu*h_vec, sigma = Sigma)
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ## Model 8
  if (index == 8) {
    ## X: t_4(0, Sigma), Y: t_4(0, Sigma), Sigma: (i,j)-th entry (0.7)^|i-j|
    
    DD = abs(outer(1:p, 1:p, "-"))
    Sigma = (0.75)^DD
    dataX = mvtnorm::rmvt(n1, sigma = Sigma, df = 4, delta = rep(0, p), type = "shifted")
    dataY = mvtnorm::rmvt(n2, sigma = Sigma, df = 4, delta = mu*h_vec, type = "shifted")
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ## Model 9
  if (index == 9) {
    ## X: t_1(0, Sigma), Y: t_1(0, Sigma), Sigma: (i,j)-th entry (0.7)^|i-j|
    
    DD = abs(outer(1:p, 1:p, "-"))
    Sigma = (0.75)^DD
    dataX = mvtnorm::rmvt(n1, sigma = Sigma, df = 1, delta = rep(0, p), type = "shifted")
    dataY = mvtnorm::rmvt(n2, sigma = Sigma, df = 1, delta = 10*mu*h_vec, type = "shifted")
    return(list(dataX = dataX, dataY = dataY))
  }
  
}

set.seed(121)
start = Sys.time()
res = simulate_tests()
end = Sys.time()
print((end - start))
print(res)
