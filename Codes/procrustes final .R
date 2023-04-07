# cv for pca with em
library(MASS);library(matlib);library(graphics);library(mdatools);library(dplyr)
p = 10;n = 500
impute = function(x, mu, sigma, indices){
  if(sum(indices) == 0){return(data)}
  inver = ginv(sigma[-indices, -indices])
  test3 = try(x[indices] <- mu[indices] + sigma[indices, -indices] %*% inver %*% (x[-indices] - mu[-indices]))
  if (inherits(test3, "try-error")) browser()
  return(x)
} #define the function to impute in one row
first_n_svd = function(sigma, n){
  SVD<-svd(sigma)
  d<-SVD$d #vector containing the singular values of cov_est
  U<-SVD$u #matrix whose columns contain the left singular vectors of cov_est
  V<-SVD$v #matrix whose columns contain the right singular vectors of cov_est
  d[-(1:n)]<-0
  L<-U%*%diag(d)%*%t(V)
  sigma<-L  #+0.01*diag(nrow(L))  
  return(sigma)
} # svd of rank n, eckart young theorem etc
error = function(hidden_data, mu, sigma, svd_dim, n , p, indices,realX ){
  
  #Now, we calculate the svd and impute in every row
  #confirim if this is what was intended or if I should add random gaussian noise
  
  #svd bit
  test = try(sigma <- first_n_svd(sigma, svd_dim))
  if (inherits(test, "try-error")) browser()
  #impute in every row
  for(i in 1:n){
    tempx = hidden_data[i,]
    tempindices =  which(indices[i,] == 1)
    if(!sum(tempindices) == 0){
      
      hidden_data[i,] <- impute(tempx, mu, sigma, tempindices)
    }
  }
  return(sum((hidden_data - realX)**2))
} #define the function to return the error after a few iterations
final_fun = function(X,iterations = 3){
  final_vec = rep(0, p)
  for(iteration in 1:iterations){
    print(iteration)
    realX<- X
    #hide some data
    hidden_data = realX
    # indices = matrix(0, ncol = 2*p, nrow = 2*n)
    indices = matrix(0, ncol = p, nrow = n)
    for(i in 1:(n)){
      indicesini = as.integer(sample(p, as.integer(sample(p,1)/2)),replace = FALSE) #nb of entries for a particular row
      hidden_data[i,indicesini] = NA
      indices[i,indicesini] = 1
    }
    #generate the validation set
    validation = pcv(realX)
    #make the relevant estimates from the validation set
    mu<-colMeans(validation)
    temp = matrix(0,nrow = p, ncol = p)
    for(i in 1:n){
      iter = validation[i,] - mu
      temp = temp + iter%*%t(iter)
    }
    sigma = temp/n
    #fll in the error vector for each possible rank
    temp = seq(1,p)
    for(i in 1:(p)){
      temp[i] = error(hidden_data, mu, sigma, i,n,p,indices, realX)}
    #print(temp)
    final_vec[min(which(temp == min(temp)))] = final_vec[min(which(temp == min(temp)))] + 1
  }
  return(final_vec)
}
#we generate the covariance matrix of rank 3 which does not have zeroes
M<-matrix(rnorm(p*p,0,1), p, p)
Sigma_full_rank<-t(M)%*%M
Eig<-eigen(Sigma_full_rank)
V<-Eig$vectors
val<-Eig$values
val[-(1:3)]<-0
realsigma<-V%*%diag(val)%*%t(V)
sigma = realsigma
#generate the sample
mu = rnorm(p)
X<-matrix(0,n,p)
for (i in 1:n){
  x<-mvrnorm(n=1, mu, sigma)
  noise<-mvrnorm(n=10,0,1)
  X[i,]<-x + 0.01*noise
}
#doing the cross validation
something<-final_fun(X,70)
print(something)
# 0  0 12  0  0  0  0  0  0 38
# 0  0  3 41  0  0  0  0  0  6
# 0  0 50  0  0  0  0  0  0  0
# 0  0 48  0  0  0  0  0  0  2
# 0   0 500   0   0   0   0   0   0   0
# 0   0 500   0   0   0   0   0   0   0
#  0   0 263 420   0   0   0   0   0  17
# 0   0 700   0   0   0   0   0   0   0
# 0   0 700   0   0   0   0   0   0   0