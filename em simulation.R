# cv for pca with em
library(MASS)
library(matlib)
library(graphics)
library(dplyr)
set.seed(1234)
p = 10
n = 200
#we generaate the data.  We intend to deal with p dimensional normal
impute = function(x, mu, sigma, indices){
  if(sum(indices) == 0){return(data)}
  inver = ginv(sigma[-indices, -indices])
  x[indices] <- mu[indices] + sigma[indices, -indices] %*% inver %*% (x[-indices] - mu[-indices])
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
  for(thing in 1:20){
    #Now, we calculate the svd and impute in every row
    #confirim if this is what was intended or if I should add random gaussian noise
    #Here, I had a matrix sigma. I took the svd with n components
    #next, I add some little noise to it.
    #next, I invert it to prevent the non-invertible mess
    
    #svd bit
    sigma <- first_n_svd(sigma, svd_dim)
    #impute in every row
    
    for(i in 1:n){
      tempx = hidden_data[i,]
      tempindices =  which(indices[i,] == 1)
      if(!sum(tempindices) == 0){
        
        hidden_data[i,] <- impute(tempx, mu, sigma, tempindices)
      }
    }

    #update the estimates
    nextmu = colSums(hidden_data)/n
    temp = matrix(0,nrow = p, ncol = p)
    for(i in 1:n){
      iter = hidden_data[i,] - mu
      temp = temp + iter%*%t(iter)
    }
    c_matrix = matrix(0, p,p)
    for(i in 1:p){
      for(j in 1:p){
        if(indices[i,j] == 1){
          c_matrix[i,j] = sigma[i,j]
        }
      }
    }
    nextsigma = temp + c_matrix
    mu = nextmu
    sigma = nextsigma/n
  }
  return(sum((hidden_data - realX)**2))
} #define the function to return the error after a few iterations
final_fun = function(X,iterations = 3, prop = 0.5){
    final_vec = rep(0, p)
    for(iteration in 1:iterations){
      realX<- X
      hidden_data = realX
     # indices = matrix(0, ncol = 2*p, nrow = 2*n)
      indices = matrix(0, ncol = p, nrow = n)
      for(i in 1:(n)){
        indicesini = as.integer(sample(p, as.integer(p*prop),replace = FALSE))  #nb of entries for a particular row
        indices[i,indicesini] = 1
      }
      mu<-c(1:(p))
      sigma = diag(seq(p,1))
      temp = seq(1,p)
      for(i in 1:(p)){
        temp[i] = error(hidden_data, mu, sigma, i,n,p,indices, realX)}

            final_vec[min(which(temp == min(temp)))] = final_vec[min(which(temp == min(temp)))] + 1
    }
    return(final_vec)
}
flatten(table)

table = matrix(nrow = 9, ncol = 10)
for(missdata in 1:9){
  for(runk in 1:10){
    
M<-matrix(rnorm(p*p,0,1), p, p)
Sigma_full_rank<-t(M)%*%M
Eig<-eigen(Sigma_full_rank)
V<-Eig$vectors
val<-Eig$values
val[-(1:runk)]<-0
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
  print(runk)
  table[missdata,runk] = which.max(final_fun(X, 100, missdata/10))
}
}






appa = """ [1]  1  2  3  4  5  6  7  8 10  8  1  2  3  4  5 10 10 10  5  5  1  2  3  4 10  4 10  2  2  1  1  2  3 10
[35]  2 10  1 10  1  1  1  2 10 10  1  1  1  1  1  1  1 10  1  1 10 10  1 10  1  1 10 10  1 10 10  1 10 10
[69] 10 10 10 10 10 10 10 10 10 10 10 10  1  1  1  1  1  1  1  1  1  1""" 







