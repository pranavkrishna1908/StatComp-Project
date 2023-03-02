datafunction = function(){
  CRTO <- read.csv("E:/EPFL/StatisticalComputation/StatComp-359822/CRTO.csv")
  DBVT <- read.csv("E:/EPFL/StatisticalComputation/StatComp-359822/DBVT.csv")
  EDAP <- read.csv("E:/EPFL/StatisticalComputation/StatComp-359822/EDAP.csv")
  AAP <- read.csv("E:/EPFL/StatisticalComputation/StatComp-359822/AAP.csv")
  AB <- read.csv("E:/EPFL/StatisticalComputation/StatComp-359822/AB.csv")
  AEO <- read.csv("E:/EPFL/StatisticalComputation/StatComp-359822/AEO.csv")
  AKR <- read.csv("E:/EPFL/StatisticalComputation/StatComp-359822/AKR.csv")
  ALSN <- read.csv("E:/EPFL/StatisticalComputation/StatComp-359822/ALSN.csv")
  AMC <- read.csv("E:/EPFL/StatisticalComputation/StatComp-359822/AMC.csv")
  AME<- read.csv("E:/EPFL/StatisticalComputation/StatComp-359822/AME.csv")
  crto = CRTO[,c(3,6)]
  dbvt = DBVT[,c(3,6)]
  edap = EDAP[,c(3,6)]
  aap = AAP[,c(3,6)]
  ab = AB[,c(3,6)]
  aeo = AEO[,c(3,6)]
  akr = AKR[,c(3,6)]
  alsn = ALSN[,c(3,6)]
  amc = AMC[,c(3,6)]
  ame = AME[,c(3,6)]
  crto = rev(log(crto[,2]/crto[,1]))[1:2050]
  dbvt = rev(log(dbvt[,2]/dbvt[,1]))[1:2050]
  edap = rev(log(edap[,2]/edap[,1]))[1:2050]
  aap = rev(log(aap[,2]/aap[,1]))[1:2050]
  ab = rev(log(ab[,2]/ab[,1]))[1:2050]
  aeo = rev(log(aeo[,2]/aeo[,1]))[1:2050]
  akr = rev(log(akr[,2]/akr[,1]))[1:2050]
  alsn = rev(log(alsn[,2]/alsn[,1]))[1:2050]
  amc = rev(log(amc[,2]/amc[,1]))[1:2050]
  ame = rev(log(ame[,2]/ame[,1]))[1:2050]
  data = cbind(crto, dbvt, edap, aap,ab, aeo, akr, alsn, amc, ame)
  data = as.data.frame(data)
  colnames(data) = c('CRTO', 'DBVT', 'EDAP', 'AAP', 'AB', 'AEO', 'AKR', 'ALSN', 'AMC', 'AME')
  return(data)
}
data = datafunction()
X = data
#em alogorithm

# cv for pca with em
library(MASS)
library(matlib)
library(graphics)
library(dplyr)
p = 10
n = 2050
impute = function(x, mu, sigma, indices){
  if(sum(indices) == 0){return(data)}
  inver = ginv(sigma[-indices, -indices])
  x[indices] <- mu[indices] + sigma[indices, -indices] %*% inver %*% t((x[-indices] - mu[-indices]))
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
  for(thing in 1:50){
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
      iter = as.vector(as.numeric(hidden_data[i,] - mu))
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
final_fun = function(X,iterations = 3){
  final_vec = rep(0, p)
  for(iteration in 1:iterations){
    print(iteration)
  realX<- X
    hidden_data = realX
    # indices = matrix(0, ncol = 2*p, nrow = 2*n)
    indices = matrix(0, ncol = p, nrow = n)
    for(i in 1:(n)){
      indicesini = as.integer(sample(p, as.integer(sample(p,1)/2)),replace = FALSE) #nb of entries for a particular row
      hidden_data[i,indicesini] = NA
      indices[i,indicesini] = 1
    }
    mu<-c(1:(p))
    M<-matrix(rnorm(p*p,0,1), p, p)
    realsigma<-t(M)%*%M
    sigma = realsigma
    temp = seq(1,p)
    for(i in 1:(p)){
      temp[i] = error(hidden_data, mu, sigma, i,n,p,indices, realX)}
    final_vec[min(which(temp == min(temp)))] = final_vec[min(which(temp == min(temp)))] + 1
    print(temp)
    }

  return(final_vec)
}
something<-final_fun(X,1)
print(something)
#weird behaviour of errors


data.pca <- prcomp(X,
                   center = TRUE,
                   scale. = TRUE)
thing = summary(data.pca)
relimp = thing$importance
props = as.vector(relimp[2,])
ggplot(data = as.data.frame(props))+
  geom_point(aes(x = as.factor(1:10), y = props))+
  labs( x = 'Index of Components', y = 'Proportion of Variability Explained')+
  theme(plot.title = element_text(hjust = 0.5), legend.position='none')




# loading library
library(ggfortify)
data.pca.plot <- autoplot(data.pca,
                          data = data)

data.pca.plot





"
> something<-final_fun(X,40)
[1] 1
 [1] 6.225861e+00 5.688068e+01 1.017916e+03 6.588370e+05 4.520820e+06 7.091415e+11 1.015682e+12 1.756427e+05 5.130132e+06 3.381779e+00
[1] 2
 [1] 1.205202e+01 8.748215e+03 4.257658e+04 3.551896e+04 1.538121e+06 1.001425e+14 2.054012e+05 1.456260e+09 4.455543e+09 3.218160e+00
[1] 3
 [1] 8.399067e+00 1.073108e+02 3.889361e+03 8.267208e+05 1.042076e+09 5.902387e+06 2.626850e+05 1.272813e+07 7.055099e+06 3.387501e+00
[1] 4
 [1] 5.899038e+00 3.882337e+05 7.259745e+03 9.647611e+04 3.006202e+06 2.204638e+07 8.887762e+17 1.512928e+08 2.982169e+05 3.434981e+00
[1] 5
 [1] 5.436877e+00 2.243384e+04 9.455703e+04 2.917248e+05 2.894514e+10 9.000529e+09 8.433409e+19 1.122245e+06 4.696488e+07 3.004981e+00
[1] 6
 [1] 1.414495e+01 8.167598e+01 2.129515e+02 5.252559e+05 8.634000e+05 2.502681e+10 1.886083e+06 8.601832e+06 1.770892e+07 3.300786e+00
[1] 7
 [1] 6.064166e+00 2.499026e+02 5.266636e+04 1.003300e+03 3.691175e+05 1.093530e+05 4.138194e+08 2.368611e+10 7.927596e+04 3.615910e+00
[1] 8
 [1] 7.590618e+00 2.079200e+02 5.529597e+04 5.248920e+03 1.097742e+06 3.109080e+15 4.335193e+09 3.380136e+06 1.733111e+06 3.235106e+00
[1] 9
 [1] 1.207476e+01 5.842069e+04 3.153926e+03 1.168354e+04 5.791447e+04 2.174670e+11 3.663464e+13 1.252279e+12 1.232828e+07 3.501719e+00
[1] 10
 [1] 9.735256e+00 1.425313e+03 2.835704e+03 1.059604e+04 6.092678e+06 1.338954e+09 4.895483e+14 1.040458e+06 1.299663e+08 3.244593e+00
[1] 11
 [1] 8.448631e+00 1.494352e+04 3.119102e+05 8.448675e+04 9.097645e+06 2.634249e+05 3.321999e+10 1.138610e+06 8.566185e+08 3.084986e+00
[1] 12
 [1] 7.224767e+01 3.018670e+02 2.351096e+03 2.015096e+05 4.731008e+05 1.187282e+08 3.605555e+12 2.044355e+09 3.204676e+09 3.299621e+00
[1] 13
 [1] 7.289560e+00 1.526975e+04 1.278769e+04 1.725375e+06 2.015796e+11 3.579443e+12 9.373465e+05 8.388058e+05 1.147048e+08 3.444840e+00
[1] 14
 [1] 5.003614e+02 3.815901e+04 4.035659e+03 2.185754e+09 6.036871e+04 7.536891e+05 2.128752e+09 1.234680e+05 6.096686e+05 3.190018e+00
[1] 15
 [1] 7.236515e+00 6.247114e+03 1.047497e+03 7.449036e+04 2.400080e+11 8.494740e+06 1.395439e+08 1.942949e+07 9.172979e+05 3.391911e+00
[1] 16
 [1] 7.069027e+00 3.672992e+06 4.130581e+04 2.728160e+05 7.803469e+07 1.953737e+07 4.433975e+05 2.403613e+07 6.206317e+05 3.335480e+00
[1] 17
 [1] 7.257610e+03 5.825991e+03 5.577534e+03 2.128745e+05 7.746279e+07 8.752685e+10 2.488046e+07 1.797349e+07 3.127404e+06 3.101840e+00
[1] 18
 [1] 7.158208e+00 1.576621e+04 1.422375e+07 2.498452e+04 1.072041e+06 4.890443e+10 2.092332e+17 1.268437e+09 3.581794e+05 3.455581e+00
[1] 19
 [1] 7.003316e+00 8.135312e+02 1.306608e+03 7.213655e+05 5.889773e+05 5.517992e+06 1.769587e+17 4.903218e+07 3.520194e+04 3.534151e+00
[1] 20
 [1] 6.577177e+00 4.960100e+01 1.456786e+04 1.101346e+05 2.176326e+09 3.109373e+07 5.375751e+06 1.996811e+06 1.019934e+08 3.589001e+00
[1] 21
 [1] 7.297255e+00 2.080490e+02 4.104210e+04 3.175295e+04 1.533182e+07 6.629741e+06 1.360181e+09 4.573906e+08 2.017130e+09 3.571246e+00
[1] 22
 [1] 1.043570e+01 1.558116e+01 6.171104e+04 4.082649e+06 2.029462e+09 6.421836e+07 3.624279e+12 3.305427e+07 9.488512e+08 3.048726e+00
[1] 23
 [1] 7.308882e+00 5.767271e+03 1.886105e+04 1.885176e+04 2.019226e+07 1.947274e+11 3.671006e+09 1.214128e+10 1.294800e+06 4.162272e+00
[1] 24
 [1] 5.773399e+00 1.731070e+02 5.993567e+03 2.201746e+04 7.134399e+06 3.714830e+06 3.515079e+06 2.616047e+08 9.170401e+05 3.325999e+00
[1] 25
 [1] 8.646853e+00 1.685102e+02 1.237877e+04 8.824715e+03 1.011977e+06 1.392163e+05 1.281930e+21 2.778962e+07 2.777574e+06 3.106436e+00
[1] 26
 [1] 7.139694e+00 1.106430e+03 8.549770e+03 1.848585e+05 4.335387e+08 2.596129e+07 1.717261e+08 9.104296e+08 1.833485e+09 3.250224e+00
[1] 27
 [1] 9.843820e+00 5.502066e+01 6.819597e+06 1.054984e+07 6.578016e+05 4.631082e+06 2.494134e+06 6.284344e+07 3.526755e+07 3.574016e+00
[1] 28
 [1] 1.078520e+01 3.606175e+02 1.819121e+03 1.356578e+05 6.054038e+07 9.210307e+22 1.772142e+04 6.982187e+05 1.994058e+06 3.255093e+00
[1] 29
 [1] 7.670115e+00 6.378010e+02 1.714185e+03 1.805096e+07 5.917166e+06 1.205498e+11 4.119366e+09 2.283652e+07 1.297866e+06 3.354391e+00
[1] 30
 [1] 1.645771e+01 2.401761e+02 4.209713e+04 2.711345e+05 3.507964e+08 2.005791e+10 6.746601e+05 3.100280e+07 7.739357e+07 3.397937e+00
[1] 31
 [1] 1.002213e+01 3.811246e+01 7.715230e+03 1.962293e+05 2.289344e+06 1.597524e+06 2.813940e+06 4.453066e+12 2.075296e+08 3.087207e+00
[1] 32
 [1] 8.103692e+00 3.217103e+02 5.527362e+05 1.251790e+05 3.622232e+07 1.224059e+09 5.830105e+26 5.320233e+06 1.569430e+06 3.308499e+00
[1] 33
 [1] 8.408934e+00 1.394817e+03 1.813781e+04 4.180679e+05 7.376338e+08 3.665746e+21 2.032243e+19 7.679082e+05 7.205019e+08 3.326184e+00
[1] 34
 [1] 6.702496e+00 3.197470e+01 1.630110e+04 2.163702e+06 5.499238e+06 1.345242e+09 6.871129e+06 1.308646e+16 1.529858e+06 3.207278e+00
[1] 35
 [1] 8.127006e+00 1.757815e+05 7.762708e+03 4.251217e+04 3.370935e+05 1.015363e+14 5.166663e+20 5.586294e+09 2.192912e+07 3.701137e+00
[1] 36
 [1] 1.113253e+01 1.381964e+04 3.687691e+05 7.250387e+04 5.136211e+05 4.683154e+12 5.142381e+10 5.668265e+08 7.739378e+06 3.273768e+00
[1] 37
 [1] 8.120513e+00 1.435600e+04 1.928762e+03 1.285443e+05 1.286135e+08 5.212631e+07 2.251987e+06 2.932279e+07 1.404403e+09 3.720176e+00
[1] 38
 [1] 8.907472e+00 1.152331e+03 1.063048e+03 4.769813e+05 6.794306e+06 3.145203e+06 9.792551e+14 1.825312e+10 1.357694e+10 3.156903e+00
[1] 39
 [1] 7.193167e+02 3.977272e+04 1.127525e+04 9.927387e+04 1.532441e+04 2.019462e+14 3.315865e+06 3.880629e+05 2.904028e+04 3.582748e+00
[1] 40
 [1] 5.823859e+00 2.521324e+02 7.127588e+02 3.407577e+05 5.180631e+05 9.188335e+05 7.204855e+07 8.773047e+11 3.008388e+08 3.370568e+00"


