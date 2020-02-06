###################################
# Algorithm
###################################

invariantPred <- function(Y,Xfull,indices,alpha){
  n_variables = dim(Xfull)[2] 
  
  l <- rep(list(0:1), n_variables)
  sets = expand.grid(l)
  n_sets = dim(sets)[1]
  
  n_envs = length(unique(indices))
  
  not_reject = integer(n_sets)
  
  for (s in 1:n_sets){
    set = sets[s,] == 1
    X = Xfull[,set,drop=FALSE]
    
    S = sum(sets[s,])
    pvalues = integer(n_envs)
    
    for (i in 1:n_envs) {
      sel = indices == i
      
      Ye = Y[sel]
      Xe = X[sel,,drop=FALSE]
      ne = length(Ye)
      
      Y_e = Y[!sel]
      X_e = X[!sel,,drop=FALSE]
      n_e = length(Y_e)
      
      # add intercept
      Xei = cbind(Xe,rep(1,dim(Xe)[1]))
      X_ei = cbind(X_e,rep(1,dim(X_e)[1]))
      
      XX = t(X_ei)%*%X_ei
      b_hat = solve(XX,t(X_ei)%*%Y_e)     
      
      d = Ye - Xei %*% b_hat 
      cov = diag(ne) + Xei %*% solve(XX,t(Xei))
      sigma_hat = var(Y_e - X_ei %*% b_hat)
      
      st = t(d)%*%solve(cov,d)/(sigma_hat * ne)
      
      pvalue = 1 - pf(st, ne, n_e - (S + 1)) 
      
      pvalues[i] = pvalue 
    }

    pvalue = min(pvalues)
    
    # reject/not reject H_{0,S}
    
    if(pvalue >= alpha/n_envs){
      not_reject[s]=1
    }
  } 
  
  # compute the S-intersection and the final parameters
  intersect_ = sapply(sets[not_reject == 1,],prod)
  S_hat = intersect_ == 1
  S = sum(intersect_)
  X = Xfull[,S_hat,drop=FALSE]
  n = dim(X)[1]
  
  # add intercept and solve the regression
  intercept = rep(1,dim(X)[1])
  Xi = cbind(X,intercept)
  XX = t(Xi)%*%Xi
  b_hat = solve(XX,t(Xi)%*%Y) 
  
  # confidence intervals
  sigma_hat = var(Y - Xi %*% b_hat)
  ci = qt(1-alpha/(2*(S+1)),n-(S+1)) * sqrt(diag(solve(XX)) * c(sigma_hat))
  
  # print results
  message(sprintf('%f CIs for the coefficients',1-alpha/(S+1)))
  for (i in 1:(S+1)){
    message(sprintf('%s : %f ± %f',attributes(ci)$names[i],b_hat[i],ci[i]))
  }
  
}

###################################
# Examples
###################################

# Graphical model
#  +--+    3     +--+   0.5    +--+    5     +--+
#  |X1| +------> |X2| +------> |Y | +------> |X4|
#  +--+          +--+          +--+          +--+
#  
#                               ^
#                               | 2
#                               |
#                               +
#
#                             +--+
#                             |X3|
#                             +--+
  

generate <- function(n,Y,X1,X2,X3,X4){
  if(missing(X1)){
    X1 = rnorm(n)  
  }
  if(missing(X2)){
    X2 = 3 * X1 + rnorm(n)
  }
  if(missing(X3)){
    X3 = rnorm(n)  
  }
  if(missing(Y)){
    Y = 0.5 * X2 + 2 * X3 + rnorm(n)  
  }
  if(missing(X4)){
    X4 = 5 * Y + rnorm(n)  
  }
  X = cbind(Y,X1,X2,X3,X4)
  return(X)
}

###################################
# Parameters
###################################

alpha = 0.05

######################################################################
# Perfect identification
# We provide one intervention per
# variable
######################################################################

n=500

# obs data
Xe1 = generate(n)
# intervention data
Xe2 = generate(n,X1=rep(1,n))
Xe3 = generate(n,X2=rep(2,n))
Xe4 = generate(n,X3=rep(3,n))
Xe5 = generate(n,X4=rep(4,n))
X = rbind(Xe1,Xe2,Xe3,Xe4,Xe5)
indices = c(rep(1,n),rep(2,n),rep(3,n),rep(4,n),rep(5,n))

Y = X[,1]
Xfull = X[,-1]

invariantPred(Y,Xfull,indices,alpha)

######################################################################
# Perfect identification
# Simultaneous interventions at causal variables
######################################################################

# obs data
Xe1 = generate(n)
# intervention data
Xe2 = generate(n,X1=rep(1,n))
Xe4 = generate(n,X3=rep(2,n),X2=rep(3,n))
Xe5 = generate(n,X4=rep(4,n))
X = rbind(Xe1,Xe2,Xe4,Xe5)
indices = c(rep(1,n),rep(2,n),rep(3,n),rep(4,n))

Y = X[,1]
Xfull = X[,-1]

invariantPred(Y,Xfull,indices,alpha)

######################################################################
# Perfect identification
# Simultaneous interventions
# We use soft-interventions to avoid getting a singular covariance matrix
# Since we only about Y|X remainig invariant, there is no problem
# having simultaneous interventions as long as Y is not intervened
######################################################################

# obs data
Xe1 = generate(n)
# intervention data
Xe2 = generate(n,X1=rnorm(n,2,4),X2=rnorm(n,2,5),X3=rnorm(n,2,6),X4=rnorm(n,2,7))
X = rbind(Xe1,Xe2)
indices = c(rep(1,n),rep(2,n))

Y = X[,1]
Xfull = X[,-1]

invariantPred(Y,Xfull,indices,alpha)

######################################################################
# Non-perfect identification
# No intervention for X4 and X2 is provided
# The null set is returned (i.e. only the intercept is used)
######################################################################

# obs data
Xe1 = generate(n)
# intervention data
Xe2 = generate(n,X1=rep(1,n))
Xe4 = generate(n,X3=rep(1,n))
X = rbind(Xe1,Xe2,Xe4)
indices = c(rep(1,n),rep(2,n),rep(3,n))

Y = X[,1]
Xfull = X[,-1]

invariantPred(Y,Xfull,indices,alpha)

######################################################################
# Non-perfect identification
# Soft-interventions cancel each other on average
# The empty set is returned since no causal predictor can be
# identified
######################################################################

# obs data
Xe1 = generate(n)
# intervention data
Xe3 = generate(n,X2=rnorm(n),X3=rnorm(n)*-1)
X = rbind(Xe1,Xe3)
indices = c(rep(1,n),rep(2,n))

Y = X[,1]
Xfull = X[,-1]

invariantPred(Y,Xfull,indices,alpha)
