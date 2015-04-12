# Variational Bayes EM for Stochastic Block Model
# Based on Latouche et al., 2010

# Set working directory
setwd("/Users/kevintee/Downloads/Predicting-Gene-Networks/src/")

# Read in files
# X is in the form: Genes vs Tumor Sample
fileName <- "data/KIRC.txt"
rawData <- read.table(fileName, header=T, sep="\t")
rawData <- rawData[,-1] # Remove first column

# Define constants
EPSILON <- 1e-6 # Threshold for termination
Q <- 10 # Number of classes
N <- 100 # Number of genes
LAMBDA <- 0.3 # Sparsity for binary matrix

rawData <- rawData[1:N,] # Just take 100 for now for speed

# Generate binary matrix via correlation cutoff
# TODO(kevintee): better way of generating binary matrix
X <- cor(t(rawData), t(rawData))
X[abs(X) < LAMBDA] <- 0 # Remove uncorrelated values
X[X != 0] <- 1 # Set all nonzero values to 1 for binary matrix

# Initialization
n <- vector(mode="integer", length=Q) + 0.5
eta <- matrix(0.5, ncol=Q, nrow=Q)
xi <- matrix(0.5, ncol=Q, nrow=Q)
oldTau <- matrix(0, ncol=Q, nrow=N)
# TODO (kevintee): better initialization based on distance
tau <- matrix(runif(Q*N), ncol=Q, nrow=N)
for(i in 1:N){
  normalize <- sum(tau[i,])
  for(q in 1:Q){
    tau[i,q] <- tau[i,q]/normalize
  }
}

iter <- 0
while(iter < 100){
  iter <- iter + 1

  # M-step: optimize q(\alpha), q(\pi)
  # n
  n <- n + apply(tau, 1, sum)

  # eta
  nonDiagX <- X
  diag(nonDiagX) <- 0
  total <- t(tau)%*%nonDiagX%*%tau
  diag(total) <- diag(total)/2
  eta <- eta + total

  # xi
  nonDiagX <- 1-X
  diag(nonDiagX) <- 0
  total <- t(tau)%*%nonDiagX%*%tau
  diag(total) = diag(total)/2
  xi <- xi + total

  # E-step: optimize each q(Z_i)
  # tau
  oldTau <- tau
  ones <- matrix(1, ncol=N, nrow=N)
  alphaTerm = digamma(replicate(Q, n)) - digamma(matrix(sum(n), ncol=Q, nrow=N))
  piTerm = ones%*%oldTau%*%(digamma(xi)-digamma(eta+xi)) + X%*%oldTau%*%(digamma(eta)-digamma(xi))
  tau <- exp(alphaTerm+piTerm)
  for(i in 1:N){
    # Normalize tau
    normalize <- sum(tau[i,])
    for(q in 1:Q){
      tau[i,q] <- tau[i,q]/normalize
    }
  } 
  print(tau[1:10,1:10])
}
