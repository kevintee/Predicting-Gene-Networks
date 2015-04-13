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
EPSILON <- 1e-0 # Threshold for termination
Q <- 10 # Number of classes
N <- 1000 # Number of genes
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
tau <- tau/rowSums(tau) # Normalize
ones <- matrix(1, ncol=N, nrow=N)

while(sum(abs(oldTau-tau)) > EPSILON){
  # M-step: optimize q(\alpha), q(\pi)
  # n
  n <- n + apply(tau, 2, sum)

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
  alphaTerm <- digamma(matrix(rep(n, N), ncol=Q, nrow=N, byrow=T)) - digamma(matrix(sum(n), ncol=Q, nrow=N))
  piTerm <- ones%*%oldTau%*%(digamma(xi)-digamma(eta+xi)) + X%*%oldTau%*%(digamma(eta)-digamma(xi))
  tau <- exp(alphaTerm+piTerm)
  tau <- tau/rowSums(tau) # Normalize
  tau[is.na(tau)] <- 1/Q
}

# Print out readable results
tau <- round(tau, digits=2)
print(tau)
