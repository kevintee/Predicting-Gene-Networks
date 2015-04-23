# Predicting Cancer Regulatory Networks
# Variational Bayes EM for Stochastic Block Model (Latouche et al., 2010)

# Set working directory
setwd("/Users/kevintee/Downloads/Predicting-Gene-Networks/src/")

# Read in files
# X is in the form: Genes vs Tumor Sample
fileName <- "data/KIRC.txt"
rawData <- read.table(fileName, header=TRUE, sep="\t")
genes <- rawData[1] # Gene Names
genes <- matrix(genes[,c("Name")]) # Convert from dataframe to matrix
rawData <- rawData[,-1] # Remove first column

# Define constants
EPSILON <- 1e-0 # Threshold for termination
Q <- 100 # Number of classes
N <- 1000 # Number of genes
LAMBDA <- 0.5 # Sparsity for binary matrix

rawData <- rawData[1:N,] # Take N for speed

# Sample covariance
X <- cov(t(rawData), t(rawData))

# Remove uncorrelated values and set all nonzero values to 1
X[abs(X) < LAMBDA] <- 0
weights <- X
X[X != 0] <- 1

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
  alphaTerm <- digamma(matrix(rep(n, N), ncol=Q, nrow=N, byrow=TRUE)) - digamma(matrix(sum(n), ncol=Q, nrow=N))
  piTerm <- ones%*%oldTau%*%(digamma(xi)-digamma(eta+xi)) + X%*%oldTau%*%(digamma(eta)-digamma(xi))
  tau <- exp(alphaTerm+piTerm)
  tau <- tau/rowSums(tau) # Normalize
  tau[is.na(tau)] <- 1/Q
}

# Write weights of nonzero gene-gene edges
weightFile <- "../results/sbm/KIRC_weights.txt"
weights <- round(weights, digits=3)
write.table(weights, file=weightFile, quote=FALSE, col.names=FALSE, row.names=FALSE)

# Get cluster assignments
clusters <- matrix(0, ncol=2, nrow=N)
for(i in 1:N){
  clusters[i,] <- c(genes[i], which.max(tau[i,]))
}

# Write clusters
clusterFile <- "../results/sbm/KIRC_cluster.txt"
write.table(clusters, file=clusterFile, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
