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
LAMBDA <- 0.5 # Sparsity for binary matrix

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
# TODO (kevintee): better initialization based on distance
tau <- matrix(runif(Q*N), ncol=Q, nrow=N)
for(i in 1:N){
  normalize <- sum(tau[i,])
  for(q in 1:Q){
    tau[i,q] <- tau[i,q]/normalize
  }
}

iter <- 0
while(iter < 5){
  iter <- iter + 1
  # E-step: optimize each q(Z_i)
  # tau
  print tau[1:10,1:10]
  oldTau <- tau
  for(i in 1:N){
    for(q in 1:Q){
      total <- 0
      for(j in 1:N){
        for(l in 1:Q){
          if(i != j){
            total <- total + oldTau[j,l]*(digamma(xi[q,l])-digamma(eta[q,l]+xi[q,l])+X[i,j]*(digamma(eta[q,l])-digamma(xi[q,l])))
          }
        }
      }
      tau[i,q] <- exp(digamma(n[q])-digamma(sum(n)) + total)
    }
    # Normalize tau
    normalize <- sum(tau[i,])
    for(q in 1:Q){
      tau[i,q] <- tau[i,q]/normalize
    }
  }
  
  
  # M-step: optimize q(\alpha), q(\pi)
  # n
  n <- n + apply(tau, 1, sum)
  
  # eta
  for(q in 1:Q){
    total <- 0
    for(i in 1:N){
      if(i > 1){
        for(j in 1:(i-1)){
          total <- total + X[i,j]*tau[i,q]*tau[j,q]
        }
      }
    }
    eta[q,q] <- eta[q,q] + total
  }
  for(q in 1:Q){
    for(l in 1:Q){
      total <- 0
      for(i in 1:N){
        for(j in 1:N){
          if(i != j){
            total <- total + X[i,j]*tau[i,q]*tau[j,l]
          }
        }
      }
      eta[q,l] <- eta[q,l] + total
    }
  }
  
  # xi
  for(q in 1:Q){
    total = 0
    for(i in 1:N){
      if(i > 1){
        for(j in 1:(i-1)){
          total <- total + (1-X[i,j])*tau[i,q]*tau[j,q]
        }
      }
    }
    xi[q,q] = xi[q,q] + total
  }
  for(q in 1:Q){
    for(l in 1:Q){
      total <- 0
      for(i in 1:N){
        for(j in 1:N){
          if(i != j){
            total <- total + (1-X[i,j])*tau[i,q]*tau[j,l]
          }
        }
      }
      xi[q,l] <- xi[q,l] + total
    }
  }
}
