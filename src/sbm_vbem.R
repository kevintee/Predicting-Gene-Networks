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
N <- 8499 # Number of genes
LAMBDA <- 0.5 # Sparsity for binary matrix

# Generate binary matrix via correlation cutoff
# TODO(kevintee): better way of generating binary matrix
X = cor(t(rawData), t(rawData))
X[abs(X) < LAMBDA] = 0 # Remove uncorrelated values
X[X != 0] = 1 # Set all nonzero values to 1 for binary matrix
