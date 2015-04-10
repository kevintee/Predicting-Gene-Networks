# Variational Bayes EM for Stochastic Block Model
# Based on Latouche et al., 2010

# Set working directory
setwd("/Users/kevintee/Downloads/Predicting-Gene-Networks/src/")

# Read in files
fileName <- "data/BRCA.txt"
allData <- read.table(fileName, header=T, sep="\t")
allData <- allData[,-1]
