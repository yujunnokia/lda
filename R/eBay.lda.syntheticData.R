# TODO: generate synthetic data as a sanity check of lda models
#       you can compare the true model with the estimated model by eye-balling the top words of each topic 
#       note that topic orderings might be different between true and estimated models
# 
# Author: Jun Yu
# Version: June, 2012
###############################################################################

rm(list=ls())

setwd("/Users/junyu/Documents/eBay/code/lda")

library("MCMCpack")
library("hash")

dyn.load("src/gibbsLDA.so")

source("R/lda.collapsed.gibbs.sampler.R")
source("R/mblda.collapsed.gibbs.sampler.R")

set.seed(8675309)

######################
# experiment settings
######################
alpha <- 0.1
eta <- 0.1
K <- 5
nDocuments <- 500
documentLengthMean <- 15
documentLengthVar <- 2
nTopWords <- 10
vocab <- as.character(0:99)
documents <- list()

###########################
# generate synthetic data
###########################
beta <- matrix(0, nrow=K, ncol=length(vocab))
for (k in 1:K) {
	beta[k,] <- rdirichlet(1, rep(eta,length(vocab)))
}
colnames(beta) <- vocab

length <- round(rnorm(nDocuments, documentLengthMean, documentLengthVar))
theta <- matrix(0, nrow=K, ncol=nDocuments)
for (i in 1:nDocuments) {
	theta[,i] <- rdirichlet(1, rep(alpha,K))
	
	if (length[i] <= 0) { length[i] <- documentLengthMean }
	
	document <- matrix(0, nrow=2, ncol=length(vocab))
	document[1,] <- 0:99
	Zs <- rmultinom(1, size = length[i], prob=theta[,i])
	for (k in 1:K) {
		words <- rmultinom(1, size = Zs[k,1], beta[k,])
		document[2,] <- document[2,] + words
	}
	
	document <- document[,- which(document[2,] == 0)]
	
	document <- as.integer(document)
	document <- matrix(document,nrow=2)
	documents[[i]] <- document
}


#######################
# train lda model
#######################
nIterations <- 1000
result <- lda.collapsed.gibbs.sampler(documents,
		K,  ## Num clusters
		vocab,
		nIterations,  ## Num iterations
		alpha,
		eta,
		compute.log.likelihood=TRUE) 

# compute word proportion matrix
word.proportions <- result$topics/(rowSums(result$topics) + 1e-05)


#############
# evaluate 
#############

# get top words from true model
model.top.words <- matrix(0, nrow=nTopWords, ncol=2*K)
for (k in 1:K) {
	word.distribution <- beta[k,]
	ordering <- order(word.distribution, decreasing = TRUE)[1:nTopWords]
	model.top.words[,2*k-1] <- colnames(beta)[ordering]
	model.top.words[,2*k] <- word.distribution[ordering]
}
cat("True model...\n")
print(model.top.words)

# get top word from estimated model
est.top.words <- matrix(0, nrow=nTopWords, ncol=2*K)
for (k in 1:K) {
	word.distribution <- word.proportions[k,]
	ordering <- order(word.distribution, decreasing = TRUE)[1:nTopWords]
	est.top.words[,2*k-1] <- colnames(word.proportions)[ordering]
	est.top.words[,2*k] <- word.distribution[ordering]
}
cat("Estimated model...\n")
print(est.top.words)

