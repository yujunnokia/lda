# TODO: generate synthetic data from mblda model as a sanity check
# 
# Author: Jun Yu
# Version: June, 2012
###############################################################################

rm(list=ls())

setwd("/Users/junyu/Documents/workspace/lda")

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
vocab <- as.character(0:99)
documents <- list()

###########################
# generate synthetic data
###########################
#beta <- matrix(0, nrow=K, ncol=length(vocab))
#for (k in 1:K) {
#	for (w in 1:length(vocab)) {
#		beta[k,w] <- rdirichlet(1, rep(eta,2))[1]
#	}
#}
#colnames(beta) <- vocab
beta <- matrix(0, nrow=K, ncol=length(vocab))
beta[1,1:10] <- sapply(rnorm(10, 0.9, 0.2), function(x) min(x,0.95))
beta[2,15:20] <- sapply(rnorm(6, 0.9, 0.2), function(x) min(x,0.95))
beta[3,35:40] <- sapply(rnorm(6, 0.9, 0.2), function(x) min(x,0.95))
beta[4,95:100] <- sapply(rnorm(6, 0.9, 0.2), function(x) min(x,0.95))
beta[5,55:60] <- sapply(rnorm(6, 0.9, 0.2), function(x) min(x,0.95))
colnames(beta) <- vocab

length <- round(rnorm(nDocuments, documentLengthMean, documentLengthVar))
theta <- matrix(0, nrow=K, ncol=nDocuments)
for (i in 1:nDocuments) {
	theta[,i] <- rdirichlet(1, rep(alpha,K))
	
	if (length[i] <= 0) { length[i] <- documentLengthMean }
	
	document <- matrix(0, nrow=2, ncol=length(vocab))
	document[1,] <- 0:99
	for (w in 1:length(vocab)) {
		z <- which(rmultinom(1, size = 1, prob=theta[,i]) == 1)
		
		word <- which(rmultinom(1, size = 1, c(beta[z,w], 1-beta[z,w])) == 1)
		if (word == 1) {
			document[2,w] <- 1	
		}
	}

	document <- as.integer(document)
	document <- matrix(document,nrow=2)
	documents[[i]] <- document
}


#######################
# train lda model
#######################
nIterations <- 1000
result <- mblda.collapsed.gibbs.sampler(documents,
		K,  ## Num clusters
		vocab,
		nIterations,  ## Num iterations
		alpha,
		eta,
		compute.log.likelihood=TRUE) 

# get word-weight matrix
normalized.topics <- result$topics_present/(result$topics_present + result$topics_absent + 1e-05)

#############
# evaluate 
#############
nTopWords <- 10
model.top.words <- matrix(0, nrow=nTopWords, ncol=2*K)
for (k in 1:K) {
	word.distribution <- beta[k,]
	ordering <- order(word.distribution, decreasing = TRUE)[1:nTopWords]
	model.top.words[,2*k-1] <- colnames(beta)[ordering]
	model.top.words[,2*k] <- word.distribution[ordering]
}

# get top word for each topic with probability
est.top.words <- matrix(0, nrow=nTopWords, ncol=2*K)
for (k in 1:K) {
	word.distribution <- normalized.topics[k,]
	ordering <- order(word.distribution, decreasing = TRUE)[1:nTopWords]
	est.top.words[,2*k-1] <- colnames(normalized.topics)[ordering]
	est.top.words[,2*k] <- word.distribution[ordering]
}

print(model.top.words)
print(est.top.words)