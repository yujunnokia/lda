# TODO: convert the click data into the format into LDA model
# 
# Author: Jun Yu
# Version: June, 2012
###############################################################################

rm(list=ls())

setwd("/Users/junyu/Documents/code/eBay/code/lda/")

dyn.load("src/gibbsLDA.so")

library("ggplot2")
library("reshape")

source("R/lda.collapsed.gibbs.sampler.R")
source("R/mblda.collapsed.gibbs.sampler.R")
source("R/top.topic.words.R")
source("R/top.topic.documents.R")

#######################
## experiment settings
#######################
ldaVersion <- "lda"  # lda or mblda
query <- "ipod"  # query term
numClicks <- "10000"  # number of days' data before endScandate
K <- 5  # number of topics
nIterations <- 5000 # number of gibbs sampling iterations
alpha <- 0.1  # alpha hyperparameter
eta <- 0.1  # eta hyperparameter

######################
## load in click data
######################
data(cora.documents)
data(cora.vocab)
load(paste("../../data/eBay/LDA/",query,"_",numClicks,"_train.RData",sep=""))

theme_set(theme_bw())  
set.seed(8675309)


if (ldaVersion == "lda") {
	# run LDA
	result <- lda.collapsed.gibbs.sampler(lda.documents,
			K,  ## Num clusters
			vocab,
			nIterations,  ## Num iterations
			alpha,
			eta,
			compute.log.likelihood=TRUE) 

	# Get the top words in the cluster
	top.words <- lda.top.topic.words(result$topics, 10, by.score=TRUE)
	top.docs  <- lda.top.topic.documents(result$document_sums, num.documents = 7, alpha = 0.1) 

	# get word-weight matrix
	word.proportions <- result$topics/(rowSums(result$topics) + 1e-05)
} else if (ldaVersion == "mblda") {
	# run MBLDA
	result <- mblda.collapsed.gibbs.sampler(mblda.documents,
			K,  ## Num clusters
			vocab,
			nIterations,  ## Num iterations
			alpha,
			eta,
			compute.log.likelihood=FALSE) 	
	
	# Get the top words in the cluster
	top.words <- mblda.top.topic.words(result$topics_present, result$topics_absent, 10, by.score=TRUE)
	
} else {
	stop("lda version is invalid...")
}


## Display document proportions
#N <- 10
#topic.proportions <- t(result$document_sums) / colSums(result$document_sums)
#sample.docs <- 1:N #sample(1:length(lda.documents),N)
#sample.docs.terms <- do.call(rbind, documents.origtitle[sample.docs])
#topic.proportions <- topic.proportions[sample.docs,]
#topic.proportions[is.na(topic.proportions)] <-  1 / K
#
#colnames(topic.proportions) <- apply(top.words, 2, paste, collapse=" ")
#topic.proportions.df <- melt(cbind(data.frame(topic.proportions),
#				document=sample.docs.terms),	#document=sample.docs.terms), #document=factor(1:N)),
#		variable_name="topic",
#		id.vars = "document")  
#
#
#qplot(topic, value, fill=document, ylab="proportion", main = query,
#				data=topic.proportions.df, geom="bar") +
#		opts(axis.text.x = theme_text(angle=90, hjust=1)) +  
#		coord_flip() +
#		facet_wrap(~ document, ncol=5)

