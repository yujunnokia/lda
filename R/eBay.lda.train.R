# TODO: train the LDA model
# 
# Author: Jun Yu
# Version: June, 2012
###############################################################################

rm(list=ls())

#setwd("/Users/junyu/Documents/workspace/lda")
#setwd("/nfs/stak/students/y/yuju/exp/eBay/code")
setwd("/Users/junyu/Documents/eBay/code/lda")

dyn.load("src/gibbsLDA.so")

source("R/lda.collapsed.gibbs.sampler.R")
source("R/mblda.collapsed.gibbs.sampler.R")

set.seed(8675309)

#######################
## experiment settings
#######################
ldaVersion <- "lda"  # lda or mblda
queries <- c("hello kitty", "fossil", "fossils", "basketball", "keyboard", "iphone", "ipod", "coach")
Ks <- c(4)  # number of topics
nIterations <- 2000 # number of gibbs sampling iterations
alpha <- eta <- 0.1  # hyperparameter
nTopWords <- 10  # number of top words in each topic

args <- commandArgs(trailingOnly = TRUE)
ldaVersion <- as.character(args[2])
if (ldaVersion != "lda" && ldaVersion != "mblda") {
	stop("lda version is invalid.\n
		please run 'Rscript eBay.lda.train.R -args lda' or \n
		'Rscript eBay.lda.train.R -args mblda'")
} 

## load queries
#queryFile <- "../../data/testQueries.tsv"
#queryData <- read.table(queryFile, head=TRUE, sep="\t", stringsAsFactors=FALSE)
#queries <- queryData$query

for (query in queries) {
	print(query)
	
	# load in train data
	load(paste("../../data/train/",query,"_train.RData",sep=""))
	
	for (K in Ks) {
		if (ldaVersion == "lda") {
			# run LDA
			result <- lda.collapsed.gibbs.sampler(lda.documents,
					K,  ## Num clusters
					vocab,
					nIterations,  ## Num iterations
					alpha,
					eta,
					compute.log.likelihood=TRUE) 
			
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
			
			# get word-weight matrix
			word.proportions <- result$topics_present/(result$topics_present + result$topics_absent + 1e-05)
		} else {
			stop("lda version is invalid...")
		}
		
#		# get top words with probability for each topic 
#		top.words <- matrix(0, nrow=nTopWords, ncol=2*K)
#		for (k in 1:K) {
#			word.distribution <- word.proportions[k,]
#			ordering <- order(word.distribution, decreasing = TRUE)[1:nTopWords]
#			top.words[,2*k-1] <- colnames(word.proportions)[ordering]
#			top.words[,2*k] <- word.distribution[ordering]
#		}
#		colnames(top.words) <- c(sapply(1:K, function(x) c(paste("Topic_",x,"_words",sep=""), paste("Topic_",x,"_prob",sep=""))))

		# word proportions
		word.proportions <- data.frame(word.proportions, check.names = FALSE)
		
		# compute topic popularity
		topic.proportions <- t(result$document_sums) / colSums(result$document_sums)
		topic.popularity <- data.frame(topic=1:K, weight=colSums(topic.proportions, na.rm=TRUE) / nrow(topic.proportions))
		
		# save LDA model
		modelFile <- paste("../../data/model/lda/",query,"_K",K,"_",ldaVersion,".RData",sep="")
		save(query, K, ldaVersion, alpha, eta, nIterations, result, file=modelFile)
		
		# save word proportion of each topic
		weightFile <- paste("../../data/model/lda/",query,"_K",K,"_",ldaVersion,"_weights.tsv",sep="")
		write.table(word.proportions, file=weightFile, quote=FALSE, row.names=FALSE, sep="\t")

		# save topic popularity
		topicPopFile <- paste("../../data/model/lda/",query,"_K",K,"_",ldaVersion,"_topicpop.tsv",sep="")
		write.table(topic.popularity, file=topicPopFile, quote=FALSE, row.names=FALSE, sep="\t")

#		# save top words
#		topWordsFile <- paste("../../data/eBay/result/K10/",query,"_K",K,"_",alpha,"_topwords_",ldaVersion,".csv",sep="")
#		write.csv(top.words, file=topWordsFile)
		
	} # K
} # query



