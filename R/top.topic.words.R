
# LDA version of top topic words
lda.top.topic.words <- function (topics, num.words = 20, by.score = FALSE) 
{
    if (by.score) {
        normalized.topics <- topics/(rowSums(topics) + 1e-05)
        scores <- apply(normalized.topics, 2, function(x) x * 
            		(log(x + 1e-05) - sum(log(x + 1e-05))/length(x)))
        apply(scores, 1, function(x) colnames(scores)[order(x, 
            decreasing = TRUE)[1:num.words]])
    }
    else {
        apply(topics, 1, function(x) colnames(topics)[order(x, 
            decreasing = TRUE)[1:num.words]])
    }
}

# MBLDA version of top topic words
mblda.top.topic.words <- function (topics_present, topics_absent, num.words = 20, by.score = FALSE) 
{
	normalized.topics <- topics_present/(topics_present + topics_absent + 1e-05)
	if (by.score) {
		scores <- apply(normalized.topics, 2, function(x) x * 
					(log(x + 1e-05) - sum(log(x + 1e-05))/length(x)))
		apply(scores, 1, function(x) colnames(scores)[order(x, 
			decreasing = TRUE)[1:num.words]])
	}
	else {
		apply(normalized.topics, 1, function(x) colnames(normalized.topics)[order(x, 
			decreasing = TRUE)[1:num.words]])
	}
}
