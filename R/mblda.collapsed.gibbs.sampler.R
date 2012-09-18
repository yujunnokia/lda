# TODO: lda Collapsed Gibbs Sampler
#
###############################################################################


mblda.collapsed.gibbs.sampler <-
function (documents, K, vocab, num.iterations, alpha, eta, initial = NULL, 
    burnin = NULL, compute.log.likelihood = FALSE, trace = 0L, 
    freeze.topics = FALSE) 
{
    if (class(vocab) == "list") {
        lengths <- as.integer(sapply(vocab, length))
        all.vocab <- do.call(c, vocab)
    }
    else {
        lengths <- as.integer(length(vocab))
        all.vocab <- vocab
    }
    retval <- structure(.Call("collapsedGibbsSamplerMBLDA", documents, 
        as.integer(K), lengths, as.integer(num.iterations), as.double(alpha), 
        as.double(eta), NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
        initial, as.integer(burnin), as.logical(compute.log.likelihood), 
        trace, as.logical(freeze.topics)), names = c("assignments", 
        "topics_present", "topics_absent", "document_sums", if (is.null(burnin)) NA else "document_expects", 
        NA, NA, NA, NA, if (compute.log.likelihood) "log.likelihoods" else NA))
    colnames(retval$topics_present) <- all.vocab
	colnames(retval$topics_absent) <- all.vocab
	
	return(retval)
}
