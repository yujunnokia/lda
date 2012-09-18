#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>


#define CHECK(VAL, TYPE) if (!is##TYPE(VAL)) { \
    error(#VAL " must be a(n) " #TYPE "."); \
  }

#define CHECKLEN(VAL, TYPE, LEN) if (!is##TYPE(VAL) || length(VAL) != LEN) { \
    error(#VAL " must be a length " #LEN " " #TYPE "."); \
  }

#define CHECKMATROW(VAL, TYPE, NROW) if (!isMatrix(VAL) || !is##TYPE(VAL) || NUMROWS(VAL) != NROW) { \
    error(#VAL " must be a matrix with " #NROW " rows of type " #TYPE "."); \
  }


#define NUMROWS(MAT) (INTEGER(GET_DIM(MAT))[0])
#define NUMCOLS(MAT) (INTEGER(GET_DIM(MAT))[1])

#define UPDATESUMS(weight) { \
  INTEGER(topics)[z + sumK * word] += weight * count; \
  INTEGER(topic_sums)[z] += weight * count; \
  int src = INTEGER(source)[c];	\
  int lt = INTEGER(local_topics)[src]; \
  INTEGER(VECTOR_ELT(document_sums, src))[z - partialK[lt]] += weight * count; \
  INTEGER(document_sources)[c] += weight * count; \
  document_lengths[src] += weight * count; \
  }


#define UPDATERTMSUMS(weight) { \
  if (dd < test_start) { \
    INTEGER(topics)[z + K * word] += weight * count; \
    INTEGER(topic_sums)[z] += weight * count; \
  } \
  INTEGER(document_sums)[z + K * dd] += weight * count; \
  }

//-----------------------
// LDA Gibbs sampler.
//-----------------------
SEXP collapsedGibbsSamplerLDA(SEXP documents,  // document: a list of 2 by N matrix
			   SEXP K_,  // number of topics
			   SEXP V_,  // size of vacabulary
			   SEXP N_,  // number of gibbs sampling iterations
			   SEXP alpha_,  // alpha hyper-parameters
			   SEXP eta_,  // eta hyper-parameters
			   SEXP annotations,  // NULL for LDA
			   SEXP beta,  // NULL for LDA
			   SEXP var_,  // NULL for LDA
			   SEXP method_,  // NULL for LDA
			   SEXP lambda_,  // NULL for LDA
			   SEXP nbeta,  // NULL for LDA
			   SEXP net_annotations,  // NULL for LDA
			   SEXP initial_,  // initialization
			   SEXP burnin_,  // a scalar integer indicating the number of Gibbs sweeps as burn-in (i.e., throw away) 
			   SEXP compute_log_likelihood_,  // a scalar logical causing the sampler to compute the log likelihood of the words. Can be helpful to evaluate convergence.
			   SEXP trace_, // output error messages
			   SEXP freeze_topics_) { // FALSE
  GetRNGstate();
  int dd;
  int ii;
  int kk;
  double var = 0;
  int logistic = 0;
  double lambda = 0;
  int burnin = -1;

  CHECKLEN(alpha_, Real, 1);
  double alpha = REAL(alpha_)[0];

  CHECKLEN(eta_, Real, 1);
  double eta = REAL(eta_)[0];

  CHECKLEN(K_, Integer, 1);
  int K = INTEGER(K_)[0];

  CHECKLEN(trace_, Integer, 1);
  int trace = INTEGER(trace_)[0];

  CHECK(V_, Integer);
  int V = 0;
  for (ii = 0; ii < length(V_); ++ii) {
    V += INTEGER(V_)[ii];
  }

  CHECKLEN(N_, Integer, 1);
  int N = INTEGER(N_)[0];

  int method = -1;

  CHECK(documents, NewList);
  int nd = length(documents);

  double* dv = NULL;
  double* wx2 = NULL;
  double* wx = NULL;

  if (!isNull(beta)) {
	error("beta must be null when annotations are empty.");
  }
  if (!isNull(var_)) {
	error("var must be null when annotations are empty.");
  }

  SEXP retval;
  PROTECT(retval = allocVector(VECSXP, 10));

  SEXP assignments;
  SEXP topics = NULL;
  SEXP topic_sums = NULL;
  SEXP document_expects = NULL;
  SEXP document_sums;
  SEXP initial = NULL;
  SEXP initial_topic_sums = NULL;
  SEXP initial_topics = NULL;
  SEXP log_likelihood = NULL;

  // allocate space to store return objects
  SET_VECTOR_ELT(retval, 0, assignments = allocVector(VECSXP, nd));
  SET_VECTOR_ELT(retval, 1, topics = allocMatrix(INTSXP, K, V));
  SET_VECTOR_ELT(retval, 2, topic_sums = allocMatrix(INTSXP, K, length(V_)));
  SET_VECTOR_ELT(retval, 3, document_sums = allocMatrix(INTSXP, K, nd));

  CHECKLEN(compute_log_likelihood_, Logical, 1);
  int compute_log_likelihood = LOGICAL(compute_log_likelihood_)[0];
  CHECKLEN(freeze_topics_, Logical, 1);
  int freeze_topics = LOGICAL(freeze_topics_)[0];
  if (compute_log_likelihood) {
    SET_VECTOR_ELT(retval, 9, log_likelihood = allocMatrix(REALSXP, 2, N));
  }
  
  if (length(burnin_) > 0) {
    CHECKLEN(burnin_, Integer, 1);
    burnin = INTEGER(burnin_)[0];
    if (burnin < 0) {
      error("burnin must be positive.");
    }

    SET_VECTOR_ELT(retval, 4, document_expects = allocMatrix(INTSXP, K, nd));
    for (ii = 0; ii < K * nd; ++ii) {
      INTEGER(document_expects)[ii] = 0;
    }
  }  // if for burnin

  // if initialized
  if (!isNull(initial_)) {
    CHECK(initial_, NewList);
    SEXP names = getAttrib(initial_, R_NamesSymbol);

    for (ii = 0; ii < length(initial_); ++ii) {
      if (!strcmp(CHAR(STRING_ELT(names, ii)), "assignments")) {
		initial = VECTOR_ELT(initial_, ii);
		CHECKLEN(initial, NewList, nd);
      } else if (!strcmp(CHAR(STRING_ELT(names, ii)), "topic_sums")) {
		initial_topic_sums = VECTOR_ELT(initial_, ii);
		if (!isInteger(initial_topic_sums) ||
			INTEGER(GET_DIM(initial_topic_sums))[0] != K ||
			INTEGER(GET_DIM(initial_topic_sums))[1] != length(V_)) {
		  error("Initial topic sums must be a K x length(V) integer matrix.");
		}
      } else if (!strcmp(CHAR(STRING_ELT(names, ii)), "topics")) {
		initial_topics = VECTOR_ELT(initial_, ii);
		if (!isInteger(initial_topics) ||
			INTEGER(GET_DIM(initial_topics))[0] != K ||
			INTEGER(GET_DIM(initial_topics))[1] != V) {
		  error("Initial topics (%d x %d) must be a %d x %d integer matrix.",
			INTEGER(GET_DIM(initial_topics))[0],
			INTEGER(GET_DIM(initial_topics))[1],
			K,
			V);
		}
      } else {
		error("Unrecognized initialization field: '%s'",
	      CHAR(STRING_ELT(names, ii)));
      }
    }
  } // end if for initial

  if ((initial_topic_sums == NULL) ^ (initial_topics == NULL)) {
    error("initial topic sums and topics must both be specified.");
  }

  // set topics to be zero matrix
  if (initial_topics == NULL) {
    for (ii = 0; ii < K * V; ++ii) {
      INTEGER(topics)[ii] = 0;
    }
  } else {
    for (ii = 0; ii < K * V; ++ii) {
      INTEGER(topics)[ii] = INTEGER(initial_topics)[ii];
    }
  }
  
  // set topic_sums to be zero matrix
  if (initial_topic_sums == NULL) {
    for (ii = 0; ii < K * length(V_); ++ii) {
      INTEGER(topic_sums)[ii] = 0;
    }
  } else {
    for (ii = 0; ii < K * length(V_); ++ii) {
      INTEGER(topic_sums)[ii] = INTEGER(initial_topic_sums)[ii];
    }
  }

  // set document_sums to be zero matrix
  for (ii = 0; ii < K * nd; ++ii) {
    INTEGER(document_sums)[ii] = 0;
  }

  // for each document, set topic assignment to be -1
  for (dd = 0; dd < nd; ++dd) {
    int ww;
    SEXP document = VECTOR_ELT(documents, dd);

    CHECKMATROW(document, Integer, 2);

	// get the number of words per doc
    int nw = INTEGER(GET_DIM(document))[1];
    SET_VECTOR_ELT(assignments, dd, allocVector(INTSXP, nw));
    SEXP zs = VECTOR_ELT(assignments, dd);

    for (ww = 0; ww < nw; ++ww) {
      int word = INTEGER(document)[ww * 2];
      int count = INTEGER(document)[ww * 2 + 1];
      if (count < 0) {
		error("Count must be positive.");
      }
      if (word >= V || word < 0) {
		error("Word (%d) must be positive and less than or "
	      "equal to the number of words (%d).", word, V);
      }
      INTEGER(zs)[ww] = -1;
    }
  }

  double* p = (double *)R_alloc(K, sizeof(double));

  double const_prior = 0;
  double const_ll = 0;

  // compute log likelihood
  if (compute_log_likelihood) {
    //                log B(\alpha)
    const_prior = (K * lgammafn(alpha) - lgammafn(alpha * K)) * nd;
    //                log B(\eta)
    const_ll = (V * lgammafn(eta) - lgammafn(eta * V)) * K;
  }

  //----------------------------------
  // start gibbs sampling iterations
  //----------------------------------
  int iteration;
  for (iteration = 0; iteration < N; ++iteration) {
    if (trace >= 1) {
      printf("Iteration %d\n", iteration);
    }
	if (iteration % 10 == 0) {
      printf("Iteration %d\n", iteration);
    }
	
	// for each doc
    for (dd = 0; dd < nd; ++dd) {
      R_CheckUserInterrupt();
      SEXP zs = VECTOR_ELT(assignments, dd); // get current topic assignments
      SEXP document = VECTOR_ELT(documents, dd);  // get documents
      int ww;
      int nw = INTEGER(GET_DIM(document))[1];  // number of word in this doc
      SEXP initial_d = NULL;

      if (initial != NULL) {
		initial_d = VECTOR_ELT(initial, dd);
		CHECKLEN(initial_d, Integer, nw);
      }

	  // for each word in the document
	  for (ww = 0; ww < nw; ++ww) {
		int* z = &INTEGER(zs)[ww]; // get topic assignment
		int word = -1;
		int count = 1;
		int* topic_wk;
		int* topic_k;
		int* document_k;
  
		word = INTEGER(document)[ww * 2];
		int partialsum = 0;
		int topic_index = -1;
		for (ii = 0; ii < length(V_); ++ii) {
		  partialsum += INTEGER(V_)[ii];
		  if (word < partialsum) {
			topic_index = ii;
		  }
		}
		if (topic_index == -1) {
		  error("Oops I did it again");
		}
		count = INTEGER(document)[ww * 2 + 1];
	
		// reduce the count from topics, topic_sums and document_sums
		// for the given topic
		if (*z != -1) {
		  topic_wk = &INTEGER(topics)[(*z) + K * word];
		  topic_k = &INTEGER(topic_sums)[*z + K * topic_index];
		  if(!freeze_topics)
		  {
			*topic_wk -= count;
			*topic_k -= count;
		  }
		  document_k = &INTEGER(document_sums)[K * dd + *z];
		  *document_k -= count;
	
		  if (*topic_wk < 0 || *topic_k < 0 || *document_k < 0) {
			error("Counts became negative for word (%d): (%d, %d, %d)",
			  word, *topic_wk, *topic_k, *document_k);
		  }
		}
	
		// compute the conditional probability
		double r = unif_rand();
		double p_sum = 0.0;
		for (kk = 0; kk < K; ++kk) {
		  if (*z == -1) {
			if (initial != NULL) {
			  if (INTEGER(initial_d)[ww] == kk) {
				p[kk] = 1;
			  } else {
				p[kk] = 0;
			  }
			} else {
			  p[kk] = 1;
			}
		  } else {
			p[kk] = (INTEGER(document_sums)[K * dd + kk] + alpha);
			p[kk] *= (INTEGER(topics)[kk + K * word] + eta);
			p[kk] /= (INTEGER(topic_sums)[kk + K * topic_index] + V * eta);
		  }
		  p_sum += p[kk];
		}
	
		if (p_sum <= 0.0) {
		  kk = K - 1;
		  error("Numerical problems.");
		}
	
		// select a topic according to the conditional probablity
		*z = -1;
		for (kk = 0; kk < K; ++kk) {
		  if (r < p[kk] / p_sum) {
			*z = kk;
			break;
		  }
		  r -= p[kk] / p_sum;
		}
	
		if (*z == -1) {
		  for (kk = 0; kk < K; ++kk) {
			printf("%g\n", p[kk]);
		  }
		  error("This should not have happened (%g).", r);
		}
	
		if(!freeze_topics)
		{
		  INTEGER(topics)[*z + K * word] += count;
		  INTEGER(topic_sums)[*z + K * topic_index] += count;
		}
		INTEGER(document_sums)[K * dd + *z] += count;
		if (burnin > -1 && iteration >= burnin) {
		  INTEGER(document_expects)[K * dd + *z] += count;
		}
	  } // end for each word
	} // end for each doc

    /*Compute the likelihoods:*/
    if (compute_log_likelihood) {
      double doc_ll = 0;
      for (dd = 0; dd < nd; ++dd) {
		double sum = alpha * K;
		for (kk = 0; kk < K; ++kk) {
		  doc_ll += lgammafn(INTEGER(document_sums)[K * dd + kk] + alpha);
		  sum += INTEGER(document_sums)[K * dd + kk];
		}
		doc_ll -= lgammafn(sum);
      }
      double topic_ll = 0;
      for (kk = 0; kk < K; ++kk) {
		double sum = eta * V;
		for (ii = 0; ii < V; ++ii) {
		  topic_ll += lgammafn(INTEGER(topics)[kk + K * ii] + eta);
		  sum += INTEGER(topics)[kk + K * ii];
		}
		topic_ll -= lgammafn(sum);
      }
      if (trace >= 2) {
		printf("ll: %g + %g - %g - %g = %g\n", doc_ll, topic_ll, const_ll, const_prior,
	       doc_ll + topic_ll - const_ll - const_prior);
      }
      REAL(log_likelihood)[2 * iteration] = doc_ll - const_prior + topic_ll - const_ll;
      REAL(log_likelihood)[2 * iteration + 1] = topic_ll - const_ll;
    } // if for compute likelihood
  } // for iteration


  PutRNGstate();
  UNPROTECT(1);
  return retval;
}

//-----------------------
// MBLDA Gibbs sampler.
//-----------------------
SEXP collapsedGibbsSamplerMBLDA(SEXP documents,  // document: a list of 2 by N matrix
			   SEXP K_,  // number of topics
			   SEXP V_,  // size of vacabulary
			   SEXP N_,  // number of gibbs sampling iterations
			   SEXP alpha_,  // alpha hyper-parameters
			   SEXP eta_,  // eta hyper-parameters
			   SEXP annotations,  // NULL for LDA
			   SEXP beta,  // NULL for LDA
			   SEXP var_,  // NULL for LDA
			   SEXP method_,  // NULL for LDA
			   SEXP lambda_,  // NULL for LDA
			   SEXP nbeta,  // NULL for LDA
			   SEXP net_annotations,  // NULL for LDA
			   SEXP initial_,  // initialization
			   SEXP burnin_,  // a scalar integer indicating the number of Gibbs sweeps as burn-in (i.e., throw away) 
			   SEXP compute_log_likelihood_,  // a scalar logical causing the sampler to compute the log likelihood of the words. Can be helpful to evaluate convergence.
			   SEXP trace_, // output error messages
			   SEXP freeze_topics_) { // FALSE
  GetRNGstate();
  int dd;
  int ii;
  int kk;
  double var = 0;
  int logistic = 0;
  double lambda = 0;
  int burnin = -1;

  CHECKLEN(alpha_, Real, 1);
  double alpha = REAL(alpha_)[0];

  CHECKLEN(eta_, Real, 1);
  double eta = REAL(eta_)[0];

  CHECKLEN(K_, Integer, 1);
  int K = INTEGER(K_)[0];

  CHECKLEN(trace_, Integer, 1);
  int trace = INTEGER(trace_)[0];

  CHECK(V_, Integer);
  int V = 0;
  for (ii = 0; ii < length(V_); ++ii) {
    V += INTEGER(V_)[ii];
  }

  CHECKLEN(N_, Integer, 1);
  int N = INTEGER(N_)[0];

  int method = -1;

  CHECK(documents, NewList);
  int nd = length(documents);

  double* dv = NULL;
  double* wx2 = NULL;
  double* wx = NULL;

  if (!isNull(beta)) {
	error("beta must be null when annotations are empty.");
  }
  if (!isNull(var_)) {
	error("var must be null when annotations are empty.");
  }

  SEXP retval;
  PROTECT(retval = allocVector(VECSXP, 10));

  SEXP assignments;
  SEXP topics_present = NULL;
  SEXP topics_absent = NULL;
  SEXP document_expects = NULL;
  SEXP document_sums;
  SEXP initial = NULL;
  SEXP initial_topic_sums = NULL;
  SEXP initial_topics = NULL;
  SEXP log_likelihood = NULL;

  // allocate space to store return objects
  SET_VECTOR_ELT(retval, 0, assignments = allocVector(VECSXP, nd));
  SET_VECTOR_ELT(retval, 1, topics_present = allocMatrix(INTSXP, K, V));
  SET_VECTOR_ELT(retval, 2, topics_absent  = allocMatrix(INTSXP, K, V));
  SET_VECTOR_ELT(retval, 3, document_sums  = allocMatrix(INTSXP, K, nd));

  CHECKLEN(compute_log_likelihood_, Logical, 1);
  int compute_log_likelihood = LOGICAL(compute_log_likelihood_)[0];
  CHECKLEN(freeze_topics_, Logical, 1);
  int freeze_topics = LOGICAL(freeze_topics_)[0];
  if (compute_log_likelihood) {
    SET_VECTOR_ELT(retval, 9, log_likelihood = allocMatrix(REALSXP, 2, N));
  }
  
  // if for burnin
  if (length(burnin_) > 0) {
    CHECKLEN(burnin_, Integer, 1);
    burnin = INTEGER(burnin_)[0];
    if (burnin < 0) {
      error("burnin must be positive.");
    }

    SET_VECTOR_ELT(retval, 4, document_expects = allocMatrix(INTSXP, K, nd));
    for (ii = 0; ii < K * nd; ++ii) {
      INTEGER(document_expects)[ii] = 0;
    }
  }  

  if ((initial_topic_sums == NULL) ^ (initial_topics == NULL)) {
    error("initial topic sums and topics must both be specified.");
  }

  // set topics to be zero matrix
  if (initial_topics == NULL) {
    for (ii = 0; ii < K * V; ++ii) {
      INTEGER(topics_present)[ii] = 0;
	  INTEGER(topics_absent)[ii] = 0;
    }
  } else {
    for (ii = 0; ii < K * V; ++ii) {
      INTEGER(topics_present)[ii] = INTEGER(initial_topics)[ii];
	  INTEGER(topics_absent)[ii] = INTEGER(initial_topics)[ii];
    }
  }
  
  
  // set document_sums to be zero matrix
  for (ii = 0; ii < K * nd; ++ii) {
    INTEGER(document_sums)[ii] = 0;
  }

  // for each document, set topic assignment to be -1
  for (dd = 0; dd < nd; ++dd) {
    int ww;
    SEXP document = VECTOR_ELT(documents, dd);

    CHECKMATROW(document, Integer, 2);

	// get the number of words per doc
    int nw = INTEGER(GET_DIM(document))[1];
    SET_VECTOR_ELT(assignments, dd, allocVector(INTSXP, nw));
    SEXP zs = VECTOR_ELT(assignments, dd);

	// for each word, set the assignment to be -1
    for (ww = 0; ww < nw; ++ww) {
      int word = INTEGER(document)[ww * 2];
      int count = INTEGER(document)[ww * 2 + 1];
      if (count < 0) {
		error("Count must be positive.");
      }
      if (word >= V || word < 0) {
		error("Word (%d) must be positive and less than or "
	      "equal to the number of words (%d).", word, V);
      }
      INTEGER(zs)[ww] = -1;
    }
  }

  double* p = (double *)R_alloc(K, sizeof(double));

  double const_prior = 0;
  double const_ll = 0;

  // compute log likelihood
  if (compute_log_likelihood) {
    //                log B(\alpha)
    const_prior = (K * lgammafn(alpha) - lgammafn(alpha * K)) * nd;
    //                log B(\eta)
    const_ll = (V * lgammafn(eta) - lgammafn(eta * V)) * K;
  }

  //----------------------------------
  // start gibbs sampling iterations
  //----------------------------------
  int iteration;
  for (iteration = 0; iteration < N; ++iteration) {
    if (trace >= 1) {
      printf("Iteration %d\n", iteration);
    }
    if (iteration % 10 == 0) {
      printf("Iteration %d\n", iteration);
    }
	
	// for each doc
    for (dd = 0; dd < nd; ++dd) {
      R_CheckUserInterrupt();
      SEXP zs = VECTOR_ELT(assignments, dd); // get current topic assignments
      SEXP document = VECTOR_ELT(documents, dd);  // get documents
      int ww;
      int nw = INTEGER(GET_DIM(document))[1];  // number of word in this doc
      SEXP initial_d = NULL;

      if (initial != NULL) {
		initial_d = VECTOR_ELT(initial, dd);
		CHECKLEN(initial_d, Integer, nw);
      }
	  
	  //printf("doc %d out of %d\n", dd, nd);
	  
	  // for each word in the document
	  for (ww = 0; ww < nw; ++ww) {
		int* z = &INTEGER(zs)[ww]; // get topic assignment
		int word = -1;
		int count = 1;
		int* topic_wk;
		int* document_k;
  
		word = INTEGER(document)[ww * 2];
		count = INTEGER(document)[ww * 2 + 1];
	
		// reduce the count from topics, topic_sums and document_sums
		// for the given topic
		if (*z != -1) {
		  if (count == 1) {
			topic_wk = &INTEGER(topics_present)[(*z) + K * word];
		  } else {
			topic_wk = &INTEGER(topics_absent)[(*z) + K * word];		  
		  }
		  if(!freeze_topics)
		  {
			*topic_wk -= 1;
		  }
		  document_k = &INTEGER(document_sums)[K * dd + *z];
		  *document_k -= 1;
	
		  if (*topic_wk < 0 || *document_k < 0) {
			error("Counts became negative for word (%d): (%d, %d)",
			  word, *topic_wk, *document_k);
		  }
		}
	
		// compute the conditional probability
		double r = unif_rand();
		double p_sum = 0.0;
		for (kk = 0; kk < K; ++kk) {
		  if (*z == -1) {
			if (initial != NULL) {
			  if (INTEGER(initial_d)[ww] == kk) {
				p[kk] = 1;
			  } else {
				p[kk] = 0;
			  }
			} else {
			  p[kk] = 1;
			}
		  } else {
			//p[kk] = (INTEGER(document_sums)[K * dd + kk] + alpha);
			//p[kk] *= (INTEGER(topics)[kk + K * word] + eta);
			//p[kk] /= (INTEGER(topic_sums)[kk + K * topic_index] + V * eta);
			p[kk] = (INTEGER(document_sums)[K * dd + kk] + alpha);
			if (count == 1) {
			  p[kk] *= (INTEGER(topics_present)[kk + K * word] + eta);
			} else {
			  p[kk] *= (INTEGER(topics_absent)[kk + K * word] + eta);
			}
			p[kk] /= (INTEGER(topics_present)[kk + K * word] + INTEGER(topics_absent)[kk + K * word] + 2 * eta);
		  }
		  p_sum += p[kk];
		}
	
		if (p_sum <= 0.0) {
		  kk = K - 1;
		  error("Numerical problems.");
		}
	
		// select a topic according to the conditional probablity
		*z = -1;
		for (kk = 0; kk < K; ++kk) {
		  if (r < p[kk] / p_sum) {
			*z = kk;
			break;
		  }
		  r -= p[kk] / p_sum;
		}
	
		if (*z == -1) {
		  for (kk = 0; kk < K; ++kk) {
			printf("%g\n", p[kk]);
		  }
		  error("This should not have happened (%g).", r);
		}
	
		// update statistics
		if(!freeze_topics)
		{
		  if (count == 1) {
			INTEGER(topics_present)[*z + K * word] += 1;
		  } else {
			INTEGER(topics_absent)[*z + K * word] += 1;
		  }
		}
		INTEGER(document_sums)[K * dd + *z] += 1;
		if (burnin > -1 && iteration >= burnin) {
		  INTEGER(document_expects)[K * dd + *z] += 1;
		}
	  } // end for each word
	} // end for each doc

    /*Compute the likelihoods:*/
    if (compute_log_likelihood) {
        /*
      double doc_ll = 0;
      for (dd = 0; dd < nd; ++dd) {
		double sum = alpha * K;
		for (kk = 0; kk < K; ++kk) {
		  doc_ll += lgammafn(INTEGER(document_sums)[K * dd + kk] + alpha);
		  sum += INTEGER(document_sums)[K * dd + kk];
		}
		doc_ll -= lgammafn(sum);
      }
      double topic_ll = 0;
      for (kk = 0; kk < K; ++kk) {
		double sum = eta * V;
		for (ii = 0; ii < V; ++ii) {
		  topic_ll += lgammafn(INTEGER(topics)[kk + K * ii] + eta);
		  sum += INTEGER(topics)[kk + K * ii];
		}
		topic_ll -= lgammafn(sum);
      }
      if (trace >= 2) {
		printf("ll: %g + %g - %g - %g = %g\n", doc_ll, topic_ll, const_ll, const_prior,
	       doc_ll + topic_ll - const_ll - const_prior);
      }
      REAL(log_likelihood)[2 * iteration] = doc_ll - const_prior + topic_ll - const_ll;
      REAL(log_likelihood)[2 * iteration + 1] = topic_ll - const_ll;
         */
    } // if for compute likelihood
  } // for iteration


  PutRNGstate();
  UNPROTECT(1);
  return retval;
}

