#modified from adaptMT

model_selection_cv <- function(x,q,indReveal,models,gamma) {
  initial=FALSE
  nfolds=5
  n <- length(q)
  mas=(1:n)[!indReveal]
  m <- length(models)
  cat(paste0("Model selection with ", nfolds, "-fold cross-validation starts!\n"))
  cat("Shrink the set of candidate models or number of folds if it is too time-consuming.")
  cat("\n")
  pb <- txtProgressBar(min = 0, max = m, style = 3, width = 50)
  #cat("\n")

  k_cv_holdout_i <- caret::createFolds(mas, k = nfolds)

  model_results <- sapply(1:m,
                          function(model_i) {
                            model <- models[[model_i]]
                            # Next loop through the folds:
                            model_cv_loglik <- sapply(1:nfolds,
                                                      function(fold_i) {
                                                        x_holdout_fold_i=x[-k_cv_holdout_i[[fold_i]],]
                                                        indReveal_holdout_fold_i=indReveal[-k_cv_holdout_i[[fold_i]]]
                                                        q_holdout_fold_i=q[-k_cv_holdout_i[[fold_i]]]
                                                        param_holdout_fold_i=initialparameter(q_holdout_fold_i,indReveal_holdout_fold_i,x_holdout_fold_i,model)
                                                        fit <- try(
                                                          param_holdout_fold_i <- modelfitting(param_holdout_fold_i,
                                                                                               x_holdout_fold_i,
                                                                                               q_holdout_fold_i,
                                                                                               indReveal_holdout_fold_i,
                                                                                               model,initial,gamma,ms=TRUE),
                                                          silent = TRUE
                                                        )
                                                        if (class(fit)[1] == "try-error"){
                                                          warning(paste0("Model ", model_i, "with fold ",
                                                                         fold_i," fails."))
                                                          NA
                                                        } else {
                                                          # Use the fitted parameters on the training
                                                          # data to compute the expected loglikelihood
                                                          # on the holdout p-values.

                                                          logl_ms(fit,x[k_cv_holdout_i[[fold_i]],],
                                                                  q[k_cv_holdout_i[[fold_i]]],
                                                                  indReveal[k_cv_holdout_i[[fold_i]]],gamma)

                                                        }
                                                      })
                            # Return sum of holdout log-likelihood:
                            loglik_sum <- sum(unlist(model_cv_loglik), na.rm = TRUE)
                            setTxtProgressBar(pb, model_i)
                            #cat("\n")

                            return(loglik_sum)
                          })
  model_results <- unlist(model_results)
  cat("\n")

  if (all(is.na(model_results)) | all(is.infinite(model_results))) {
    stop("All models fail.")
  }

  # Otherwise which was the best model:
  best_model_i <- which.max(model_results)
  cat(paste0("Selected model parameter choice: ", best_model_i, "!"))
  cat("\n")

  best_model <- models[[best_model_i]]


  return(list(model = best_model,
              loglinfo = model_results))

}
