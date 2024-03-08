# modified from adaptMT
safe_gam <- function(formula, family, data, weights = NULL, fit=TRUE, G=NULL, sp=NULL,
                      ...){
    options(warn = -1)

    formula <- as.formula(formula)
    # fit is the param in gam, let res be the results of gam
    
    if (family$link %in% c("inverse", "log")){
        # fit <- try(mgcv::gam(formula, family, data, weights, ...),
        #            silent = TRUE)
        res <- try(mgcv::gam(formula, family, data, weights, fit=fit, G=G, sp=sp, ...),
                   silent = TRUE)
        if (class(res)[1] == "try-error"){
            mod_mat <- model.matrix(formula, data = data)
            p <- ncol(mod_mat) - 1
            start <- c(1, rep(0, p))
            res <- mgcv::gam(formula, family, data, weights, fit=fit, G=G, sp=sp,
                             start = start, ...)
            cat("gam fitting error\n")
        }
    } else {
        res <- mgcv::gam(formula, family, data, weights, fit=fit, G=G, sp=sp, ...)
    }
    
    # for prefit
    
    if (!fit) {
        return(res)
    }

    param <- as.numeric(
        predict(res, type = "response")
        )

    # to change, sp
    info <- list(sp = res$sp)

    options(warn = 0)

    return(list(param = param, info = info))
}

safe_glmnet <- function(x, y, family, weights = NULL,
                        ...){
  options(warn = -1)
  
  if (class(family)[1] == "family"){
    family <- family$family
  }
  
  if (family %in% c("gaussian", "binomial", "poisson", "multinomial", "cox","mgaussian")){
    if (is.null(weights)){
      res <- glmnet::cv.glmnet(x, y,
                               family = family, ...)
    } else {
      weights <- pminmax(weights, 1e-5, 1-1e-5)
      res <- glmnet::cv.glmnet(x, y, weights,
                               family = family, ...)
    }
  } else if (family == "Gamma"){
    if (is.null(weights)){
      res <- HDtweedie::cv.HDtweedie(x, y, p = 2,
                                     standardize = TRUE,
                                     ...)
    } else {
      weights <- pminmax(weights, 1e-5, 1-1e-5)
      res <- HDtweedie::cv.HDtweedie(x, y, p = 2,
                                     weights = weights,
                                     standardize = TRUE,
                                     ...)
    }
  }
  
  param <- as.numeric(
    predict(res, newx = x, s = "lambda.min",
            type = "response")
  )
  
  info <- "empty"
  
  options(warn = 0)
  # Return the model fit for the user to be able to access as well:
  return(list(param = param, info = info))
}


# Function for training xgboost model:
safe_xgboost <- function(data, family, weights = NULL, eta, max_depth, xgb_model=NULL,...) {
  
  options(warn = -1)
  
  # For now we assume that the provided family argument is just a string
  # indicating what type of function to use for the xgboost. This will decide
  # what the objective argument is
  if (class(family)[1] == "family"){
    family <- family$family
  }

  
  # With the way xgboost is set up, each of the possible options will require
  # separate instances, can just use a named vector to access each:
  family_list <- c("gaussian" = "reg:linear",
                   "binomial" = "binary:logistic",
                   "poisson" = "count:poisson",
                   "multinomial" = "multi:softmax",
                   "cox" = "survival:cox",
                   "Gamma" = "reg:gamma")
  
  # If weights are provided then make sure they are within the bounds
  if (!is.null(weights)) {
    weights <- pminmax(weights, 1e-5, 1 - 1e-5)
  }
  
  y=data[,1]
  x=as.matrix(data[,-1])
  
  # Fit the model using the provided family type as the objective with the
  # additional passed in arguments:
  res <- xgboost::xgboost(data = x, label = y, weight = weights, xgb_model=xgb_model,
                          objective = family_list[[family]],nrounds=5,
                          eta=eta,max_depth=max_depth,verbose=0)
  
  # Get the predicted values using the boosting model
  param <- as.numeric(
    predict(res, newdata = x)
  )
  
  info <- res
  
  options(warn = 0)
  # Return the model fit for the user to be able to access as well:
  return(list(param = param, info = info))
}
