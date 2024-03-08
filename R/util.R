func_input_type <- function(fun){
  argnames <- formalArgs(fun)
  if ("formula" %in% argnames){
    return("formula")
  } else if ("x" %in% argnames){
    return("xy")
  } else if ("X" %in% argnames){
    return("Xy")
  }
}

find_newname <- function(names_vec){
  name <- "aaa"
  while (name %in% names_vec){
    name <- paste0(name, "a")
  }
  return(name)
}

complete_formula <- function(formula, response_name){
  if (is.null(formula)){
    stop("No formula is found. Please specify a formula ")
  }
  completed_formula <- NULL
  for (i in 1:length(formula)) {
    gformula <- as.character(formula[i])
    gformula <- tail(gformula, 1)
    gformula <- tail(strsplit(gformula, "~")[[1]], 1)
    gformula <- paste0(" ", gformula)
    ## completed_formula <- as.formula(
    ##     paste(response_name, "~", formula),
    ##     env = environment(args$formula))
    completed_formula <- c(completed_formula,paste0(response_name, " ~", gformula))
  }
  return(completed_formula)
}

complete_args <- function(x, response, #fun,
                          args = NULL,
                          weights = NULL){
  # input_type <- func_input_type(fun)
  # if (!input_type %in% c("formula", "xy", "Xy")){
  #   stop("Wrong input type.")
  # }

  response_name <- find_newname(colnames(x))
  
  if(args$type=="gam"){
    # if (input_type == "formula"){
    if (is.null(args) || !"formula" %in% names(args)){
      stop("Formula is not found. Please specify a formula for the fitting function.")
    }
    data <- cbind(data.frame(response), x)
    colnames(data)[1] <- response_name
    args$formula <-  complete_formula(args$formula, response_name)
    data_args <- c(list(data = data), args)
  }
  
  if(args$type=="glmnet" || args$type=="xgboost"){
    data <- cbind(data.frame(response), x)
    colnames(data)[1] <- response_name
    data_args <- c(list(data = data), args)
  }

  data_args <- c(data_args, list(weights = weights))

  return(data_args)
}
