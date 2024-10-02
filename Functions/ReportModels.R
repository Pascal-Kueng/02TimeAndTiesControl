# Main Functions in this file (3):

###################################
# 1. my_brm
###################################
# Description:
# This function fits a Bayesian regression model using the `brms` package, with the option to include multiple imputed datasets.
# If multiple imputation is specified, the function fits the model across all imputed datasets using parallel processing. 
# The function can save the model object to a file and provides feedback on model convergence based on Rhat values.
#
# Arguments:
# - data: The primary dataset to be used for modeling.
# - imputed_data: A list of imputed datasets (required if mi = TRUE).
# - mi: A logical value indicating whether to perform multiple imputation. Default is FALSE.
# - file: Optional. A string specifying the file path to save the model object.
# - ...: Additional arguments to be passed to the `brm` or `brm_multiple` function.
#
# Example:
# ```r
# model <- my_brm(data = my_data, imputed_data = my_imputed_data, mi = TRUE, file = "my_model")
# ```

###################################
# 1.5 evaluate_model
###################################


###################################
# 2. summarize_brms
###################################
# Description:
# This function summarizes the fixed and random effects from a Bayesian regression model fitted using the `brms` package.
# It can provide a detailed summary or a shortened version, and allows for the exponentiation of coefficients if needed.
# The function also formats the output for easy interpretation, including the addition of significance stars based on credible intervals.
#
# Arguments:
# - model: The fitted Bayesian model object from the `brms` package.
# - short_version: A logical value indicating whether to return a shortened summary. Default is FALSE.
# - exponentiate: A logical value indicating whether to exponentiate the coefficients. Default is FALSE.
# - model_rows_fixed: Optional. A vector specifying which rows of the fixed effects summary to include.
# - model_rows_random: Optional. A vector specifying which rows of the random effects summary to include.
# - model_rownames_fixed: Optional. A vector specifying the row names for the fixed effects summary.
# - model_rownames_random: Optional. A vector specifying the row names for the random effects summary.
#
# Example:
# ```r
# summary <- summarize_brms(model = my_model, short_version = TRUE, exponentiate = TRUE)
# ```

###################################
# 3. report_side_by_side
###################################
# Description:
# This function generates a side-by-side comparison of summaries from multiple Bayesian models fitted using the `brms` package.
# It supports the comparison of both fixed and random effects across models, and can exponentiate coefficients based on the model type.
#
# Arguments:
# - ...: One or more fitted Bayesian model objects from the `brms` package.
# - model_rows_fixed: Optional. A vector specifying which rows of the fixed effects summary to include across all models.
# - model_rows_random: Optional. A vector specifying which rows of the random effects summary to include across all models.
# - model_rownames_fixed: Optional. A vector specifying the row names for the fixed effects summary across all models.
# - model_rownames_random: Optional. A vector specifying the row names for the random effects summary across all models.
#
# Example:
# ```r
# comparison <- report_side_by_side(model1, model2, model_rows_fixed = c("Intercept", "Variable1"))
# ```



# Function for modelling where we can specify whether to use MI or not. 
my_brm <- function(data, imputed_data = NULL, mi = FALSE, file = NULL, ...) {
  if (mi) {
    if (is.null(imputed_data)) {
      stop("Imputed data not provided.")
    }
    if (!is.null(file)) {
      file <- paste0(file, "_imputed")
    }
    library(future)
    plan(multisession, workers = parallel::detectCores(logical = FALSE))
    
    model_object <- brm_multiple(data = imputed_data, file = file, ...)
    print("Rhats of all imputed models range:")
    print(range(model_object$rhats))
    if(max(model_object$rhats) > 1.1) {
      message("Model might not have converged")
    } else {
      message("Model likely converged")
    }
    return(model_object)
  } else {
    return(brm(data = data, file = file, ...))
  }
}


# Functions to facilitate reporting
summarize_brms <- function(model, short_version = FALSE, 
                           exponentiate = FALSE, 
                           model_rows_fixed = NULL, 
                           model_rows_random = NULL,
                           model_rownames_fixed = NULL,
                           model_rownames_random = NULL) {
  
  summ_og <- summary(model)
  summ <- summ_og$fixed
  rands <- summ_og$random[[1]][grep('sd\\(', rownames(summ_og$random[[1]])), ]
  
  ar1sterr <- summ_og$cor_pars
  sdresid <- summ_og$spec_pars
  
  rands <- rbind(rands, ar1sterr, sdresid)
  
  # Get rows that we need
  if (is.null(model_rows_fixed)) {
    model_rows_fixed <- rownames(summ)
  }
  if (is.null(model_rows_random)) {
    model_rows_random <- rownames(rands)
  }
  
  summ <- summ[model_rows_fixed, ]
  rands <- format(round(rands[model_rows_random, ], 2), nsmall = 2)
  
  Rhat <- format(round(summ$Rhat, 3), nsmall = 3)
  
  summ$Est.Error <- format(round(summ$Est.Error, 2), nsmall=2)
  summ$Bulk_ESS <- format(round(summ$Bulk_ESS, 2), nsmall=2)
  summ$Tail_ESS <- format(round(summ$Tail_ESS, 2), nsmall=2)
  
  #exponentiate if needed
  if (exponentiate) {
    summ$Estimate <- exp(summ$Estimate)
  }
  
  # Add stars based on Credible Interval
  significance <- ifelse(
    (summ$`l-95% CI` > 0 & summ$`u-95% CI` > 0) | 
      (summ$`l-95% CI` < 0 & summ$`u-95% CI` < 0), 
    '*', 
    '')
  
  # exponentiate CI
  if (exponentiate) {
    summ$`u-95% CI`<- exp(summ$`u-95% CI`)
    summ$`l-95% CI` <- exp(summ$`l-95% CI`)
  }
  
  
  
  # Round CI
  summ$`l-95% CI` <- format(round(summ$`l-95% CI`, 2), nsmall = 2)
  summ$`u-95% CI` <- format(round(summ$`u-95% CI`, 2), nsmall = 2)
  
  
  summ$Estimate <- ifelse( 
    is.na(summ$Estimate), 
    NA, 
    paste0(
      format(round(summ$Estimate,2), nsmall=2), 
      significance
    )
  )
  
  
  # Rename "Estimate" column correctly
  if (exponentiate) {
    if ('bernoulli' %in% model$family | 'cumulative' %in% model$family) {
      correct_name <- 'OR'
    } else if ('negbinomial' %in% model$family) {
      correct_name <- 'IRR'
    } else {
      correct_name <- 'exp(Est.)'
      warning(
        "Coefficients were exponentiated. Double check if this was intended."
      )
    }
    
  } else {
    correct_name <- 'b'
  }
  
  colnames(summ)[1] <- correct_name
  colnames(rands)[1] <- correct_name
  
  
  # Add rhat rounded to 3 digits back
  summ$Rhat <- Rhat
  
  # Add back all rownames
  if (!is.null(model_rownames_fixed)) {
    rownames(summ) <- model_rownames_fixed
  } else if(!is.null(model_rows_fixed)) {
    rownames(summ) <- model_rows_fixed
  }
  
  if (!is.null(model_rownames_random)) {
    rownames(rands) <- model_rownames_random
  } else if (!is.null(model_rows_random)) {
    rownames(rands) <- model_rows_random
  }
  
  # Add random effects to fixed effects df. 
  fullrep <- rbind(summ, rands)
  
  # if we have multiple imputed datasets, Rhat values are not meaningful, so we remove them
  summ$Rhat <- NA
  
  # We are finished, if we want the full table
  if (!short_version) {
    return(fullrep[ , c(1, 3, 4, 5, 6, 7)])
  }
  
  
  fullrep$`95% CI` <-  ifelse(
    is.na(fullrep[correct_name]) |
      is.na(fullrep$`u-95% CI`) |
      grepl("NA", fullrep$`u-95% CI`) |
      grepl("^\\s*NA\\s*$", fullrep[correct_name]),
    NA,
    paste0("[", fullrep$`l-95% CI`, ", ", fullrep$`u-95% CI`, "]")
  )
  
  return(fullrep[ , c(1, 8)])
}



DHARMa.check_brms <- function(model,        
                       integer = FALSE,   # integer response? (TRUE/FALSE)
                       plot = TRUE,       
                       ...) {
  
  mdata <- brms::standata(model)
  if (!"Y" %in% names(mdata))
    stop("Cannot extract the required information from this brms model")
  
  dharma.obj <- DHARMa::createDHARMa(
    simulatedResponse = t(brms::posterior_predict(model, ndraws = 1000)),
    observedResponse = mdata$Y, 
    fittedPredictedResponse = apply(
      t(brms::posterior_epred(model, ndraws = 1000, re.form = NA)),
      1,
      median),
    integerResponse = integer,
    seed = 123, 
    )
  
  if (isTRUE(plot)) {
    plot(dharma.obj, ...)
  }
  invisible(dharma.obj)
}









DHARMa.check_brms <- function(model,        
                              integer = FALSE,   
                              plot = TRUE,       
                              debug = FALSE,     # Add debug parameter
                              ...) {
  mdata <- brms::standata(model)
  if (!"Y" %in% names(mdata))
    stop("Cannot extract the required information from this brms model")
  
  # Check if model is ordinal
  is_ordinal <- model$family$family %in% c("cumulative", "sratio", "cratio", "acat")
  
  # Get posterior predictions
  simulatedResponse <- brms::posterior_predict(model, ndraws = 1000)
  
  if (is_ordinal) {
    # For ordinal models, get mode of predictions efficiently
    modes <- apply(simulatedResponse, c(1, 2), function(x) {
      tab <- tabulate(x)
      which.max(tab)
    })
    simulatedResponse <- t(modes)
  } else {
    simulatedResponse <- t(simulatedResponse)
  }
  
  # Get fitted values (posterior expectations)
  fittedValues <- if (is_ordinal) {
    epred <- brms::posterior_epred(model, ndraws = 1000, re.form = NA)
    apply(epred, 3, function(x) which.max(colMeans(x)))
  } else {
    apply(t(brms::posterior_epred(model, ndraws = 1000, re.form = NA)), 1, median)
  }
  
  # Ensure all components have the same length
  n <- length(mdata$Y)
  if (nrow(simulatedResponse) != n) {
    stop("Adjusting simulatedResponse to match observed data length")
  }
  if (length(fittedValues) != n) {
    stop("Adjusting fittedValues to match observed data length")
  }
  
  tryCatch({
    dharma.obj <- DHARMa::createDHARMa(
      simulatedResponse = simulatedResponse,
      observedResponse = mdata$Y, 
      fittedPredictedResponse = fittedValues,
      integerResponse = integer,
      seed = 123
    )
    
    if (isTRUE(plot)) {
      plot(dharma.obj, ...)
    }
    
    invisible(dharma.obj)
  }, error = function(e) {
    cat("Error in createDHARMa:\n")
    cat("simulatedResponse dim:", dim(simulatedResponse), "\n")
    cat("observedResponse length:", length(mdata$Y), "\n")
    cat("fittedPredictedResponse length:", length(fittedValues), "\n")
    stop(e)
  })
}













DHARMa.check_brms.all <- function(model, integer = FALSE, ...) {
  model.check <- DHARMa.check_brms(model, integer = integer, plot = FALSE)
  plot(model.check)
  try(testDispersion(model.check))
  try(testZeroInflation(model.check))
  try(testOutliers(model.check))
}


check_brms <- function(model, check_loo = TRUE, integer = FALSE, ...) {
  message('Checking Convergence')
  rstan::check_hmc_diagnostics(model$fit)
  plot(model, ask = FALSE)
  
  message('Checking Residuals')
  DHARMa.check_brms.all(model, integer = integer)
  
  message('Checking Fit')
  plot(pp_check(model, type = 'ecdf_overlay'))
  plot(pp_check(model))
  if (check_loo) {
    loo(model)
    
  }
}



# Function to report all models side by side.
report_side_by_side <- function(..., model_rows_random = NULL, model_rows_fixed = NULL,
                                model_rownames_fixed = NULL, model_rownames_random = NULL) {
  models <- list(...)
  model_names <- sapply(substitute(list(...))[-1], deparse)
  names(models) <- model_names
  
  side_by_side <- NULL
  for (i in seq_along(models)) {
    model <- models[[i]]
    model_name <- model_names[i]
    print(model_name)
    if ('bernoulli' %in% model$family | 'negbinomial' %in% model$family | 'cumulative' %in% model$family | grepl('log', model_name)) {
      exponentiate <- TRUE
    } else {
      exponentiate <- FALSE
    }
    model_summary <- summarize_brms(
      model, short_version = TRUE, 
      exponentiate = exponentiate, 
      model_rows_random = model_rows_random,
      model_rows_fixed = model_rows_fixed,
      model_rownames_fixed = model_rownames_fixed,
      model_rownames_random = model_rownames_random
    )
    
    colnames(model_summary) <- paste(colnames(model_summary), model_name)
    
    if (is.null(side_by_side)) {
      side_by_side <- model_summary
    } else {
      side_by_side <- cbind(side_by_side, model_summary)
    }
  }
  return(side_by_side)
}

