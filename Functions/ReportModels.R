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




# New untested version:
summarize_brms <- function(model, 
                           short_version = FALSE,
                           exponentiate = FALSE,
                           model_rows_fixed = NULL,
                           model_rows_random = NULL,
                           model_rownames_fixed = NULL,
                           model_rownames_random = NULL) {
  
  # Extract summaries
  summ_og <- summary(model)
  fixed_effects <- summ_og$fixed
  random_effects <- summ_og$random[[1]][grep('sd\\(', rownames(summ_og$random[[1]])), ]
  random_effects <- rbind(random_effects, summ_og$cor_pars, summ_og$spec_pars)
  
  # Add p_direction to fixed effects
  p_dir <- as.data.frame(bayestestR::p_direction(
    model,
    effects = 'fixed', 
    component = 'conditional'
  ))
  
  p_dir <- p_dir$pd
  
  if (length(p_dir) != nrow(fixed_effects)) {
    stop("Number of variables in p_direction and fixed_effects do not match.")
  } 
  
  fixed_effects$p_direction <- round(p_dir, 3)
  random_effects$p_direction <- NA
  
  
  
  
  
  
  # Select rows
  model_rows_fixed <- model_rows_fixed %||% rownames(fixed_effects)
  model_rows_random <- model_rows_random %||% rownames(random_effects)
  
  fixed_effects <- fixed_effects[model_rows_fixed, ]
  random_effects <- format(round(random_effects[model_rows_random, ], 2), nsmall = 2)
  
  # Format fixed effects
  format_number <- function(x, digits = 2) format(round(x, digits), nsmall = digits)
  
  fixed_effects$Est.Error <- format_number(fixed_effects$Est.Error)
  fixed_effects$Bulk_ESS <- format_number(fixed_effects$Bulk_ESS)
  fixed_effects$Tail_ESS <- format_number(fixed_effects$Tail_ESS)
  Rhat <- format_number(fixed_effects$Rhat, 3)
  
  # Add significance stars
  is_significant <- function(low, high) {
    (low > 0 & high > 0) | (low < 0 & high < 0)
  }
  significance <- ifelse(
    is_significant(fixed_effects$`l-95% CI`, fixed_effects$`u-95% CI`),
    '*', 
    ''
  )
  
  # Handle exponentiation
  if (exponentiate) {
    fixed_effects$Estimate <- exp(fixed_effects$Estimate)
    fixed_effects$`u-95% CI` <- exp(fixed_effects$`u-95% CI`)
    fixed_effects$`l-95% CI` <- exp(fixed_effects$`l-95% CI`)
  }
  
  # Format CIs
  fixed_effects$`l-95% CI` <- format_number(fixed_effects$`l-95% CI`)
  fixed_effects$`u-95% CI` <- format_number(fixed_effects$`u-95% CI`)
  
  # Format estimates with significance
  fixed_effects$Estimate <- ifelse(
    is.na(fixed_effects$Estimate),
    NA,
    paste0(format_number(fixed_effects$Estimate), significance)
  )
  
  # Determine correct column name
  correct_name <- if (exponentiate) {
    if (model$family[[1]] %in% c('bernoulli', 'cumulative')) {
      'OR'
    } else if (model$family[[1]] == 'negbinomial') {
      'IRR'
    } else {
      warning("Coefficients were exponentiated. Double check if this was intended.")
      'exp(Est.)'
    }
  } else {
    'b'
  }
  
  # Rename columns
  colnames(fixed_effects)[1] <- colnames(random_effects)[1] <- correct_name
  
  # Handle row names
  if (!is.null(model_rownames_fixed)) rownames(fixed_effects) <- model_rownames_fixed
  if (!is.null(model_rownames_random)) rownames(random_effects) <- model_rownames_random
  
  # Combine fixed and random effects
  fixed_effects$Rhat <- Rhat
  full_results <- rbind(fixed_effects, random_effects)
  
  if (!short_version) {
    return(full_results[, c(1, 3:8)])
  }
  
  # Create short version with CI
  full_results$`95% CI` <- ifelse(
    is.na(full_results[[correct_name]]) | 
      is.na(full_results$`u-95% CI`) | 
      grepl("NA", full_results$`u-95% CI`) | 
      grepl("^\\s*NA\\s*$", full_results[[correct_name]]),
    NA,
    paste0("[", full_results$`l-95% CI`, ", ", full_results$`u-95% CI`, "]")
  )
  
  return(full_results[, c(1, 8)])
}


#### untested new version
# Updated function to include Bayes Factor column
summarize_brms <- function(model, 
                           short_version = FALSE,
                           exponentiate = FALSE,
                           model_rows_fixed = NULL,
                           model_rows_random = NULL,
                           model_rownames_fixed = NULL,
                           model_rownames_random = NULL,
                           plot_which_bayes_factor = c()) {
  
  # Extract summaries
  summ_og <- summary(model)
  fixed_effects <- summ_og$fixed
  random_effects <- summ_og$random[[1]][grep('sd\\(', rownames(summ_og$random[[1]])), ]
  random_effects <- rbind(random_effects, summ_og$cor_pars, summ_og$spec_pars)
  
  # Add p_direction to fixed effects
  p_dir <- as.data.frame(bayestestR::p_direction(
    model,
    effects = 'fixed', 
    component = 'conditional'
  ))
  
  p_dir <- p_dir$pd
  
  if (length(p_dir) != nrow(fixed_effects)) {
    stop("Number of variables in p_direction and fixed_effects do not match.")
  } 
  
  fixed_effects$p_direction <- format(round(p_dir, 3), nsmall = 2)
  random_effects$p_direction <- NA
  
  # Calculate Bayes Factor for fixed effects
  bayesfac <- bayestestR::bayesfactor(
    model,
    effects = "fixed"
  )
  fixed_effects$Bayes_Factor <- ifelse(
    exp(bayesfac$log_BF) > 100, 
    '>100', 
    sprintf("%.3f", exp(bayesfac$log_BF))
  )
  
  # Add evidence interpretation using case_when for clarity
  fixed_effects <- fixed_effects %>%
    mutate(
      Evidence = case_when(
        exp(bayesfac$log_BF) > 100          ~ "Overwhelming Evidence",
        exp(bayesfac$log_BF) > 30           ~ "Very Strong Evidence",
        exp(bayesfac$log_BF) > 10           ~ "Strong Evidence",
        exp(bayesfac$log_BF) > 3            ~ "Moderate Evidence",
        exp(bayesfac$log_BF) > 1            ~ "Weak Evidence",
        exp(bayesfac$log_BF) > 0.3          ~ "Weak Evidence for Null",
        exp(bayesfac$log_BF) > 0.1          ~ "Moderate Evidence for Null",
        exp(bayesfac$log_BF) > 0.03         ~ "Strong Evidence for Null",
        TRUE                                ~ "Very Strong Evidence for Null"
      )
    )
  
  random_effects$Bayes_Factor <- NA
  random_effects$Evidence <- NA
  # Select rows
  model_rows_fixed <- model_rows_fixed %||% rownames(fixed_effects)
  model_rows_random <- model_rows_random %||% rownames(random_effects)
  
  fixed_effects <- fixed_effects[model_rows_fixed, ]
  random_effects <- format(round(random_effects[model_rows_random, ], 2), nsmall = 2)
  
  # Format fixed effects
  format_number <- function(x, digits = 2) format(round(x, digits), nsmall = digits)
  
  fixed_effects$Est.Error <- format_number(fixed_effects$Est.Error)
  fixed_effects$Bulk_ESS <- format_number(fixed_effects$Bulk_ESS)
  fixed_effects$Tail_ESS <- format_number(fixed_effects$Tail_ESS)
  Rhat <- format_number(fixed_effects$Rhat, 3)
  
  # Add significance stars
  is_significant <- function(low, high) {
    (low > 0 & high > 0) | (low < 0 & high < 0)
  }
  significance <- ifelse(
    is_significant(fixed_effects$`l-95% CI`, fixed_effects$`u-95% CI`),
    '*', 
    ''
  )
  
  # Handle exponentiation
  if (exponentiate) {
    fixed_effects$Estimate <- exp(fixed_effects$Estimate)
    fixed_effects$`u-95% CI` <- exp(fixed_effects$`u-95% CI`)
    fixed_effects$`l-95% CI` <- exp(fixed_effects$`l-95% CI`)
  }
  
  # Format CIs
  fixed_effects$`l-95% CI` <- format_number(fixed_effects$`l-95% CI`)
  fixed_effects$`u-95% CI` <- format_number(fixed_effects$`u-95% CI`)
  
  # Format estimates with significance
  fixed_effects$Estimate <- ifelse(
    is.na(fixed_effects$Estimate),
    NA,
    paste0(format_number(fixed_effects$Estimate), significance)
  )
  
  # Determine correct column name
  correct_name <- if (exponentiate) {
    if (model$family[[1]] %in% c('bernoulli', 'cumulative')) {
      'OR'
    } else if (model$family[[1]] == 'negbinomial') {
      'IRR'
    } else {
      warning("Coefficients were exponentiated. Double check if this was intended.")
      'exp(Est.)'
    }
  } else {
    'b'
  }
  
  # Rename columns
  colnames(fixed_effects)[1] <- colnames(random_effects)[1] <- correct_name
  
  # Handle row names
  if (!is.null(model_rownames_fixed)) rownames(fixed_effects) <- model_rownames_fixed
  if (!is.null(model_rownames_random)) rownames(random_effects) <- model_rownames_random
  
  # Combine fixed and random effects
  fixed_effects$Rhat <- Rhat
  full_results <- rbind(fixed_effects, random_effects)
  
  if (!short_version) {
    return(full_results[, c(1, 3:10)])
  }
  
  # Create short version with CI
  full_results$`95% CI` <- ifelse(
    is.na(full_results[[correct_name]]) | 
      is.na(full_results$`u-95% CI`) | 
      grepl("NA", full_results$`u-95% CI`) | 
      grepl("^\\s*NA\\s*$", full_results[[correct_name]]),
    NA,
    paste0("[", full_results$`l-95% CI`, ", ", full_results$`u-95% CI`, "]")
  )
  
  return(full_results[, c(1, 8)])
}



#### Check Models

pp_check_transformed <- function(model, transform = log1p) {
  # only for univariate models
  outcome_variable <- strsplit(as.character(formula(model)), " ~ ")[[1]][1]
  pred <- brms::posterior_predict(model)
  
  # Use get to dynamically reference the outcome variable in the model's data
  bayesplot::ppc_dens_overlay(
    y = transform(model$data[[outcome_variable]]),
    yrep = transform(pred[1:10, ])
  )
}


# Check brms models for convergence etc. 
check_brms <- function(
    model, 
    log_pp_check = FALSE, # a function needs to be passed!
    transform = log1p
) { 
  rstan::check_hmc_diagnostics(model$fit)
  plot(model, ask = FALSE, nvariables = 3)
  plot(pp_check(model, type = 'ecdf_overlay'))
  plot(pp_check(model))
  
  if (log_pp_check) {
    plot(pp_check_transformed(model, transform = transform))
  }
  loo(model)
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
    seed = 123
    )
  
  if (isTRUE(plot)) {
    plot(dharma.obj, ...)
  }
  invisible(dharma.obj)
}



DHARMa.check_brms_ordinal <- function(model, plot = TRUE, debug = FALSE, ...) {
  
  # Extract data and posterior predictions
  mdata <- brms::standata(model)
  
  if(debug) {
    cat("Model family:", family(model)$family, "\n")
    cat("Response variable levels:", levels(factor(mdata$Y)), "\n")
  }
  
  # Get posterior predictions
  latent_predictions <- try({
    brms::posterior_epred(model, ndraws = 1000)
  })
  
  if(inherits(latent_predictions, "try-error")) {
    stop("Error in getting posterior predictions. Check if model is ordinal.")
  }
  
  if(debug) {
    cat("Dimensions of latent_predictions:", dim(latent_predictions), "\n")
    cat("Class of latent_predictions:", class(latent_predictions), "\n")
  }
  
  # Convert categorical responses to numeric
  observed_numeric <- as.numeric(factor(mdata$Y))
  
  if(debug) {
    cat("Range of observed_numeric:", range(observed_numeric), "\n")
  }
  
  # For ordinal models, we need to convert the 3D array to a matrix
  # Each row will be a simulation, each column an observation
  sim_response <- try({
    # Convert predictions to probabilities for each category
    probs <- aperm(latent_predictions, c(1, 2, 3))
    # Sample categories based on these probabilities
    samples <- apply(probs, 1:2, function(x) sample(seq_along(x), size = 1, prob = x))
    # Transpose to match DHARMa expectations
    t(samples)
  })
  
  if(inherits(sim_response, "try-error")) {
    stop("Error in processing predictions. Check dimensions.")
  }
  
  if(debug) {
    cat("Dimensions of simulated response:", dim(sim_response), "\n")
  }
  
  # Create fitted predicted response (expected category)
  fitted_pred <- try({
    # Calculate expected category for each observation
    probs_mean <- apply(latent_predictions, 2:3, mean)
    apply(probs_mean, 1, function(x) sum(seq_along(x) * x))
  })
  
  if(inherits(fitted_pred, "try-error")) {
    stop("Error in creating fitted predictions. Check dimensions.")
  }
  
  # Create DHARMa object
  dharma.obj <- try({
    DHARMa::createDHARMa(
      simulatedResponse = sim_response,
      observedResponse = observed_numeric,
      fittedPredictedResponse = fitted_pred,
      integerResponse = TRUE,  # Changed to TRUE for ordinal responses
      seed = 123
    )
  })
  
  if(inherits(dharma.obj, "try-error")) {
    stop("Error in creating DHARMa object. Check all inputs.")
  }
  
  if (isTRUE(plot)) {
    plot(dharma.obj, ...)
  }
  
  invisible(dharma.obj)
}




DHARMa.check_brms.all <- function(model, integer = FALSE, outliers_type = 'default', ...) {
  
  if ("factor" %in% class(model$data[[1]]) 
      & 'ordered' %in% class(model$data[[1]])) {
    model.check <- DHARMa.check_brms_ordinal(model, plot = FALSE, ...)
  } else {
    model.check <- DHARMa.check_brms(model, integer = integer, plot = FALSE, ...)
  }
  
  plot(model.check)
  try(testDispersion(model.check))
  try(testZeroInflation(model.check))
  try(testOutliers(model.check, type = outliers_type))
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
    if (model$family[[1]] %in% c(
      'bernoulli', 
      'negbinomial',
      'cumulative', 
      'hurdle_lognormal',
      'hurdle_poisson',
      'hurdle_negbinomial',
      'hurdle_cumulative',
      'lognormal',
      'skewnormal') | grepl('log', model_name)) {
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

