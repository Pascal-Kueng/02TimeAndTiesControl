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

###########################################################################
####################### REPORT MODELS #####################################
###########################################################################


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





############################################################################
############################ PLOT EFFECTS ##################################
############################################################################

conditional_spaghetti <- function(
    model,                    # a brmsfit object
    effects,                  # character vector of variable names
    group_var = NULL,         # character of variable name
    n_groups = NULL,          # if NULL all slopes are plotted, else n random samples are plotted
    seed = 45,                # seed for random samples
    robust = TRUE,            # TRUE takes the median, FALSE the mean from the samples of the posterior
    plot_full_range = FALSE,  # Plot over full predictor range if TRUE instead of range with data
    x_limits = NULL,          # vector with lower and upper bound of x-axis (optional)
    y_limits = NULL,          # vector with lower and upper bound of y-axis (optional)
    x_label = NULL,           # character
    y_label = NULL,           # character
    transform_fn = NULL       # function to transform values after inverse link
) {
  if (!inherits(model, 'brmsfit')) {
    stop("Only brmsfit objects supported")
  }
  
  # Check for required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required but not installed.")
  }
  
  # Set up labels
  if (is.null(x_label)) {
    x_label <- effects
  }
  
  # Determine default y_label based on family and link
  if (is.null(y_label)) {
    if ((model$family$family %in% c("binomial", "bernoulli")) && model$family$link == "logit") {
      y_label <- 'Probability'
    } else if (model$family$family == "gaussian" && model$family$link == "identity") {
      y_label <- 'Outcome'
    } else if ((model$family$family %in% c("poisson", "negbinomial")) && model$family$link == "log") {
      y_label <- 'Count'
    } else {
      y_label <- 'Response'
    }
  }
  
  # Define inverse link functions for different links
  reverse_link <- function(x) { x } # Default to identity link
  
  # Prepare transformations based on link and family
  if (model$family$link == 'logit') {
    reverse_link <- function(x) {
      exp(x) / (1 + exp(x))
    }
  } else if (model$family$link == 'probit') {
    reverse_link <- function(x) {
      pnorm(x)
    }
  } else if (model$family$link == 'log') {
    reverse_link <- function(x) {
      exp(x)
    }
  } else if (model$family$link == 'identity') {
    reverse_link <- function(x) {
      x
    }
  } else if (model$family$link == 'inverse') {
    reverse_link <- function(x) {
      1 / x
    }
  } else if (model$family$link == 'cloglog') {
    reverse_link <- function(x) {
      1 - exp(-exp(x))
    }
  } else {
    stop("Link function not recognized or not supported.")
  }
  
  # Extract posterior samples for fixed effects
  posterior_samples <- posterior::as_draws_df(model) %>% as.data.frame()
  
  for (i in seq_along(effects)) {
    e <- effects[[i]]
    
    # Define a range of values for the predictor across all data
    predictor_range <- seq(min(model$data[[e]], na.rm = TRUE), max(model$data[[e]], na.rm = TRUE), length.out = 100)
    
    # Extract fixed intercept and slope samples
    fixed_intercept_samples <- posterior_samples$b_Intercept
    fixed_slope_samples <- posterior_samples[[paste0("b_", e)]]
    
    # Handle cases where the effect is not in the fixed effects
    if (is.null(fixed_slope_samples)) {
      stop(paste("Effect", e, "not found in fixed effects."))
    }
    
    # Choose mean or median based on the 'robust' argument
    if (robust) {
      # Use median for robust summaries
      fixed_intercept_central <- median(fixed_intercept_samples)
      fixed_slope_central <- median(fixed_slope_samples)
    } else {
      # Use mean for traditional summaries
      fixed_intercept_central <- mean(fixed_intercept_samples)
      fixed_slope_central <- mean(fixed_slope_samples)
    }
    
    
    # Generate fixed effect predictions across the predictor range
    fixed_predictions <- data.frame(x_value = predictor_range)
    linear_predictor <- fixed_intercept_central + fixed_slope_central * predictor_range
    
    # Apply inverse link and custom transformation if provided
    if (is.null(transform_fn)) {
      fixed_predictions$response <- reverse_link(linear_predictor)
    } else {
      fixed_predictions$response <- transform_fn(reverse_link(linear_predictor))
    }
    
    # Calculate credible intervals for the fixed effect line
    ci_bounds_fixed <- sapply(predictor_range, function(x_val) {
      lp <- fixed_intercept_samples + fixed_slope_samples * x_val
      if (is.null(transform_fn)) {
        reverse_link(lp)
      } else {
        transform_fn(reverse_link(lp))
      }
    })
    fixed_predictions$lower <- apply(ci_bounds_fixed, 2, quantile, probs = 0.025)
    fixed_predictions$upper <- apply(ci_bounds_fixed, 2, quantile, probs = 0.975)
    
    # Initialize ggplot
    p <- ggplot2::ggplot() +
      # Add CI ribbon for the fixed effect
      ggplot2::geom_ribbon(
        data = fixed_predictions,
        ggplot2::aes(x = x_value, ymin = lower, ymax = upper),
        fill = "navy", alpha = 0.1
      ) +
      # Add fixed effect line
      ggplot2::geom_line(
        data = fixed_predictions,
        ggplot2::aes(x = x_value, y = response),
        color = "navy",
        size = 1.3
      ) +
      ggplot2::labs(
        title = paste("Fixed Effects:", x_label[[i]]),
        x = x_label[[i]],
        y = y_label
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      )
    
    # If group_var is provided, add random effects
    if (!is.null(group_var)) {
      # Build the pattern for random effects variable names
      random_effects_prefix <- paste0("r_", group_var)
      random_effects_pattern <- paste0("r_", group_var, "\\[(.*),(.+)]")
      random_effects <- posterior_samples %>%
        dplyr::select(dplyr::starts_with(random_effects_prefix)) %>%
        tidyr::pivot_longer(
          cols = dplyr::everything(),
          names_to = c("group_level", ".value"),
          names_pattern = random_effects_pattern
        ) %>%
        dplyr::mutate(group_level = as.character(group_level)) %>%
        dplyr::select(group_level, Intercept, dplyr::starts_with(e)) # Extract only intercepts and slopes for effect `e`
      
      # Handle cases where random effects for 'e' are not present
      random_slope_name <- e
      if (!random_slope_name %in% names(random_effects)) {
        # If random slope for 'e' is not present, set it to zero
        random_effects[[random_slope_name]] <- 0
      }
      
      # Determine which group levels to include
      unique_groups <- unique(random_effects$group_level)
      if (!is.null(n_groups)) {
        if (n_groups < length(unique_groups)) {
          if (!is.null(seed)) {
            set.seed(seed)  # Set seed for reproducibility
          }
          selected_groups <- sample(unique_groups, n_groups)
        } else {
          selected_groups <- unique_groups
        }
      } else {
        selected_groups <- unique_groups
      }
      
      # Filter random effects to include only selected group levels
      random_effects <- random_effects %>% dplyr::filter(group_level %in% selected_groups)
      
      # Generate individual-level predictions based on the robust setting
      individual_predictions <- do.call(rbind, lapply(selected_groups, function(j) {
        # Filter for the posterior samples of this specific individual
        individual_samples <- random_effects %>% dplyr::filter(group_level == j)
        
        # Compute the combined intercept and slope for this individual
        if (robust) {
          # Use median for robust summaries
          individual_intercept <- fixed_intercept_central + median(individual_samples$Intercept)
          individual_slope <- fixed_slope_central + median(individual_samples[[random_slope_name]])
        } else {
          # Use mean for traditional summaries
          individual_intercept <- fixed_intercept_central + mean(individual_samples$Intercept)
          individual_slope <- fixed_slope_central + mean(individual_samples[[random_slope_name]])
        }
        
        # Determine predictor range for this group
        if (plot_full_range) {
          # Use the full predictor range
          predictor_range_group <- predictor_range
        } else {
          # Get the data for this group
          group_data <- model$data[model$data[[group_var]] == j, ]
          
          # Check if there is data for this group
          if (nrow(group_data) == 0) {
            return(NULL)
          }
          
          # Get the range of predictor 'e' for this group
          min_e <- min(group_data[[e]], na.rm = TRUE)
          max_e <- max(group_data[[e]], na.rm = TRUE)
          
          # Handle case where min and max are the same
          if (min_e == max_e) {
            predictor_range_group <- min_e
          } else {
            predictor_range_group <- seq(min_e, max_e, length.out = 100)
          }
        }
        
        # Generate individual-level predictions
        lp <- individual_intercept + individual_slope * predictor_range_group
        if (is.null(transform_fn)) {
          response <- reverse_link(lp)
        } else {
          response <- transform_fn(reverse_link(lp))
        }
        data.frame(
          group_level = j,
          x_value = predictor_range_group,
          response = response
        )
      }))
      
      # Remove any NULL elements (in case some groups had no data)
      individual_predictions <- individual_predictions %>% dplyr::filter(!is.na(response))
      
      # Add individual lines to the plot
      p <- p +
        ggplot2::geom_line(
          data = individual_predictions,
          ggplot2::aes(x = x_value, y = response, group = group_level),
          color = "grey50",
          size = 0.3,
          alpha = 0.8
        ) +
        ggplot2::labs(
          title = paste("Conditional Fixed and Random Effects:", x_label[[i]])
        )
    }
    
    # Apply x and y limits if specified
    if (!is.null(x_limits)) {
      p <- p + ggplot2::xlim(x_limits)
    }
    if (!is.null(y_limits)) {
      p <- p + ggplot2::ylim(y_limits)
    }
    
    # Plot final result
    print(p)
    gc()
  }
}






####################



# Main function
conditional_spaghetti <- function(
    model,                    # a brmsfit object
    effects,                  # character vector of variable names
    group_var = NULL,         # character of variable name
    n_groups = NULL,          # if NULL all slopes are plotted, else n random samples are plotted
    seed = 45,                # seed for random samples
    robust = TRUE,            # TRUE takes the median, FALSE the mean from the samples of the posterior
    plot_full_range = FALSE,  # Plot over full predictor range if TRUE instead of range with data
    x_limits = NULL,          # vector with lower and upper bound of x-axis (optional)
    y_limits = NULL,          # vector with lower and upper bound of y-axis (optional)
    x_label = NULL,           # character
    y_label = NULL,           # character
    transform_fn = NULL       # function to transform values after inverse link
) {
  if (!inherits(model, 'brmsfit')) {
    stop("Only brmsfit objects supported")
  }
  
  # Check for required packages
  required_packages <- c("dplyr", "tidyr", "tidyverse", "ggplot2", "posterior")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Determine if the model is a hurdle or zero-inflated model
  is_hurdle_model <- grepl("^hurdle_", model$family$family)
  is_zi_model <- grepl("^zero_inflated_", model$family$family)
  
  # Call appropriate plotting function based on model type
  if (is_hurdle_model || is_zi_model) {
    plots <- plot_hurdle_model(
      model = model,
      effects = effects,
      group_var = group_var,
      n_groups = n_groups,
      seed = seed,
      robust = robust,
      plot_full_range = plot_full_range,
      x_limits = x_limits,
      y_limits = y_limits,
      x_label = x_label,
      y_label = y_label,
      transform_fn = transform_fn
    )
  } else {
    plots <- plot_general_model(
      model = model,
      effects = effects,
      group_var = group_var,
      n_groups = n_groups,
      seed = seed,
      robust = robust,
      plot_full_range = plot_full_range,
      x_limits = x_limits,
      y_limits = y_limits,
      x_label = x_label,
      y_label = y_label,
      transform_fn = transform_fn
    )
  }
  
  # Return the plots
  return(plots)
}



# Function for general models
plot_general_model <- function(
    model,
    effects,
    group_var,
    n_groups,
    seed,
    robust,
    plot_full_range,
    x_limits,
    y_limits,
    x_label,
    y_label,
    transform_fn
) {
  plots_list <- list()
  
  # Extract posterior samples for fixed effects
  posterior_samples <- posterior::as_draws_df(model) %>% as.data.frame()
  
  for (i in seq_along(effects)) {
    e <- effects[[i]]
    
    # Define a range of values for the predictor across all data
    predictor_range <- seq(min(model$data[[e]], na.rm = TRUE), max(model$data[[e]], na.rm = TRUE), length.out = 100)
    
    # Set up labels
    x_lab <- ifelse(is.null(x_label), e, x_label[[i]])
    y_lab <- ifelse(is.null(y_label), 'Response', y_label)
    
    # Define inverse link function
    inverse_link_functions <- list(
      logit = plogis,
      probit = pnorm,
      cloglog = function(x) 1 - exp(-exp(x)),
      identity = function(x) x,
      log = exp,
      inverse = function(x) 1 / x
    )
    reverse_link <- inverse_link_functions[[model$family$link]]
    if (is.null(reverse_link)) {
      stop("Link function not recognized or not supported.")
    }
    
    # Extract fixed intercept and slope samples
    fixed_intercept_samples <- posterior_samples$b_Intercept
    fixed_slope_samples <- posterior_samples[[paste0("b_", e)]]
    
    # Handle cases where the effect is not in the fixed effects
    if (is.null(fixed_slope_samples)) {
      stop(paste("Effect", e, "not found in fixed effects."))
    }
    
    # Choose mean or median based on the 'robust' argument
    if (robust) {
      # Use median for robust summaries
      fixed_intercept_central <- median(fixed_intercept_samples)
      fixed_slope_central <- median(fixed_slope_samples)
    } else {
      # Use mean for traditional summaries
      fixed_intercept_central <- mean(fixed_intercept_samples)
      fixed_slope_central <- mean(fixed_slope_samples)
    }
    
    # Generate fixed effect predictions across the predictor range
    linear_predictor <- fixed_intercept_central + fixed_slope_central * predictor_range
    
    # Apply inverse link and custom transformation if provided
    response <- reverse_link(linear_predictor)
    if (!is.null(transform_fn)) {
      response <- transform_fn(response)
    }
    fixed_predictions <- data.frame(
      x_value = predictor_range,
      response = response
    )
    
    # Calculate credible intervals
    ci_bounds_fixed <- sapply(predictor_range, function(x_val) {
      lp <- fixed_intercept_samples + fixed_slope_samples * x_val
      response_samples <- reverse_link(lp)
      if (!is.null(transform_fn)) {
        response_samples <- transform_fn(response_samples)
      }
      response_samples
    })
    fixed_predictions$lower <- apply(ci_bounds_fixed, 2, quantile, probs = 0.025)
    fixed_predictions$upper <- apply(ci_bounds_fixed, 2, quantile, probs = 0.975)
    
    # Initialize ggplot
    p <- ggplot() +
      # Add CI ribbon for the fixed effect
      geom_ribbon(
        data = fixed_predictions,
        aes(x = x_value, ymin = lower, ymax = upper),
        fill = "navy", alpha = 0.1
      ) +
      # Add fixed effect line
      geom_line(
        data = fixed_predictions,
        aes(x = x_value, y = response),
        color = "navy",
        size = 1.3
      ) +
      labs(
        title = paste("Fixed Effects:", x_lab),
        x = x_lab,
        y = y_lab
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    # If group_var is provided, add random effects
    if (!is.null(group_var)) {
      # Handle random effects as before
      random_effects_prefix <- paste0("r_", group_var)
      random_effects_pattern <- paste0("r_", group_var, "\\[(.*),(.+)]")
      random_effects <- posterior_samples %>%
        select(starts_with(random_effects_prefix)) %>%
        pivot_longer(
          cols = everything(),
          names_to = c("group_level", ".value"),
          names_pattern = random_effects_pattern
        ) %>%
        mutate(group_level = as.character(group_level)) %>%
        select(group_level, Intercept, starts_with(e))
      
      # Handle cases where random effects for 'e' are not present
      random_slope_name <- e
      if (!random_slope_name %in% names(random_effects)) {
        # If random slope for 'e' is not present, set it to zero
        random_effects[[random_slope_name]] <- 0
      }
      
      # Determine which group levels to include
      unique_groups <- unique(random_effects$group_level)
      if (!is.null(n_groups)) {
        if (n_groups < length(unique_groups)) {
          if (!is.null(seed)) {
            set.seed(seed)  # Set seed for reproducibility
          }
          selected_groups <- sample(unique_groups, n_groups)
        } else {
          selected_groups <- unique_groups
        }
      } else {
        selected_groups <- unique_groups
      }
      
      # Filter random effects to include only selected group levels
      random_effects <- random_effects %>% filter(group_level %in% selected_groups)
      
      # Generate individual-level predictions
      individual_predictions <- do.call(rbind, lapply(selected_groups, function(j) {
        # Filter for the posterior samples of this specific group
        individual_samples <- random_effects %>% filter(group_level == j)
        
        # Compute the combined intercept and slope for this group
        if (robust) {
          # Use median for robust summaries
          individual_intercept <- fixed_intercept_central + median(individual_samples$Intercept, na.rm = TRUE)
          individual_slope <- fixed_slope_central + median(individual_samples[[random_slope_name]], na.rm = TRUE)
        } else {
          # Use mean for traditional summaries
          individual_intercept <- fixed_intercept_central + mean(individual_samples$Intercept, na.rm = TRUE)
          individual_slope <- fixed_slope_central + mean(individual_samples[[random_slope_name]], na.rm = TRUE)
        }
        
        # Determine predictor range for this group
        if (plot_full_range) {
          # Use the full predictor range
          predictor_range_group <- predictor_range
        } else {
          # Get the data for this group
          group_data <- model$data[model$data[[group_var]] == j, ]
          
          # Check if there is data for this group
          if (nrow(group_data) == 0) {
            return(NULL)
          }
          
          # Get the range of predictor 'e' for this group
          min_e <- min(group_data[[e]], na.rm = TRUE)
          max_e <- max(group_data[[e]], na.rm = TRUE)
          
          # Handle case where min and max are the same
          if (min_e == max_e) {
            predictor_range_group <- min_e
          } else {
            predictor_range_group <- seq(min_e, max_e, length.out = 100)
          }
        }
        
        # Generate individual-level predictions
        lp <- individual_intercept + individual_slope * predictor_range_group
        response <- reverse_link(lp)
        if (!is.null(transform_fn)) {
          response <- transform_fn(response)
        }
        data.frame(
          group_level = j,
          x_value = predictor_range_group,
          response = response
        )
      }))
      
      # Remove any NULL elements (in case some groups had no data)
      individual_predictions <- individual_predictions %>% filter(!is.na(response))
      
      # Add individual lines to the plot
      p <- p +
        geom_line(
          data = individual_predictions,
          aes(x = x_value, y = response, group = group_level),
          color = "grey50",
          size = 0.3,
          alpha = 0.8
        ) +
        labs(
          title = paste("Conditional Fixed and Random Effects:", x_lab)
        )
    }
    
    # Apply x and y limits if specified
    if (!is.null(x_limits)) {
      p <- p + xlim(x_limits)
    }
    if (!is.null(y_limits)) {
      p <- p + ylim(y_limits)
    }
    
    # Store the plot
    plots_list[[e]] <- p
  }
  
  return(plots_list)
}




# Function for hurdle and zero-inflated models
# Function for hurdle and zero-inflated models (corrected)
plot_hurdle_model <- function(
    model,
    effects,
    group_var,
    n_groups,
    seed,
    robust,
    plot_full_range,
    x_limits,
    y_limits,
    x_label,
    y_label,
    transform_fn
) {
  plots_list <- list()
  
  # Extract posterior samples for fixed effects
  posterior_samples <- posterior::as_draws_df(model) %>% as.data.frame()
  
  for (i in seq_along(effects)) {
    e <- effects[[i]]
    
    # Define a range of values for the predictor across all data
    predictor_range <- seq(min(model$data[[e]], na.rm = TRUE), max(model$data[[e]], na.rm = TRUE), length.out = 100)
    
    # Set up labels
    x_lab <- ifelse(is.null(x_label), e, x_label[[i]])
    y_lab <- ifelse(is.null(y_label), 'Response', y_label)
    
    # Identify the hurdle or zero-inflation component
    if (grepl("^hurdle_", model$family$family)) {
      hurdle_component <- "hu"
      model_type <- "hurdle"
    } else {
      hurdle_component <- "zi"
      model_type <- "zero_inflated"
    }
    
    # Get link functions for both components
    count_link <- model$family$link
    if (!is.null(model$formula$auxiliary[[hurdle_component]]$link)) {
      hurdle_link <- model$formula$auxiliary[[hurdle_component]]$link
    } else {
      # Default link for hurdle/zi component is 'logit' in brms
      hurdle_link <- "logit"
    }
    
    # Define inverse link functions
    inverse_link_functions <- list(
      logit = plogis,
      probit = pnorm,
      cloglog = function(x) 1 - exp(-exp(x)),
      identity = function(x) x,
      log = exp,
      inverse = function(x) 1 / x
    )
    reverse_link_count <- inverse_link_functions[[count_link]]
    reverse_link_hurdle <- inverse_link_functions[[hurdle_link]]
    if (is.null(reverse_link_count) || is.null(reverse_link_hurdle)) {
      stop("Unsupported link function in one of the model components.")
    }
    
    # Extract fixed intercept and slope samples for count component
    fixed_intercept_samples_count <- posterior_samples$b_Intercept
    fixed_slope_samples_count <- posterior_samples[[paste0("b_", e)]]
    
    # Extract fixed intercept and slope samples for hurdle component
    fixed_intercept_samples_hurdle <- posterior_samples[[paste0("b_", hurdle_component, "_Intercept")]]
    fixed_slope_samples_hurdle <- posterior_samples[[paste0("b_", hurdle_component, "_", e)]]
    
    # Handle cases where the effect is not in the fixed effects
    if (is.null(fixed_slope_samples_count)) {
      stop(paste("Effect", e, "not found in fixed effects for count component."))
    }
    if (is.null(fixed_slope_samples_hurdle)) {
      # If the effect is not in the hurdle component, set slope to zero
      fixed_slope_samples_hurdle <- rep(0, length(fixed_intercept_samples_hurdle))
    }
    
    # Choose mean or median based on the 'robust' argument
    if (robust) {
      # Use median for robust summaries
      fixed_intercept_central_count <- median(fixed_intercept_samples_count)
      fixed_slope_central_count <- median(fixed_slope_samples_count)
      fixed_intercept_central_hurdle <- median(fixed_intercept_samples_hurdle)
      fixed_slope_central_hurdle <- median(fixed_slope_samples_hurdle)
    } else {
      # Use mean for traditional summaries
      fixed_intercept_central_count <- mean(fixed_intercept_samples_count)
      fixed_slope_central_count <- mean(fixed_slope_samples_count)
      fixed_intercept_central_hurdle <- mean(fixed_intercept_samples_hurdle)
      fixed_slope_central_hurdle <- mean(fixed_slope_samples_hurdle)
    }
    
    # Generate linear predictors
    linear_predictor_count <- fixed_intercept_central_count + fixed_slope_central_count * predictor_range
    linear_predictor_hurdle <- fixed_intercept_central_hurdle + fixed_slope_central_hurdle * predictor_range
    
    # Apply inverse link functions
    mu_count <- reverse_link_count(linear_predictor_count)
    prob_zero <- reverse_link_hurdle(linear_predictor_hurdle)  # P(Y = 0)
    
    # For expected value, we need (1 - P(Y = 0))
    prob_positive <- 1 - prob_zero  # P(Y > 0)
    
    # Expected value
    expected_value <- prob_positive * mu_count
    
    # Apply custom transformation if provided
    if (!is.null(transform_fn)) {
      mu_count <- transform_fn(mu_count)
      prob_zero <- transform_fn(prob_zero)
      prob_positive <- transform_fn(prob_positive)
      expected_value <- transform_fn(expected_value)
    }
    
    # Create data frames for each component
    fixed_predictions_count <- data.frame(
      x_value = predictor_range,
      response = mu_count
    )
    fixed_predictions_hurdle <- data.frame(
      x_value = predictor_range,
      response = prob_positive
    )
    fixed_predictions_combined <- data.frame(
      x_value = predictor_range,
      response = expected_value
    )
    
    # Calculate credible intervals for each component
    ci_bounds_count <- sapply(predictor_range, function(x_val) {
      lp_count <- fixed_intercept_samples_count + fixed_slope_samples_count * x_val
      mu_count_samples <- reverse_link_count(lp_count)
      if (!is.null(transform_fn)) {
        mu_count_samples <- transform_fn(mu_count_samples)
      }
      mu_count_samples
    })
    fixed_predictions_count$lower <- apply(ci_bounds_count, 2, quantile, probs = 0.025)
    fixed_predictions_count$upper <- apply(ci_bounds_count, 2, quantile, probs = 0.975)
    
    ci_bounds_hurdle <- sapply(predictor_range, function(x_val) {
      lp_hurdle <- fixed_intercept_samples_hurdle + fixed_slope_samples_hurdle * x_val
      prob_zero_samples <- reverse_link_hurdle(lp_hurdle)
      prob_positive_samples <- 1 - prob_zero_samples
      if (!is.null(transform_fn)) {
        prob_positive_samples <- transform_fn(prob_positive_samples)
      }
      prob_positive_samples
    })
    fixed_predictions_hurdle$lower <- apply(ci_bounds_hurdle, 2, quantile, probs = 0.025)
    fixed_predictions_hurdle$upper <- apply(ci_bounds_hurdle, 2, quantile, probs = 0.975)
    
    ci_bounds_combined <- sapply(predictor_range, function(x_val) {
      lp_count <- fixed_intercept_samples_count + fixed_slope_samples_count * x_val
      lp_hurdle <- fixed_intercept_samples_hurdle + fixed_slope_samples_hurdle * x_val
      mu_count_samples <- reverse_link_count(lp_count)
      prob_zero_samples <- reverse_link_hurdle(lp_hurdle)
      prob_positive_samples <- 1 - prob_zero_samples
      expected_value_samples <- prob_positive_samples * mu_count_samples
      if (!is.null(transform_fn)) {
        expected_value_samples <- transform_fn(expected_value_samples)
      }
      expected_value_samples
    })
    fixed_predictions_combined$lower <- apply(ci_bounds_combined, 2, quantile, probs = 0.025)
    fixed_predictions_combined$upper <- apply(ci_bounds_combined, 2, quantile, probs = 0.975)
    
    # Create plots for each component
    p_count <- ggplot() +
      geom_ribbon(
        data = fixed_predictions_count,
        aes(x = x_value, ymin = lower, ymax = upper),
        fill = "blue", alpha = 0.1
      ) +
      geom_line(
        data = fixed_predictions_count,
        aes(x = x_value, y = response),
        color = "blue",
        size = 1.3
      ) +
      labs(
        title = paste("Count Component:", x_lab),
        x = x_lab,
        y = "Expected Count Given Positive Outcome"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    p_hurdle <- ggplot() +
      geom_ribbon(
        data = fixed_predictions_hurdle,
        aes(x = x_value, ymin = lower, ymax = upper),
        fill = "green", alpha = 0.1
      ) +
      geom_line(
        data = fixed_predictions_hurdle,
        aes(x = x_value, y = response),
        color = "green",
        size = 1.3
      ) +
      labs(
        title = paste(ifelse(model_type == "hurdle", "Probability of Positive Outcome", "Probability of Not Zero-Inflated"), ":", x_lab),
        x = x_lab,
        y = "Probability"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    p_combined <- ggplot() +
      geom_ribbon(
        data = fixed_predictions_combined,
        aes(x = x_value, ymin = lower, ymax = upper),
        fill = "purple", alpha = 0.1
      ) +
      geom_line(
        data = fixed_predictions_combined,
        aes(x = x_value, y = response),
        color = "purple",
        size = 1.3
      ) +
      labs(
        title = paste("Combined Expected Value:", x_lab),
        x = x_lab,
        y = y_lab
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    # Store the plots
    plots_list[[paste0(e, "_count")]] <- p_count
    plots_list[[paste0(e, "_hurdle")]] <- p_hurdle
    plots_list[[paste0(e, "_combined")]] <- p_combined
  }
  
  return(plots_list)
}




##############################################################################
############################ Check Models ####################################
##############################################################################
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

