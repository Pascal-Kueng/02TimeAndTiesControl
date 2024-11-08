# Function for hurdle and zero-inflated models with adjustable layout
plot_hurdle_model <- function(
    model,
    effects,
    group_var = NULL,
    n_groups = NULL,
    seed = NULL,
    robust = FALSE,
    plot_full_range = TRUE,
    x_limits = NULL,
    y_limits = NULL,
    x_label = NULL,
    y_label = NULL,
    transform_fn = NULL,
    use_pr_notation = FALSE,  # Option to switch between 'P' and 'Pr'
    filter_quantiles = NULL,
    font_family = 'Segoe UI'
) {
  # Load required packages
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(posterior)
  library(patchwork)
  library(grid)  
  library(ggridges)
  library(extrafont)
  
  plots_list <- list()
  

  # Extract posterior samples for fixed effects
  posterior_samples <- posterior::as_draws_df(model) %>% as.data.frame()
  
  # Extract sigma samples if necessary
  if (grepl("lognormal", model$family$family)) {
    sigma_samples <- posterior_samples$`sigma`
  }
  
  for (i in seq_along(effects)) {
    e <- effects[[i]]
    
    # Define a range of values for the predictor across all data
    predictor_range <- seq(
      min(model$data[[e]], na.rm = TRUE),
      max(model$data[[e]], na.rm = TRUE),
      length.out = 100
    )
    
    
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
      if (grepl("lognormal", model$family$family)) {
        sigma_central <- median(sigma_samples)
      }
    } else {
      # Use mean for traditional summaries
      fixed_intercept_central_count <- mean(fixed_intercept_samples_count)
      fixed_slope_central_count <- mean(fixed_slope_samples_count)
      fixed_intercept_central_hurdle <- mean(fixed_intercept_samples_hurdle)
      fixed_slope_central_hurdle <- mean(fixed_slope_samples_hurdle)
      if (grepl("lognormal", model$family$family)) {
        sigma_central <- mean(sigma_samples)
      }
    }
    
    # Generate linear predictors
    linear_predictor_count <- fixed_intercept_central_count + fixed_slope_central_count * predictor_range
    linear_predictor_hurdle <- fixed_intercept_central_hurdle + fixed_slope_central_hurdle * predictor_range
    
    # Apply inverse link functions
    if (grepl("lognormal", model$family$family)) {
      mu_count <- exp(linear_predictor_count)
      mu_count_expected <- exp(linear_predictor_count + 0.5 * sigma_central^2)
    } else {
      mu_count <- reverse_link_count(linear_predictor_count)
      mu_count_expected <- mu_count
    }
    prob_zero <- reverse_link_hurdle(linear_predictor_hurdle)  # P(Y = 0)
    prob_positive <- 1 - prob_zero  # P(Y > 0)
    
    # Expected value
    expected_value <- prob_positive * mu_count_expected
    
    # Apply custom transformation if provided
    if (!is.null(transform_fn)) {
      mu_count <- transform_fn(mu_count)
      mu_count_expected <- transform_fn(mu_count_expected)
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
      lp <- fixed_intercept_samples_count + fixed_slope_samples_count * x_val
      if (grepl("lognormal", model$family$family)) {
        response_samples <- exp(lp)
      } else {
        response_samples <- reverse_link_count(lp)
      }
      if (!is.null(transform_fn)) {
        response_samples <- transform_fn(response_samples)
      }
      response_samples
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
      if (grepl("lognormal", model$family$family)) {
        mu_count_samples <- exp(lp_count + 0.5 * sigma_samples^2)
      } else {
        mu_count_samples <- reverse_link_count(lp_count)
      }
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
    
    
    # Set up labels
    ## Plot titles
    if (grepl("lognormal", model$family$family)) {
      positive_component_title <- "Non-Zero Component"
    } else {
      positive_component_title <- "Count Component"
    }
    
    ## Axes Titles
    x_lab <- ifelse(is.null(x_label), e, x_label[[i]])
    
    outcome_names <- NULL
    single_outcome_name <- NULL
    if (!is.null(y_label) & length(y_label) == 3) {
      outcome_names <- y_label
    } else {
      single_outcome_name <- ifelse(is.null(y_label), "Y", y_label)
    }
    
    if (!is.null(single_outcome_name)) {
      positive_component_ylabel <- bquote(E*""[.(single_outcome_name) ~ "|" ~ .(single_outcome_name) ~ ">" ~ 0])
      # Use 'P' or 'Pr' based on the parameter
      prob_label <- if (use_pr_notation) {
        bquote(Pr(.(single_outcome_name) ~ ">" ~ 0))
      } else {
        bquote(P(.(single_outcome_name) ~ ">" ~ 0))
      }
      combined_ylabel <- bquote(E*""[.(single_outcome_name)])
    } else {
      positive_component_ylabel <- outcome_names[[2]]
      prob_label <- outcome_names[[1]]
      combined_ylabel <- outcome_names[[3]]
    }
    # Create plots for each component
    p_hurdle <- ggplot() +
      geom_ribbon(
        data = fixed_predictions_hurdle,
        aes(x = x_value, ymin = lower, ymax = upper),
        fill = "#33a02c", alpha = 0.24
      ) +
      geom_line(
        data = fixed_predictions_hurdle,
        aes(x = x_value, y = response),
        color = "#33a02c",
        linewidth =1.8
      ) +
      labs(
        title = "Hurdle Component",
        x = x_lab,
        y = prob_label
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(t = 10, b = 10)),
        axis.title.x = element_text(size = 12, margin = margin(t = 10)),
        axis.title.y = element_text(size = 12, margin = margin(r = 10, l = 10)),
        axis.text = element_text(size = 10.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        panel.border = element_blank(),         # Removes the main plot panel border
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(color = "grey80", linetype = "dotted"), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = 'grey45')
      )
    
    p_count <- ggplot() +
      geom_ribbon(
        data = fixed_predictions_count,
        aes(x = x_value, ymin = lower, ymax = upper),
        fill = "#1f78b4", alpha = 0.24
      ) +
      geom_line(
        data = fixed_predictions_count,
        aes(x = x_value, y = response),
        color = "#1f78b4",
        linewidth =1.8
      ) +
      labs(
        title = positive_component_title,
        x = x_lab,
        y = positive_component_ylabel
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(t = 10, b = 10)),
        axis.title.x = element_text(size = 12, margin = margin(t = 10)),
        axis.title.y = element_text(size = 12, margin = margin(r = 10, l = 10)),
        axis.text = element_text(size = 10.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        panel.border = element_blank(),         # Removes the main plot panel border
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(color = "grey80", linetype = "dotted"), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = 'grey45')
      )
    
    p_combined <- ggplot() +
      geom_ribbon(
        data = fixed_predictions_combined,
        aes(x = x_value, ymin = lower, ymax = upper),
        fill = "#6a3d9a", alpha = 0.24
      ) +
      geom_line(
        data = fixed_predictions_combined,
        aes(x = x_value, y = response),
        color = "#6a3d9a",
        linewidth =1.8
      ) +
      labs(
        title = "Combined Expected Value",
        x = x_lab,
        y = combined_ylabel
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(t = 10, b = 10)),
        axis.title.x = element_text(size = 12, margin = margin(t = 10)),
        axis.title.y = element_text(size = 12, margin = margin(r = 10, l = 10)),
        axis.text = element_text(size = 10.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        panel.border = element_blank(),         # Removes the main plot panel border
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(color = "grey80", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = 'grey45')
      )
    
    # Include the random effects code
    if (!is.null(group_var)) {
      # Extract random effects for count component
      random_effects_prefix_count <- paste0("r_", group_var)
      random_effects_pattern_count <- paste0("r_", group_var, "\\[(.*),(.+)]")
      random_effects_count <- posterior_samples %>%
        select(.draw, starts_with(random_effects_prefix_count)) %>%
        select(-matches(paste0("__", hurdle_component))) %>%  # Exclude hurdle component random effects
        pivot_longer(
          cols = -c(.draw),
          names_to = c("group_level", ".value"),
          names_pattern = random_effects_pattern_count
        ) %>%
        mutate(group_level = as.character(group_level))
      
      # Random effects for hurdle component
      random_effects_prefix_hurdle <- paste0("r_", group_var, "__", hurdle_component)
      random_effects_pattern_hurdle <- paste0("r_", group_var, "__", hurdle_component, "\\[(.*),(.+)]")
      random_effects_hurdle <- posterior_samples %>%
        select(.draw, starts_with(random_effects_prefix_hurdle)) %>%
        pivot_longer(
          cols = -c(.draw),
          names_to = c("group_level", ".value"),
          names_pattern = random_effects_pattern_hurdle
        ) %>%
        mutate(group_level = as.character(group_level))
      
      # Adjust column names
      names(random_effects_count)[-(1:2)] <- paste0(names(random_effects_count)[-(1:2)], "_count")
      names(random_effects_hurdle)[-(1:2)] <- paste0(names(random_effects_hurdle)[-(1:2)], "_hurdle")
      
      # Merge random effects
      random_effects <- full_join(random_effects_count, random_effects_hurdle, by = c("group_level", ".draw"))
      
      # Handle missing random effects
      random_slope_name <- e
      if (!paste0(random_slope_name, "_count") %in% names(random_effects)) {
        random_effects[[paste0(random_slope_name, "_count")]] <- 0
      }
      if (!paste0(random_slope_name, "_hurdle") %in% names(random_effects)) {
        random_effects[[paste0(random_slope_name, "_hurdle")]] <- 0
      }
      if (!"Intercept_count" %in% names(random_effects)) {
        random_effects$Intercept_count <- 0
      }
      if (!"Intercept_hurdle" %in% names(random_effects)) {
        random_effects$Intercept_hurdle <- 0
      }
      
      # Select groups to include
      unique_groups <- unique(random_effects$group_level)
      if (!is.null(n_groups)) {
        if (n_groups < length(unique_groups)) {
          if (!is.null(seed)) {
            set.seed(seed)  # For reproducibility
          }
          selected_groups <- sample(unique_groups, n_groups)
        } else {
          selected_groups <- unique_groups
        }
      } else {
        selected_groups <- unique_groups
      }
      
      # Filter random effects for selected groups
      random_effects <- random_effects %>% filter(group_level %in% selected_groups)
      
      # Generate individual-level predictions
      individual_predictions <- do.call(rbind, lapply(selected_groups, function(j) {
        # Filter for the posterior samples of this specific group
        individual_samples <- random_effects %>% filter(group_level == j)
        
        # Compute the combined intercept and slope for this group
        if (robust) {
          # Use median for robust summaries
          individual_intercept_count <- fixed_intercept_central_count + median(individual_samples$Intercept_count, na.rm = TRUE)
          individual_slope_count <- fixed_slope_central_count + median(individual_samples[[paste0(e, "_count")]], na.rm = TRUE)
          individual_intercept_hurdle <- fixed_intercept_central_hurdle + median(individual_samples$Intercept_hurdle, na.rm = TRUE)
          individual_slope_hurdle <- fixed_slope_central_hurdle + median(individual_samples[[paste0(e, "_hurdle")]], na.rm = TRUE)
        } else {
          # Use mean for traditional summaries
          individual_intercept_count <- fixed_intercept_central_count + mean(individual_samples$Intercept_count, na.rm = TRUE)
          individual_slope_count <- fixed_slope_central_count + mean(individual_samples[[paste0(e, "_count")]], na.rm = TRUE)
          individual_intercept_hurdle <- fixed_intercept_central_hurdle + mean(individual_samples$Intercept_hurdle, na.rm = TRUE)
          individual_slope_hurdle <- fixed_slope_central_hurdle + mean(individual_samples[[paste0(e, "_hurdle")]], na.rm = TRUE)
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
        # For count component
        lp_count <- individual_intercept_count + individual_slope_count * predictor_range_group
        if (grepl("lognormal", model$family$family)) {
          mu_count <- exp(lp_count)
          mu_count_expected <- exp(lp_count + 0.5 * sigma_central^2)
        } else {
          mu_count <- reverse_link_count(lp_count)
          mu_count_expected <- mu_count
        }
        
        # For hurdle component
        lp_hurdle <- individual_intercept_hurdle + individual_slope_hurdle * predictor_range_group
        prob_zero <- reverse_link_hurdle(lp_hurdle)
        prob_positive <- 1 - prob_zero
        
        # Expected value
        expected_value <- prob_positive * mu_count_expected
        
        # Apply custom transformation if provided
        if (!is.null(transform_fn)) {
          mu_count <- transform_fn(mu_count)
          prob_positive <- transform_fn(prob_positive)
          expected_value <- transform_fn(expected_value)
        }
        
        # Return data frames for each component
        data.frame(
          group_level = j,
          x_value = predictor_range_group,
          mu_count = mu_count,
          prob_positive = prob_positive,
          expected_value = expected_value
        )
      }))
      
      # Remove any NULL elements
      individual_predictions <- individual_predictions %>% filter(!is.na(expected_value))
      
      # Add individual lines to the plots
      # For hurdle component
      p_hurdle <- p_hurdle +
        geom_line(
          data = individual_predictions,
          aes(x = x_value, y = prob_positive, group = group_level),
          color = "#33a02c",
          linewidth = 0.25,
          alpha = 0.50
        )
      
      # For positive outcome component
      p_count <- p_count +
        geom_line(
          data = individual_predictions,
          aes(x = x_value, y = mu_count, group = group_level),
          color = "#1f78b4",
          linewidth = 0.25,
          alpha = 0.50
        )
      
      # For combined expected value
      p_combined <- p_combined +
        geom_line(
          data = individual_predictions,
          aes(x = x_value, y = expected_value, group = group_level),
          color = "#6a3d9a",
          linewidth = 0.25,
          alpha = 0.50
        )
    }  # End of random effects code
    
    # Apply x and y limits if specified
    if (!is.null(x_limits)) {
      p_count <- p_count + xlim(x_limits) 
      p_hurdle <- p_hurdle + xlim(x_limits)
      p_combined <- p_combined + xlim(x_limits)
    }
    if (!is.null(y_limits)) {
      p_count <- p_count + ylim(y_limits)
      p_combined <- p_combined + ylim(y_limits)
    }
    
    
    
    ## create posterior density plot
    # Create additional density plot for transformed slopes
    # Prepare the posterior density plot using the current effect
    draws_df <- posterior_samples %>%
      as_draws_df() %>%
      mutate(
        combined_Intercept = (1 - plogis(posterior_samples[[paste0("b_", hurdle_component, "_Intercept")]])) * exp(posterior_samples$b_Intercept),
        predictions = 
          (1 - plogis(posterior_samples[[paste0("b_", hurdle_component, "_Intercept")]] + posterior_samples[[paste0("b_", hurdle_component, "_", e)]])) *
          exp(posterior_samples$b_Intercept + posterior_samples[[paste0("b_", e)]]),
        combined_slope = predictions / combined_Intercept,
        hurdle_slope = 1 / exp(posterior_samples[[paste0("b_", hurdle_component, "_", e)]]),
        nonzero_slope = exp(posterior_samples[[paste0("b_", e)]])
      ) %>%
      select(c(hurdle_slope, nonzero_slope, combined_slope))
    
    # Prepare the data for the current effect
    plot_data <- draws_df %>%
      pivot_longer(cols = c(hurdle_slope, nonzero_slope, combined_slope),
                   names_to = "slope_type",
                   values_to = "value") %>%
      mutate(
        slope_type = factor(
          slope_type, 
          levels = c("combined_slope", "nonzero_slope", "hurdle_slope"),
          labels = c("Combined", "Non-Zero Component", "Hurdle Component")
        )
      )
    
    # Calculate the quantiles for each slope type
    if (!is.null(filter_quantiles)) {
      quantiles <- plot_data %>%
        group_by(slope_type) %>%
        summarize(
          lower = quantile(value, 1 - filter_quantiles),
          upper = quantile(value, filter_quantiles)
        )
      
      # Filter the data to include only values within the 99.99% quantile range
      plot_data <- plot_data %>%
        left_join(quantiles, by = "slope_type") %>%
        filter(value >= lower & value <= upper) %>%
        select(-lower, -upper)
    }
    
    # Create the vertical density plot
    p_density <- ggplot(plot_data, aes(x = value, y = slope_type, height = after_stat(density))) +
      geom_density_ridges_gradient(
        color = 'black',
        linewidth = 0.25,
        aes(fill = after_stat(x > 1)),
        scale = 4,
        rel_min_height = 0.0001,
        gradient_lwd = 0.1,
        quantile_lines = TRUE, 
        quantiles = 0.5  # Adding lines for the median
      ) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
      scale_y_discrete(expand = c(0.01, 0)) +
      scale_x_continuous(expand = c(0.01, 0)) +
      scale_fill_manual(
        values = c("lightcoral", "steelblue2"),
        name = "Effect Direction",
        labels = c("Negative (<1)", "Positive (>1)")
      ) +
      theme_ridges(font_size = 12, grid = TRUE) +  # Reduced font size for a cleaner look
      theme(
        panel.grid.major = element_blank(),  # Remove major grid lines for less clutter
        panel.grid.minor.x = element_line(color = "grey80", linetype = "dotted"),  # Keep only minor grid lines on x-axis
        panel.grid.minor.y = element_blank(),  # Remove y-axis grid lines
        axis.title.x = element_text(hjust = 0.5, size = 12, margin = margin(t = 10)),  # Slightly smaller x-axis title
        axis.text.x = element_text(size = 10.5),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(t = 10, b = 10)),  # Slightly smaller plot title
        plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 10.5),
        legend.position = 'right',  
        legend.title = element_text(face = "bold"),
        legend.background = element_blank(),  # Remove legend background for a cleaner look
        legend.key.size = unit(0.7, "cm")
      ) +
      labs(
        x = "Slopes Representing Multiplicative Changes:\nOdds Ratios for Hurdle Component and Expected Values for other Components",
        y = NULL,
        title = "Posterior Density of Transformed Slopes",
        subtitle = "Median lines shown for each component"
      )
    
    
    
    # Arrange plots for vertical option
    design <- "
      AAACCCCCCC
      AAACCCCCCC
      AAACCCCCCC
      BBBCCCCCCC
      BBBCCCCCCC
      BBBCCCCCCC
      DDDDDDDDDD
      DDDDDDDDDD
      DDDDDDDDDD
      DDDDDDDDDD
    "
    combined_plot <- 
      p_hurdle + p_count + p_combined + free(p_density) + 
      plot_layout(design = design, widths = 1) +
      plot_annotation(
        title = 'This is the Title',
        subtitle = 'These 4 plots will reveal yet-untold secrets about our beloved data-set',
        caption = 'By Pascal Küng',
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic")
        )
      ) & theme(
        title = element_text(family = font_family),
        text = element_text(family = font_family) # Adjust size as needed
      )
  
    # Store the combined plot
    plots_list[[e]] <- combined_plot
  }  # End of effects loop
  
  return(plots_list)
}