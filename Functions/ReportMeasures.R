
report_measures <- function(data, measures, ICC = TRUE, cluster_var = NULL) {
  
  if (ICC & any(sapply(data[,measures], is.factor))) {
    warning("Cannot add ICCs when Factors are in Dataframe.")
    ICC <- FALSE
  }
  
  if (is.null(cluster_var) & ICC == TRUE) {
    stop("Cluster Variable needed if ICC = TRUE")
  }
  
  measures_table <- as.data.frame(report::report_table(data[, measures]))
  
  measures_table <- measures_table %>%
    mutate(across(.cols = c('Mean', 'SD'),
                  .fns = ~format(round(.x, 2), nsmall = 2)))
  
  # compute percent missing. 
  if (is.null(measures_table$percentage_Missing)) {
    if (is.null(measures_table$n_Missing)) {
      measures_table$Missing <- NA
    } else {
      measures_table$Missing <- ifelse(
        is.na(measures_table$n_Missing), 
        NA,
        paste0(round(measures_table$n_Missing / measures_table$n_Obs), '%')
      )
    }
  } else {
    measures_table$Missing <- ifelse(
      is.na(measures_table$percentage_Missing),
      NA,
      paste0(round(measures_table$percentage_Missing), '%')
    )
  }
  
  measures_table$Range <- ifelse(
    is.na(measures_table$Min) | is.na(measures_table$Max),
    NA,
    paste0(format(round(measures_table$Min, 2), nsmall = 2), '-', format(round(measures_table$Max, 2), nsmall = 2))
  )
  
  
  
  finished_df <- measures_table
  
  if (ICC) {
    cors <- wbCorr(data[,measures], cluster_var)
    
    finished_df <- cbind(
      measures_table[, c('Variable', 'n_Obs', 'Missing', 'Mean', 'SD', 'Range')],
      ICC = format(round(get_icc(cors)$ICC, 2), nsmall = 2)
    )
  } else {
    measures_table$percentage_Obs <- ifelse(
      is.na(measures_table$percentage_Obs),
      NA,
      paste0(round(measures_table$percentage_Obs), '%')
    )
    
    finished_df <- measures_table[, c('Variable', 'Level', 'n_Obs', 'percentage_Obs', 'Missing', 'Mean', 'SD', 'Range')]
  }
  
  # Make them all Text for easier copying to excel. 
  finished_df[] <- lapply(finished_df, as.character)
  
  rownames(finished_df) <- NULL
  
  return(finished_df)
}
