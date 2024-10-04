
# Data Preparation

prepare_data <- function(df, recode_pushing = TRUE, use_mi = FALSE) {
    
  
  
  predictors_to_center <- c("persuasion",
                            'pressure',
                            'pushing', 
                            'weartime',
                            'barriers',
                            'support',
                            'comf',
                            'reas')
  
  
  predictors_not_to_center <- c('day', 'userID', 'coupleID', # independent variables that don't need to be centered.
                                'aff', # outcomes (not centered)
                                'pa_sub',
                                'pa_obj',                               
                                'reactance', 
                                'plan',
                                'isWeekend',
                                'studyGroup',
                                'got_JITAI',
                                'skilled_support',
                                'ss_pa'
  )
  
  
  
  all_variables <- c(predictors_not_to_center, predictors_to_center)
  
  cat('Constructing scales\n')
  df <- df %>% 
    mutate(ss_affect1 = ss_affect1 + 1, # re-code items from 0-5 to 1-6 scale.
           ss_affect2 = ss_affect2 + 1, 
           ss_affect3 = ss_affect3 + 1, 
           ss_affect4 = ss_affect4 + 1, 
           aff = (ss_affect1 + ss_affect3 + (7 - ss_affect2) + (7 - ss_affect4)) / 4, #Affect Scale
           
           pa_sub = case_when(
             is.na(ss_pa) ~ NA, # If ss_pa is NA, pa_min_total becomes NA
             ss_pa == 0 ~ 0,          # If ss_pa is 0, pa_min_total becomes 0
             ss_pa == 1 ~ ss_pa_min_collab + ss_pa_min_solo # If ss_pa is 1, it retains its value
           ),
           
           pa_obj = ifelse(wear_minutes > 600, minutes_mvpa_non_filtered, NA),
           
           day = day/54,
           
           barriers = (ss_barr_1 + ss_barr_2 + ss_barr_3 + ss_barr_4 + ss_barr_5 + ss_barr_6 + ss_barr_7) / 7,
           
           plan = ifelse(
             is.na(ss_pa_no), 
             ifelse(ss_pa_yes_0 == 1 | ss_pa_yes_1 == 1 | ss_pa_yes_2 == 1 | ss_pa_yes_3 == 1, 1,0),
             ifelse(ss_pa_no ==1, 1, 0)
           ), 
           
           got_JITAI = ifelse(((userID %% 2 != 0) & (trig_target_slot_0 %in% c("A", "C") | trig_target_slot_1 %in% c("A", "C") | trig_target_slot_2 %in% c("A", "C"))) |
                                ((userID %% 2 == 0) & (trig_target_slot_0 %in% c("B", "C") | trig_target_slot_1 %in% c("B", "C") | trig_target_slot_2 %in% c("B", "C"))), 1, 0),
           
           got_JITAI = ifelse(is.na(got_JITAI), 0, got_JITAI),
           
           skilled_support = ifelse(
             (studyGroup == 3 & ((day * 54 + 1) > 28)) | (studyGroup != 3 & ((day * 54 + 1) > 7)), 
             1, 0
           )
           
    )%>% 
    rename(persuasion = sp_psc_more,
           pressure = sp_nsc_more,
           pushing = sp_push_plan,
           reactance = ss_reactance,
           weartime = wear_minutes,
           
           support = sp_emo_pleasure,
           comf = sp_emo_comf,
           reas = sp_emo_reass
    ) %>%
    arrange(
      coupleID,
      day
    )
  
  
  
  
  df_full <- df
  df <- df %>%
    select(all_of(all_variables))
  
  
  # pushing can only occur if there was a plan! Therefore, we set it to zero if there was no plan. 
  if (recode_pushing) {
    cat('Re-coding pusing\n')
    df$pushing[df$plan == 0] <- 0
    
    df_full$pushing[df_full$plan == 0] <- 0
  }

  
  
  
  
  
  ## Multiple Imputation
  
  if (use_mi) {
    cat('Imputing\n')
    #* Load packages
    library(mice)
    library(miceadds)
    library(micemd)
    library(VIM)
    library(naniar)
    library(ggmice)
    library(semTools)
    
    df$coupleID <- as.integer(df$coupleID)
    
    # Inspect missing patterns
    md.pattern(df)
    aggr(df, numbers = TRUE, sortVars = TRUE)
    
    # Percentage of missing data 
    missing <- pct_miss(df) / 100
    
    # Calculate number of imputed datasets according to the quadaradic rule 
    numImp <- 1 + 0.5*(missing / 0.05)^2
    numImp <- ceiling(numImp)
    
    
    # Define predictor matrix. Set hh_id2 as Level-2 cluster variable 
    pred <- quickpred(df)
    
    # include random effects for all variables. 
    pred[pred == 1] <- 2
    
    # Define clustering variable: couple ID 
    pred[ ,"coupleID"] <- -2
    plot_pred(pred)
    
    # Define methods
    meth <- make.method(df)
    meth[meth == 'pmm'] <- '2l.lmer' #2l.lmer
    meth["plan"] <- '2l.bin'
    
    #* Run Multiple imputation
    imp <- mice::mice(
      data = df, 
      pred = pred, 
      method = meth,
      m = numImp, maxit = 30,
      seed = 7777)              
    
    #plot(imp)
    densityplot(imp)
    
    
    # Create list of imputed data sets
    implist <- mice::complete(imp, action = "all")
  } else {
    implist <- NULL
  }
  
  
  
  ## Rehsape Data 
  cat('reshaping data (4field)\n')
  reshape_4field <- function(data, cluster, time, transform) {
    #' Data = dataframe, 
    #' cluster = string, Name of cluster variable, 
    #' time = string, name of time variable
    #' transform = list of strings, names of variables to reshape
    
    library(tidyverse)
    
    result <- as_tibble(data) %>%
      group_by(across(all_of(c(cluster, time)))) %>%
      mutate(across(all_of(transform), 
                    list(
                      self = ~.x,
                      partner = ~if_else(row_number() == 1, last(.x), first(.x))
                    ),
                    .names = "{.col}_{.fn}"
      )) %>%
      ungroup()
    
    return(result)
  }
  
  
  df_double <- reshape_4field(
    df, 
    "coupleID", 
    "day", 
    c("persuasion", "pressure", "pushing", 
      "weartime", "barriers", "plan", "support", "got_JITAI")
  )
  
  df_double_full <- reshape_4field(
    df_full, 
    "coupleID", 
    "day", 
    c("persuasion", "pressure", "pushing", 
      "weartime", "barriers", "plan", "support", "got_JITAI")
  )
  
  
  
  ## center Data
  cat('centering data within and between\n')
  center_df <- function(df) {
    df_centred <- df %>%
      bmlm::isolate(by = "userID",
                    value = colnames(df)[grepl('_self|_partner', colnames(df))], 
                    which = "both")
    
    return(df_centred)
  }
  
  # Prepare imputed datasets, if we imputed.
  if (use_mi) {
    for(i in 1:length(implist)){
      implist[[i]] <- center_df(implist[[i]])[[1]]
    }
  }
  
  # And prepare original dataset in the same manner
  df_double <- center_df(df_double)
  

  
  
  # Recode factors
  recode_factors <- function(dataframe) {
    dataframe$got_JITAI <- factor(dataframe$got_JITAI, levels = c(0,1), labels = c("No JITAI", "JITAI received"))
    
    if (!is.null(dataframe$got_JITAI_self)) {
      dataframe$got_JITAI_self <- factor(dataframe$got_JITAI_self, levels = c(0,1), labels = c("No JITAI", "JITAI received"))
    }
    
    # weekend as factor
    dataframe$isWeekend <- factor(dataframe$isWeekend, levels = c(0,1), labels = c("Weekday", "Weekend"))
    
    # Plan as a factor
    dataframe$plan <- factor(dataframe$plan, levels = c(0,1), labels = c("No plan", "Plan"))
    
    # Study Group as a Factor
    dataframe$studyGroup <- factor(dataframe$studyGroup - 1, levels = c(0,1,2), labels = c("Allways inerventions", "First 3 weeks interventions", "last 3 weeks interventions"))
    
    # Skilled support as factor
    dataframe$skilled_support <- factor(dataframe$skilled_support, levels = c(0,1), labels = c("Days before Intervention", "Days after Intervention"))
    
    return(dataframe)
  }
  
  df_double <- recode_factors(df_double)
  df_full <- recode_factors(df_full)
  
  
  return(list(df_double, df_full))
}
  




