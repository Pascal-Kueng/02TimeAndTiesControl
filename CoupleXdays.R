library(tidyverse)

df_merged <- openxlsx::read.xlsx(file.path('long.xlsx'))


usable_couples_count <- df_merged %>% 
  # 1. Filter for valid data first
  filter(locale == "de_CH") %>% 
  
  # 2. Group by Couple AND Day to find overlaps
  group_by(coupleID, day) %>% 
  
  # 3. Keep only days where 2 distinct users (both partners) are present
  filter(n_distinct(userID) == 2) %>% 
  
  # 4. Count how many valid 'shared days' each couple has
  group_by(coupleID) %>% 
  summarise(n_shared_days = n_distinct(day)) %>% 
  
  # 5. Filter for your threshold (10+)
  filter(n_shared_days >= 35) %>% 
  
  # 6. Count the remaining couples
  nrow()

print(usable_couples_count)


report::report_table(df_merged[, c('pre_interruptions', 'pre_completionTime')])

report::report_table(df_merged[, c('interruptions', 'completionTime')])
