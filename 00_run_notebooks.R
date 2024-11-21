


SHUTDOWN <- T



# Load necessary libraries
library(rmarkdown)
library(fs)


# Define the list of .Rmd files you want to knit
notebooks <- c(
  "01_FinalModels.Rmd"
  , "02_SensitivityAr1.Rmd"
  , "03_SensitivityPushingSeparate.Rmd"
  , "04_SensitivityExchangeProcessesCovariates.Rmd"
  , "04_SensitivityOtherCovariates.Rmd"
)


# Create a log file to record any failures, overwrite if it exists
log_file <- file(file.path("knit_failures.log"), open = "wt")
writeLines(paste("Knit failures log -", Sys.time(), "\n"), log_file)

# Function to knit Rmd file and log failures
knit_rmd <- function(file) {
  tryCatch({
    render(file)
    TRUE
  }, error = function(e) {
    writeLines(paste0("Failed to knit: ", file, "\nError: \n", e$message, "\n\n\n"), log_file)
    FALSE
  })
}

# Get list of specified .Rmd files in the directory
rmd_files <- file.path(notebooks)

# Check if the specified files exist
rmd_files <- rmd_files[file_exists(rmd_files)]

# Knit each file and log failures
results <- sapply(rmd_files, knit_rmd)

# Close the log file
close(log_file)

# shutdown sequence

if (SHUTDOWN) {
  if (.Platform$OS.type == "windows") {
    shell("shutdown -s -t 300")  # 300 seconds delay before shutdown
  } else if (.Platform$OS.type == "unix") {
    system("shutdown -h +5")  # 5 minutes delay before shutdown
  }
} else {
  shell("shutdown /a")
}

