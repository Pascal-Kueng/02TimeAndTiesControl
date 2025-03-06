

# Load the rmarkdown package
library(rmarkdown)

# Get the full paths for all .Rmd files in that directory
rmd_files <- list.files(getwd(), pattern = "\\.Rmd$", full.names = TRUE)

#rmd_files[3: length(rmd_files)]

# Loop through each file and knit it
for (rmd in rmd_files) {
  message("Rendering: ", rmd)
  try(render(input = rmd, output_format = "all"))  
}

message("All .Rmd files have been successfully rendered.")

