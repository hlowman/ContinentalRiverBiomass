## 206 rivers data output from Teton
## Created: April 22, 2022
## Heili Lowman

# The following code will do some preliminary processing of the output
# of the 206 site dataset sent to Teton in April 2022.

# The visualization code has been moved to TetonOutput_exploration_re.R.

# Load packages
lapply(c("calecopal", "cowplot",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan", "here",
         "ggrepel", "patchwork", "grid","gridExtra"), require, character.only=T)

#### Data Import & Processing ####

### NOTE: THIS CODE MUST BE RUN ON THE SERVER.
### If I load this on my desktop, it crashes (cannot allocate vector size xyz).

# Import Teton run results.
#data_out <- readRDS("data_teton/teton_206rivers_output_Ricker_2022_04_26.rds") # doesn't work even on server
data_out43 <- readRDS("data_teton/teton_43rivers_output_Ricker_2022_04_27.rds") # took ~20 minutes
###

# Examine for chains mixing - too big doesn't work
# launch_shinystan(data_out43)

# Pulling out parameters separately for now.

# Extract only r data from the model
#data_out_r <- extract(data_out, c("r", "rsite"))
data_out_r43 <- extract(data_out43, c("r", "rsite"))
# so, this is an array in the format (rows, columns, matrices)

# Export data.
#saveRDS(data_out_r, "data_working/teton_206rivers_r_all_iterations_042922.rds")
saveRDS(data_out_r43, "data_working/teton_43rivers_r_all_iterations_042922.rds")

# Extract only c data from the model
#data_out_c <- extract(data_out, c("c", "csite"))
data_out_c43 <- extract(data_out43, c("c", "csite"))

# Export data.
#saveRDS(data_out_c, "data_working/teton_206rivers_c_all_iterations_042922.rds")
saveRDS(data_out_c43, "data_working/teton_43rivers_c_all_iterations_042922.rds")

# End of script.
