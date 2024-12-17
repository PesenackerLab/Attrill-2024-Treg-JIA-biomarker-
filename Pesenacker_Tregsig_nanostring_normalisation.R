# # Install packages, if necessary
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("tidyverse", "here", "janitor", "NACHO", "formattable"))

# Load Packages----
library(tidyverse)
library(janitor)
library(formattable)
library(NACHO)

### Change path to find the data
path <- "path/to/primary/directory"
setwd(path)

rccdata <- load_rcc('./data/raw_rcc',
                    './data/annotations.csv',
                    'file_name',
                    housekeeping_genes = NULL,
                    housekeeping_predict = FALSE,
                    housekeeping_norm = FALSE,
                    normalisation_method = "GEO",
                    n_comp = 5)

# Graphs for the different QC and Norm parameters
# However, careful as limits might be different
QCparameters <- visualise(rccdata)

# Added extra column for positive control normalisation flags
positive_control_flag <- (rccdata[["nacho"]][["Positive_factor"]] > 3 | rccdata[["nacho"]][["Positive_factor"]] <0.3)

poscontrolflag <- rccdata[["nacho"]] %>%
  mutate(positive_control_flag)

# Visualise table in Viewer - flags come up as red- NOT WORKING, NEED TO FIX
poscontrolflag <- formattable(poscontrolflag, positive_control_flag = color_tile('green','red'))

# Write CSV with flags
write_csv(poscontrolflag, './output/poscontrolflag.csv')


# Write Csv----
csvdata <- write_csv(rccdata[["nacho"]], 
                     './output/normdata.csv')

# Data in readable format
posnormdata <- pivot_wider(csvdata,
                   names_from = 'Name', 
                   values_from = Count_Norm, 'file_name')

write_csv(posnormdata, './output/posnormdata.csv')

# Total sum normalisation----

# Remove control data + added row with total
# Change column index here to only those of controls
totalcount <- posnormdata %>%
  select(-c(44,45,46,47,48,49,50,51,52,53,54,55,56,57))
  
total <- rowSums(totalcount[,-1])

totalcount <- totalcount %>%
  mutate('total' = total)

str(totalcount)

# Calculate normalisation factor
norm_factor <- 5000 / totalcount$total
total * norm_factor

# Multiply values by norm_factor
total_norm <- norm_factor * totalcount[,-1]


# Add normalisation factor to table
total_norm <- total_norm %>%
  mutate('norm_factor' = norm_factor)

# Total count normalisation flag
totalcount_norm_flag <- (total_norm$norm_factor > 10 |total_norm$norm_factor <0.1)

total_norm <- total_norm %>%
  mutate(totalcount_norm_flag = totalcount_norm_flag) %>%
  mutate('file_name' = totalcount$file_name)

# Add annotations----
annot <- read_csv('./data/annotations.csv')
total_norm <- total_norm %>%
  mutate (annot)

# Add information about positive flags
poscontrol <- poscontrolflag %>%
  select(file_name, Positive_factor, positive_control_flag)
  
full_norm <- total_norm %>%
  inner_join(poscontrol, by = 'file_name')

# Remove duplicates due to the way the table was arranged at the start
full_norm <- full_norm[!duplicated(full_norm), ]

# Export file with both normalisation data
full_norm <- write_csv(full_norm, ('./output/full_norm.csv'))

# Remove samples that don't pass QC----
no_flag <- full_norm %>% slice(-c(22, 40, 51, 52, 96, 98, 115, 117, 149))

# Export file without any of the flagged up samples
write_csv(no_flag, ('./output/no_flag.csv'))

#Log2 transformation----

# Alternative approach to log2 transformation: for each entry in the no_flag table, add 1 before applying log2.
# this avoids counts of 0 and of 2 mapping to 1.

# Convert 0 to 1, to keep absolute zero counts after applying log.
log2_no_flag <- read_csv(('./output/no_flag.csv'))
log2_no_flag[log2_no_flag == 0] <- 1

# log2 tranform data
log2_no_flag <- rapply(log2_no_flag, f = log2, classes = c("numeric", "integer"), how = "replace")
log2_noflag_data <- log2_no_flag %>% select(-c(49,50,51,54,55))

write_csv(log2_noflag_data, './output/log2_noflag_data.csv')
