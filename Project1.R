library(tidyverse)
library(readxl)

############ DATA EXPLORATION #############

# Load species metadata
animals <- read_excel("speciesTable.xlsx")

# 1. How many species in this dataset?
print(head(animals))
print(nrow(animals))


# 2. 
