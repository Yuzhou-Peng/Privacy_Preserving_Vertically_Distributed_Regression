library(readxl)
library(tidyverse)

nmda <- read_excel("D:/Research_Vertically_Federated_Learning/NMDA_with_population_sex_ratio_copy.xlsx")



# 2) Convert onset month to Date
df <- nmda %>%
  mutate(
    onset_month = as.Date(paste0(`Date of onset1`, "-01"))
  )

# 3) Automatically extract ALL chemical columns
chem_cols <- names(df)[str_detect(names(df), "M[0-3]$")]

# Quick check
chem_cols

summary_by_onset <- df %>%
  group_by(onset_month) %>%
  summarise(
    n_cases = n(),
    across(
      all_of(chem_cols),
      ~ mean(as.numeric(.x), na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  arrange(onset_month)

