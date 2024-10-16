#################################################################
# ECON2120 PS3 Part II
#################################################################

# Clear workspace
rm(list = ls())

#===============================================================#
# Setup
#===============================================================#

# Change globals as needed
file_dir <- "C:/Users/austi/Desktop/PhD/Coursework/ECON2120 Metrics I/Problem Sets/PS3/"

# Import libraries
library(tidyverse)
library(R.matlab)
library(fixest)

# Load data 
input_file <- paste(file_dir, "firms13.mat", sep = "")
raw_data <- readMat(input_file) %>% as.data.frame()

# Variable labels for tables
table_dict = c(
  "firm_id" = "Firm ID",
  "ldsa" = "Log Output",
  "lemp" = "Log Labor",
  "lcap" = "Log Capital",
  "ldsa_diff1" = "Log Output (First Diff.)",
  "lemp_diff1" = "Log Labor (First Diff.)",
  "lcap_diff1" = "Log Capital (First Diff.)",
  "ldsa_dev" = "Log Output (Deviations)",
  "lemp_dev" = "Log Labor (Deviations)",
  "lcap_dev" = "Log Capital (Deviations)"
)

# Function to export regression tables
# Customize table style here
etable_export <- function(regression_model){
  
  # Display regression model in console
  print(etable(regression_model, dict = table_dict))
  
  # Export regression table
  etable(
    regression_model, 
    dict = table_dict, # table_dict is a global variable defined in this section
    file = paste(
      file_dir, # file_dir is a global variable defined in this section
      deparse1(substitute(regression_model)), # name of model object will be name of file
      ".tex", 
      sep = ""
    ),
    style.tex = style.tex("aer", stats.title = "\\midrule"),
    tabular = "X",
    signif.code = NA,
    digits = 3,
    digits.stats = 3,
    replace = TRUE
  )  
}

#===============================================================#
# Q1
#===============================================================#

# Reshape data from the firm level to firm-X-year level
q1_data <- raw_data %>%
  ungroup() %>%
  # Add firm identifier 
  mutate(firm_id = row_number()) %>%
  # Reshape non-identifier columns from the firm level to the firm-x-year level
  pivot_longer(
    cols = -c("firm_id"),
    names_to = c(".value", "year"),
    names_sep = "\\." # note: need to escape full stop when declaring it as separator
  ) %>%
  # Year as a numeric variable 
  mutate(year = as.numeric(year))

# Assert that we have 441 firms x 13 years (fully balanced)
stopifnot(length(unique(q1_data$firm_id)) == 441)
stopifnot(length(unique(q1_data$year)) == 13)
stopifnot(nrow(q1_data) == 441 * 13)

# Regress log output (ldsa) on a constant, log labor (lemp), and log capital (lcap)
# Variables are already in logs; constant is included by default
q1_regression <- feols(fml = ldsa ~ lemp + lcap, data = q1_data)
etable_export(q1_regression)

#===============================================================#
# Q2
#===============================================================#

q2_data <- q1_data %>%
  ungroup() %>%
  # Sort dataset by firm ID and year
  arrange(firm_id, year) %>%
  # Within each firm, construct 
  # (1) t-1 lag
  # (2) first difference
  # of each variable  
  group_by(firm_id) %>%
  mutate(
    across(
      .cols = all_of(c("ldsa", "lcap", "lemp")),
      .fns = list(
        lag1 = dplyr::lag, 
        diff1 = ~ . - dplyr::lag(.)
      ), 
      .names = "{.col}_{.fn}"
    )
  ) %>%
  ungroup()

# Regress first-differenced log output (ldsa_diff1) on a constant, 
# first-differenced log labor (lemp_diff1), and first-differenced log capital (lcap_diff1)
q2_regression <- feols(fml = ldsa_diff1 ~ lemp_diff1 + lcap_diff1, data = q2_data)
etable_export(q2_regression)

#===============================================================#
# Q4
#===============================================================#

q4_data <- q1_data %>%
  ungroup() %>%
  # Within each firm, construct 
  # (1) mean
  # (2) deviation from mean
  # of each variable
  group_by(firm_id) %>%
  mutate(
    across(
      .cols = all_of(c("ldsa", "lcap", "lemp")),
      .fns = list(
        mean = mean, 
        dev = ~ . - mean(.)
      ), 
      .names = "{.col}_{.fn}"
    )
  ) %>%
  ungroup()

# Regress deviations in log output (ldsa_dev) on a constant, 
# deviations in log labor (lemp_dev), and deviations in log capital (lcap_dev)
q4_regression <- feols(fml = ldsa_dev ~ lemp_dev + lcap_dev, data = q4_data)
etable_export(q4_regression)

# Confirm equivalence of point estimates from 
# the within-transformation regression and a regression with firm fixed effects
q1_regression_with_firm_fes <- feols(fml = ldsa ~ lemp + lcap | firm_id, data = q1_data)
etable(q1_regression_with_firm_fes, dict = table_dict)
