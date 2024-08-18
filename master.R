### The impact of extreme heat on cardiovascular morbidity in Herefordshire and Worcestershire, England
### McEvoy L 2024

## PRELIMS -------------

# Clear environment
rm(list = ls())

# Set working directory
setwd("C:\\Users\\sirsa\\OneDrive\\Documents\\2024Mcevoy")


# Install packages
list.of.packages <- c("dplyr", "lubridate", "survival", "tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
library(lubridate)
library(dplyr)
library(survival)
library(tidyr)

## READ IN DATA ---------

# Clear patient data

data <- read.csv("example_extract.csv") %>%
  # Raname required columns
  rename(
    id = Event.ID,
    age = STARTAGE_CALC,
    age.group = Quinary.Age.Group..Text.,
    sex = SEX,
    date.event = EPISTART,
    icd.code = DIAG_3_01,
    diagnosis = SHORT.DESCRIPTION..3.DIGIT.,
    county = RESCTY_ONS
  ) %>%
  # Select required columns
  select(id, age, age.group, sex, date.event, icd.code, diagnosis, county) %>%
  # Convert date to date format
  mutate(date.event = dmy(date.event)) %>%
  # COnvert ID to text format
  mutate(id = as.character(id)) %>%
  # Label county
  mutate(county = ifelse(county == "E10000034", "Worcestershire", ifelse(county == "E99999999", "Herefordshire", "Other"))) %>%
  filter(county != "Other") %>%
  # Label sex
  mutate(sex = ifelse(sex == 1, "Male", ifelse(sex ==2, "Female", "Unknown")))  %>%
  filter(sex != "Unknown") %>%
  # Rescale age so the analysis gives the OR for every 10 years
  mutate(age.scaled = age / 10)


# Ensure age group as factor level in appropriate order
age.levels <- c("00-04", "05-09", "10-14", "15-19", "20-24", "25-29", 
                "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", 
                "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90+")

data$age.group <- factor(data$age.group, levels = age.levels, ordered = TRUE)

# Ensure sex and county are as factors
data$sex <- factor(data$sex, levels = c("Male", "Female", "Unknown"))
data$county <- factor(data$county, levels = c("Worcestershire", "Herefordshire", "Other"))  

# Get heatwave periods

heatwave <- read.csv("heatwave_periods.csv") %>%
  # Select required columns
  select(1,2) %>%
  # Convert date to date format
  mutate(date = dmy(date)) 

# Define ICD disease subgroups

other.heart.disease <- paste0("I", seq(30,52)) # Other forms of heart disease
hypertension <- paste0("I", seq(10,15)) # Hypertensive diseases
ihd <- paste0("I", seq(20,25))  # Ischaemic heart diseases
stroke <- paste0("I", seq(60,69))  # Cerebrovascular diseases
  
# Filter for disease subgroup (unhash the line and edit as needed)

# data <- data %>% filter(icd.code %in% hypertension) # Options hypertension, other.heart.disease, ihd, stroke 

## PROCESS DATA ---------

# Function to create control periods based on the same day of the week, month, and year
create_control_windows <- function(event_date) {
  # Get all dates in the same month and year as the event date
  start_of_month <- floor_date(event_date, "month")
  end_of_month <- ceiling_date(event_date, "month") - days(1)
  same_month_dates <- seq(from = start_of_month, to = end_of_month, by = "day")
  
  # Select dates that are the same day of the week as the event date
  control_dates <- same_month_dates[weekdays(same_month_dates) == weekdays(event_date)]
  
  # Exclude the event date itself
  control_dates <- control_dates[control_dates != event_date]
  
  return(control_dates)
}

# Apply function to create a long data frame with event and control dates
expanded.data <- data %>%
  rowwise() %>%
  mutate(
    date = list(c(date.event, create_control_windows(date.event))),
    case_control = list(c("Case", rep("Control", length(create_control_windows(date.event)))))
  ) %>%
  unnest(cols = c(date, case_control)) %>%
  ungroup() %>%
  select(id, date, age.scaled, age.group, sex, county, case_control)

# Merge the heatwave data

final.data <- expanded.data %>%
  # Add heatwave dates
  left_join(heatwave, by = "date") %>%
  # If heatwave date, label as 1, else 0
  mutate(heatwave = ifelse(heat.period.status == 2, 1, 0)) %>%
  # Keep required columns
  select(-heat.period.status) %>%
  # Arrange by patient id then by date
  arrange(id, date)



## UNIVARIATE ANALYSIS -----------

# Perform conditional logistic regression
model <- clogit(case_control == "Case" ~ heatwave + strata(id), data = final.data)

# Extract the model coefficients summary
coef_summary <- summary(model)$coefficients

# Calculate the Odds Ratios (OR)
OR <- exp(coef_summary[, "coef"])

# Calculate the 95% Confidence Intervals (CI)
CI_lower <- exp(coef_summary[, "coef"] - 1.96 * coef_summary[, "se(coef)"])
CI_upper <- exp(coef_summary[, "coef"] + 1.96 * coef_summary[, "se(coef)"])

# Extract the p-values
p_values <- coef_summary[, "Pr(>|z|)"]

# Combine the OR, CI, and p-values into a data frame
OR_CI_p <- data.frame(
  Variable = rownames(coef_summary),
  OR = OR,
  CI_lower = CI_lower,
  CI_upper = CI_upper,
  p_value = p_values
)

# Print the OR, CI, and p-values
print(OR_CI_p)


## MULTIVARIATE ANALYSIS -----------

# Perform conditional logistic regression
model <- clogit(case_control == "Case" ~ heatwave + age.scaled + sex + county + strata(id), data = final.data)

# Extract the model coefficients summary
coef_summary <- summary(model)$coefficients

# Calculate the Odds Ratios (OR)
OR <- exp(coef_summary[, "coef"])

# Calculate the 95% Confidence Intervals (CI)
CI_lower <- exp(coef_summary[, "coef"] - 1.96 * coef_summary[, "se(coef)"])
CI_upper <- exp(coef_summary[, "coef"] + 1.96 * coef_summary[, "se(coef)"])

# Extract the p-values
p_values <- coef_summary[, "Pr(>|z|)"]

# Combine the OR, CI, and p-values into a data frame
OR_CI_p <- data.frame(
  Variable = rownames(coef_summary),
  OR = OR,
  CI_lower = CI_lower,
  CI_upper = CI_upper,
  p_value = p_values
)

# Print the OR, CI, and p-values
print(OR_CI_p)

# Reference groups
    # Sex vs Males
    # county vs Worcestershire
    # Age scales for every 10-year increase in age
