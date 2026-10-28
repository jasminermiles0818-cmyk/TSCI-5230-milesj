#'---
#' title: "TSCI 5050: Introduction to Data Science"
#' author: 'Author One ^1^, Author Two ^1^'
#' abstract: |
#'  | Provide a summary of objectives, study design, setting, participants,
#'  | sample size, predictors, outcome, statistical analysis, results,
#'  | and conclusions.
#' documentclass: article
#' description: 'Manuscript'
#' clean: false
#' self_contained: true
#' number_sections: false
#' keep_md: true
#' fig_caption: true
#' output:
#'  html_document:
#'    toc: true
#'    toc_float: true
#'    code_folding: show
#' ---
#'

# Chunk Options & Environment Setup ----
# Configure knitr chunk options and load required libraries for data import,
# string/date manipulation, and matrix reshaping.
#+ init, echo=FALSE, message=FALSE, warning=FALSE
# init ----
# This part does not show up in your rendered report, only in the script,
# because we are using regular comments instead of #' comments
debug <- 0

knitr::opts_chunk$set(
  echo = debug > -1,
  warning = debug > 0,
  message = debug > 0,
  class.output = "scroll-20",
  attr.output = 'style="max-height: 150px; overflow-y: auto;"'
)


library(rio)
# simple command for importing and exporting
library(pander)
# format tables
#library(printr); # set limit on number of lines printed
library(dplyr)
#add dplyr library
library(lubridate)#date manipulation
library(stringr) #string manipulation
library(tidyr) #pivot_wider function
library(purrr)
options(max.print = 500)

panderOptions('table.split.table', Inf)
panderOptions('table.split.cells', Inf)



# Data Ingestion ----
# Bulk-import all CSV files from the target directory into a named list of data frames.
datalocation <- "~/Downloads/archive/"
list.files(datalocation, full.names = TRUE)


# Your two data frames
dat <- sapply(list.files(datalocation, full.names = TRUE), import, simplify = FALSE) %>%
  setNames(., basename(names(.)))


# Subset & Age Calculation ----
# Filter medical records for acute/viral conditions and calculate patient age at the time of encounter.
acute_viral <- dat$conditions.csv %>%
  
  filter(grepl("acute|viral", DESCRIPTION, ignore.case = TRUE))



# Match patients and calculate age at encounter

acute_viral_age <- acute_viral %>%
  
  left_join(dat$patients.csv %>% select(Id, BIRTHDATE),
            by = c("PATIENT" = "Id")) %>%
  
  mutate(
    START = as.Date(START),
    
    BIRTHDATE = as.Date(BIRTHDATE),
    age_at_encounter = time_length(interval(BIRTHDATE, START), unit = "years")
  )


# Acute Viral Pharyngitis Time-Series Trend----
#  Aggregate monthly occurrences of pharyngitis to run a baseline linear trend model.
temp <- filter(dat$conditions.csv,
               DESCRIPTION == "Acute viral pharyngitis (disorder)") %>%
  mutate(month = floor_date(START, unit = "month")) %>%
  group_by(month) %>%  summarise(count = n())



# Top Conditions Trend Analysis ----
# Identify conditions with the highest positive growth slope post-2023 and isolate the top 25 conditions.
conditionslope <- mutate(dat$conditions.csv, month = floor_date(START, unit = "month")) %>%
  group_by(month, CODE, DESCRIPTION) %>%
  summarise(count = n()) %>%
  group_by(CODE, DESCRIPTION) %>%
  filter(year(month) >= 2023 & length(unique(month)) > 10) %>%
  summarise(events = lm(count ~ month)$coefficients[2]) %>%
  arrange(desc(events))

plot(conditionslope$events, type = "l")
abline(v = 25, col = "green")

topconditionslope <- head(conditionslope, 25)$CODE


# Code Mapping & Vector Lookups ----
# Create a key-value lookup vector mapping condition codes to readable descriptions.
codemap <- dat$conditions.csv[c("CODE", "DESCRIPTION")] %>%
  unique() %>%
  with(setNames(DESCRIPTION, CODE))


# Matrix Preparation for Co-Occurrence ----
# Reshape patient and encounter data into binary presence/absence matrices for the top conditions.
paitentcodes <- filter(dat$conditions.csv, CODE %in% topconditionslope) %>%
  distinct(PATIENT, CODE) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = CODE,
              values_from = present,
              values_fill = 0)

encountercodes <- filter(dat$conditions.csv, CODE %in% topconditionslope) %>%
  distinct(ENCOUNTER, CODE) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = CODE,
              values_from = present,
              values_fill = 0)

npatients <- nrow(dat$patients.csv)
nencounters <- nrow(dat$encounters.csv)
codecombos <- combn(topconditionslope, 2, simplify = FALSE) #2 gives you combinations of pairs, then FALSE gives you all those possible combinations of the pairs


# Lift Metric & Co-Occurrence Analysis ----
# Define helper function to compute statistical lift between condition pairs, then visualize co-occurrence via heatmaps.
codecombos <- combn(topconditionslope, 2, simplify = FALSE)

fn_lift <- function(xx,
                    codesource = paitentcodes,
                    denom = npatients) {
  counta <- sum(codesource[[xx[1]]])
  countb <- sum(codesource[[xx[2]]])
  expected <- counta * countb / denom
  observed <- sum(codesource[[xx[1]]] * codesource[[xx[2]]])
  out <- if (expected == 0) {
    1
  } else {
    observed / expected
  }
  data.frame(cnd_a = xx[],
             cnd_b = rev(xx[]),
             lift = out)
}
paitentliftmatrix <- map(codecombos, fn_lift) %>% list_rbind() %>%
  mutate(cnd_a = codemap[cnd_a], cnd_b = codemap[cnd_b]) %>%
  xtabs(lift ~ cnd_a + cnd_b, data = .)
heatmap(
  log1p(paitentliftmatrix),
  symm = T,
  scale = "none",
  col = hcl.colors(50, "RdBu", rev = TRUE)
)

# Patient-level co-occurrence heatmap
paitentliftmatrix <- map(codecombos, fn_lift) %>%
  list_rbind() %>%
  mutate(cnd_a = codemap[cnd_a], cnd_b = codemap[cnd_b]) %>%
  xtabs(lift ~ cnd_a + cnd_b, data = .)

heatmap(
  log1p(paitentliftmatrix),
  symm = TRUE,
  scale = "none",
  col = hcl.colors(50, "RdBu", rev = TRUE)
)

# Encounter-level co-occurrence heatmap
encounterliftmatrix <- map(codecombos, fn_lift, codesource = encountercodes, denom = nencounters) %>%
  list_rbind() %>%
  mutate(cnd_a = codemap[cnd_a], cnd_b = codemap[cnd_b]) %>%
  xtabs(lift ~ cnd_a + cnd_b, data = .)

e_heat <- heatmap(
  log1p(encounterliftmatrix),
  symm = TRUE,
  scale = "none",
  col = hcl.colors(50, "RdBu", rev = TRUE)
)
systemic_conditions <- colnames(encounterliftmatrix)[e_heat$colInd][18:25]
local_conditions <- colnames(encounterliftmatrix)[e_heat$colInd][-(3:5)] %>% setdiff(systemic_conditions)



#Summarizing conditions by encounter -----
conditions_by_encounter <- group_by(dat$conditions.csv, ENCOUNTER, PATIENT) %>%
  summarize(
    local = length(intersect(local_conditions, DESCRIPTION))
    ,
    systemic = length(intersect(systemic_conditions, DESCRIPTION))
    ,
    first_condition_start = min(START, na.rm = TRUE)
    ,
    last_condition_start = max(START, na.rm = TRUE)
    ,
    conditon_type = case_when(
      local > 0 &
        systemic > 0 ~ "both",
      #When encounter meet both conditions
      local > 0 ~ "local",
      #when encounter only has local conditons
      systemic > 0 ~ "systemic",
      #when encounter only has systemic conditions
      TRUE ~ "neither"
    )
  ) #When encounter meets neither condition
#Most encounters have more than one condition, at least 4
subset(conditions_by_encounter, systemic > 0)$systemic %>% hist()

#Most encounters have 1 local condition, none have more than 3
subset(conditions_by_encounter, local > 0)$local %>% hist()

sum(conditions_by_encounter$systemic > 0)#Number of encounters with systemic conditions
sum(conditions_by_encounter$local > 0) #Number of encounters with local conditions
sum((conditions_by_encounter$systemic > 0) * (conditions_by_encounter$local >
                                                0)) #Number of both systemic and local conditions

#Numeric plot of the histogram created for local conditions
with(conditions_by_encounter, table(local))
#Numeric plot of the histogram created for systemic conditions
with(conditions_by_encounter, table(systemic))
#Contingency table of both local and systemic encounters
with(conditions_by_encounter, table(local, systemic))

#Counts of different types of encounters
conditions_by_encounter$conditon_type %>% table %>% prop.table() #prop.table changes from raw data to percentage
subset(
  conditions_by_encounter,
  first_condition_start != last_condition_start &
    conditon_type != 'neither'
) #!= shows inequality  #Make sure we can use first conditions start for it to show when the encounter took place

#Creating Conditons by Patient
conditions_by_patient <- group_by(conditions_by_encounter, PATIENT) %>%
  summarise(
    nlocal = sum(local > 0)
    ,
    nsystemic = sum(systemic > 0)
    ,
    date_systemic = if (nsystemic > 0) {
      min(first_condition_start[systemic > 0])
    } else {
      NA
    }
    %>% as.Date(),
    #date of the first systemic condition first
    date_local = if (nlocal > 0) {
      min(first_condition_start[local > 0])
    } else {
      NA
    }
    %>% as.Date(),
    #date of first local condition first
    condition_order = case_when(
      nlocal > 0 & nsystemic == 0 ~ "local only",
      nsystemic > 0 & nlocal == 0 ~ "systemic only",
      nlocal > 0 &
        nsystemic > 0 & date_local < date_systemic ~ "local first",
      nlocal > 0 &
        nsystemic > 0 &
        date_systemic < date_local ~ "systemic first",
      nlocal > 0 &
        nsystemic > 0 & date_local == date_systemic ~ "same date",
      TRUE ~ "neither"
    ),
    
    interval =  as.numeric(difftime(date_systemic, date_local, units = "days"))   #debug = {if(n_local+n_systemic > 0){browser()
    #Adding patient demographics to the condition summary
  ) %>% left_join(dat$patients.csv, by = c("PATIENT" = "Id"))

#Creating histogram for each patient----
subset(conditions_by_patient, nlocal > 0)$nlocal %>%
  hist(
    main = "Local Encounters per Patient",
    xlab = "Number of local encounters",
    ylab = "Number of patients",
    col = "lightgreen"
  )
#head(conditions_by_patient) in the console allows you to view first six patients
#names(conditions_by_patient) in the console allows you to see all available columns

# History of disease spread -----------------------------------------------
hist(conditions_by_patient$date_local, breaks = 100)

#Getting first date for local and systemic encounters
min(conditions_by_patient$date_local, na.rm = TRUE)
min(conditions_by_patient$date_systemic, na.rm = TRUE)

summarise(
  conditions_by_patient,
  datelocal = min(date_local, na.rm = TRUE)
  ,
  datesystemic = min(date_systemic, na.rm = TRUE)
)

#Grouping when the first patient encounter happened by State
group_by(conditions_by_patient, STATE) %>% summarise(
  datelocal = min(date_local, na.rm = TRUE)
  ,
  datesystemic = min(date_systemic, na.rm = TRUE)
) %>% arrange(datelocal, datesystemic)

#Getting more specific, first patient encounter within Texas by city
subset(conditions_by_patient, STATE == "Texas" &
         !is.na(date_local)) %>% group_by(CITY) %>% summarise(
           datelocal = min(date_local, na.rm = TRUE)
           ,
           datesystemic = min(date_systemic, na.rm = TRUE)
         ) %>% arrange(datelocal, datesystemic)

#Getting encounters
localencounters <- subset(conditions_by_encounter, local > 0)$ENCOUNTER
systemicencounters <- subset(conditions_by_encounter, systemic > 0)$ENCOUNTER

#Getting patient
localpatients <- subset(conditions_by_encounter, local > 0)$PATIENT %>% unique()
systemicpatients <- subset(conditions_by_encounter, systemic > 0)$PATIENT %>%
  unique()

#Procedures.csv ----
subset(dat$procedures.csv, ENCOUNTER %in% localencounters)$DESCRIPTION %>% table %>% sort
subset(dat$procedures.csv, ENCOUNTER %in% systemicencounters)$DESCRIPTION %>% table %>% sort

#Observations.csv
subset(dat$observations.csv, ENCOUNTER %in% localencounters)$DESCRIPTION %>% table %>% sort
subset(dat$observations.csv, ENCOUNTER %in% systemicencounters)$DESCRIPTION %>% table %>% sort

#Number of records within each encounter for local and systemic conditions
subset(dat$observations.csv, ENCOUNTER %in% localencounters) %>% nrow()


# Summarize a selected dataset to compare its overall records with records linked to local and systemic encounters. Count both rows and unique
# for each group. This helps show how much information is recorded
# per encounter and how that differs between the groups.
fn_rows_encounters <- function(xx)
{
  summarise(
    dat[[xx]],
    source = xx,
    rows = n(),
    rows_local = sum(ENCOUNTER %in% localencounters, na.rm = TRUE),
    rows_systemic = sum(ENCOUNTER %in% systemicencounters, na.rm = TRUE),
    encounters = length(unique(ENCOUNTER)),
    encounters_local = length(intersect(ENCOUNTER, localencounters)),
    encounters_systemic = length(intersect(ENCOUNTER, systemicencounters)),
    rows_per_encounter = rows / encounters,
    local_rows_per_encounter = rows_local / encounters_local,
    systemic_rows_per_encounter = rows_systemic / encounters_systemic
  )
}

#Counts of conditions types by patient phenotype by all datasets in dat
fn_rows_encounters("observations.csv")
names(dat) %>% setdiff(c("patients.csv", "encounters.csv")) %>% sapply(fn_rows_encounters, simplify = FALSE) %>% bind_rows()

fn_rows_patients <- function(xx)
{
  summarise(
    dat[[xx]],
    source = xx,
    rows = n(),
    rows_local = sum(PATIENT %in% localpatients, na.rm = TRUE),
    rows_systemic = sum(PATIENT %in% systemicpatients, na.rm = TRUE),
    patients = length(unique(PATIENT)),
    patients_local = length(intersect(PATIENT, localpatients)),
    patients_systemic = length(intersect(PATIENT, systemicpatients)),
    rows_per_patients = rows / patients,
    local_rows_per_patients = rows_local / patients_local,
    systemic_rows_per_patients = rows_systemic / patients_systemic
  )
}

#Counts of conditions types by patient phenotype by all datasets in dat
fn_rows_patients("observations.csv")
names(dat) %>% setdiff(c("patients.csv", "encounters.csv")) %>%
  sapply(fn_rows_patients, simplify = FALSE) %>% bind_rows()


# Extract medications for each patient and encounter, excluding missing medication codes.
transmute(dat$medications.csv,
          PATIENT,
          ENCOUNTER,
          CODE,
          DESCRIPTION,
          source = "medications") %>% filter(!is.na(CODE))
# Extract prescription reasons for each patient and encounter, excluding missing reason codes.
transmute(
  dat$medications.csv,
  PATIENT,
  ENCOUNTER,
  CODE = REASONCODE,
  DESCRIPTION = REASONDESCRIPTION,
  source = "medications_reason"
) %>% filter(!is.na(CODE))


fn_standardize <- function(xx) {
  if ("SSN" %in% names(xx)){return(NULL)}
  if ("Id" %in% names(xx) & !"ENCOUNTER" %in% names(xx)) {xx <- mutate(xx,ENCOUNTER=Id)}
  out <- transmute(xx,
                   PATIENT,
                   ENCOUNTER,
                   CODE= as.character(CODE),
                   DESCRIPTION,
                   subtype = "CODE") %>% filter(!is.na(CODE))
  if (length(intersect(colnames(xx), c("REASONCODE", "REASONDESCRIPTION"))) ==
      2) {
    #alternative if(all(c("REASONCODE","REASONDESCRIPTION") %in% names(xx))
    out <- bind_rows(
      out,
      transmute(
        xx,
        PATIENT,
        ENCOUNTER,
        CODE = as.character(REASONCODE),
        DESCRIPTION = as.character(REASONDESCRIPTION),
        subtype = "reason"
      ) %>% filter(!is.na(CODE))
    )}
  #browser()
  return(out)
}

#Create a fact table of standardized codes and descriptions
fact_table <- sapply(dat, fn_standardize, simplify= FALSE) %>% bind_rows(.id= "type") %>% filter(ENCOUNTER %in% dat$encounters.csv$Id)

#Determine the number of encounters with a specified code
summarise(fact_table, encounter_no = length(unique(ENCOUNTER)),
  encounter_w_code =length(unique(ENCOUNTER[CODE == 24199005])),
  systemic=length(systemicencounters),
  systemic_w_code= length(intersect(ENCOUNTER[CODE=="24199005"],systemicencounters)))
  
#Create a function to determine number of encounters for each code (medication, allergy, care plan, etc.)
fn_encounter_codes <- function(xx){
  encounter_w_code <- length(unique(fact_table$ENCOUNTER[fact_table$CODE == xx]))
  systemic_w_code <- length(intersect(fact_table$ENCOUNTER[fact_table$CODE==xx], systemicencounters))
  out <- c(systemic_w_code = systemic_w_code, encounter_w_code= encounter_w_code)
  return(out)
  }

#test if fn_encounter_code works on a test code (agitation)
fn_encounter_codes("24199005")

# Apply fn_encounter_codes to all unique codes from fact_table
encounter_code_tallies <- sapply(unique(fact_table$CODE),fn_encounter_codes, simplify=FALSE) %>% bind_rows %>% View()

#Assignment: find additional columns in one of the df in dat, not captured in fn_standardize
#Examining the first recorded allergy reaction for each patiient and encounter
dat$allergies.csv %>% 
  select(PATIENT, ENCOUNTER, REACTION1, DESCRIPTION1, SEVERITY1)%>%
  filter(!is.na(REACTION1)) %>% 
  head()
#Examining the second recorded allergy 
dat$allergies.csv %>% 
  select(PATIENT, ENCOUNTER, REACTION2, DESCRIPTION2, SEVERITY2)%>%
  filter(!is.na(REACTION2)) %>% 
  head()
# Count how many allergy records have a first and second reaction code
dat$allergies.csv %>%
  summarise(
    first_reaction = sum(!is.na(REACTION1)),
    second_reaction = sum(!is.na(REACTION2))
  )
# Standardize codes and descriptions from each data frame
fn_standardize <- function(xx) {
  
  # Skip the patient demographic data frame
  if ("SSN" %in% names(xx)) {
    return(NULL)
  }
  
  # Use Id as the encounter identifier when ENCOUNTER is not present
  if ("Id" %in% names(xx) & !"ENCOUNTER" %in% names(xx)) {
    xx <- mutate(xx, ENCOUNTER = Id)
  }
  
  # Collect the main code and description from the data frame
  out <- transmute(
    xx,
    PATIENT,
    ENCOUNTER,
    CODE = as.character(CODE),
    DESCRIPTION,
    subtype = "CODE"
  ) %>%
    filter(!is.na(CODE))
  
  # Add reason codes and descriptions when those columns exist
  if (all(c("REASONCODE", "REASONDESCRIPTION") %in% names(xx))) {
    out <- bind_rows(
      out,
      transmute(
        xx,
        PATIENT,
        ENCOUNTER,
        CODE = as.character(REASONCODE),
        DESCRIPTION = as.character(REASONDESCRIPTION),
        subtype = "reason"
      ) %>%
        filter(!is.na(CODE))
    )
  }
  
# Testing function on allergies -----
  # Add the first allergy reaction when those columns exist
  if (all(c("REACTION1", "DESCRIPTION1") %in% names(xx))) {
    out <- bind_rows(
      out,
      transmute(
        xx,
        PATIENT,
        ENCOUNTER,
        CODE = as.character(REACTION1),
        DESCRIPTION = as.character(DESCRIPTION1),
        subtype = "allergy_reaction_1"
      ) %>%
        filter(!is.na(CODE))
    )
  }
  
  # Add the second allergy reaction when those columns exist
  if (all(c("REACTION2", "DESCRIPTION2") %in% names(xx))) {
    out <- bind_rows(
      out,
      transmute(
        xx,
        PATIENT,
        ENCOUNTER,
        CODE = as.character(REACTION2),
        DESCRIPTION = as.character(DESCRIPTION2),
        subtype = "allergy_reaction_2"
      ) %>%
        filter(!is.na(CODE))
    )
  }
  
  # Return all the standardized rows collected above
  return(out)
}
# Test the function on allergies only
allergy_test <- fn_standardize(dat$allergies.csv)

# Count the rows produced for each type of information
table(allergy_test$subtype)
# Apply the updated function to all data frames and combine the results
fact_table <- sapply(dat, fn_standardize, simplify = FALSE) %>%
  bind_rows(.id = "type")

#Group records by code and description
#overall and within the systemic and local encounter lists, remove grouping after calc.
#Add total number of encounters for comparsion 
encounter_code_tallies_2 <- fact_table %>%
  group_by(CODE, DESCRIPTION, type) %>%
  summarise(encounter_w_code = n_distinct(ENCOUNTER),
            systemic_w_code = sum(unique(ENCOUNTER) %in% systemicencounters),
            local_w_code = sum(unique(ENCOUNTER) %in% localencounters)
  ) %>% 
  mutate(total_encounters=nrow(dat$encounters.csv),
         total_systemic_encounters = length(systemicencounters),
         total_local_encounters = length(localencounters),
         fraction_total_w_code= encounter_w_code/total_encounters,
         fraction_systemic_w_code= systemic_w_code/total_systemic_encounters,
         fraction_local_w_code= local_w_code/total_local_encounters,
         local_no_code = total_local_encounters - local_w_code,
         nonlocal_w_code = encounter_w_code - local_w_code,
         nonlocal_no_code = total_encounters - encounter_w_code - total_local_encounters + local_w_code,
         #Haldane corrections
         haldane_local_w_code = local_w_code + 0.5,
         haldane_local_no_code = local_no_code + 0.5,
         haldane_nonlocal_w_code = nonlocal_w_code + 0.5,
         haldane_nonlocal_no_code = nonlocal_no_code + 0.5,
         #Proportions
         local_code_prevalence= haldane_local_w_code/(haldane_local_w_code+haldane_local_no_code),
         nonlocal_code_absence=  haldane_nonlocal_no_code/(haldane_nonlocal_w_code+haldane_nonlocal_no_code),
         local_likelihood_ratio= local_code_prevalence/(1-nonlocal_code_absence),
         local_OR= haldane_local_w_code*haldane_nonlocal_no_code/(haldane_local_no_code*haldane_nonlocal_w_code)
         ) %>% arrange(desc(local_likelihood_ratio)
         )
         

# |           | Has Code                                          | No Code                                                                                   |
# |-----------|------------------------------------------------   |-------------------------------------------------------------------------------------------|
# | Is Local  | local_w_code                                      | local_no_code = total_local_encounters - local_w_code                                      |
# | Not Local | nonlocal_w_code = encounter_w_code - local_w_code | nonlocal_no_code = total_encounters - encounter_w_code - total_local_encounters + local_w_code |
filter(fact_table, CODE == "59408-5")
filter(dat$observations.csv, CODE == "59408-5") %>% with(browser())



# Validating total encounters between dataframes
n_distinct(fact_table$ENCOUNTER) # number of encounters from fact_table
n_distinct(dat$encounters.csv$Id) # number of encounters from encounters.csv where Id is ENCOUNTER
discrepant_encounters <- setdiff(fact_table$ENCOUNTER,c(dat$encounters.csv$Id,""))

#plot(conditionslope$events, type="l") #plot the ranked slopes to visually spot cutoff/elbow
#abline(v=25,col="blue") #25 reasonable cutoff


