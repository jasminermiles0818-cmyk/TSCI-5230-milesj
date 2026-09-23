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
#+ init, echo=FALSE, message=FALSE, warning=FALSE
# init ----
# This part does not show up in your rendered report, only in the script,
# because we are using regular comments instead of #' comments
debug <- 0;
knitr::opts_chunk$set(echo=debug>-1, warning=debug>0, message=debug>0, class.output="scroll-20", attr.output='style="max-height: 150px; overflow-y: auto;"');

library(rio);# simple command for importing and exporting
library(pander); # format tables
#library(printr); # set limit on number of lines printed
library(dplyr); #add dplyr library
library(lubridate)#date manipulation
library(stringr) #string manipulation
library(tidyr) #pivot_wider function
library(purrr)
options(max.print=500);
panderOptions('table.split.table',Inf); panderOptions('table.split.cells',Inf);
datalocation <- "~/Downloads/archive/"
list.files(datalocation,full.names = T)
#TEST <- import("C:/Users/nahid/OneDrive/Desktop/class/careplans.csv")
dat <- sapply(list.files(datalocation, full.names = TRUE), import,simplify = FALSE) %>% 
  setNames(.,basename(names(.)))
# data ingestion ----

# Your two data frames




# First, subset to acute/viral conditions

acute_viral <- dat$conditions.csv %>%
  
  filter(grepl("acute|viral", DESCRIPTION, ignore.case = TRUE))



# data ingestion ----
# Match patients and calculate age at encounter

acute_viral_age <- acute_viral %>%
  
  left_join(
    dat$patients.csv %>% select(Id, BIRTHDATE),
    by = c("PATIENT" = "Id")
  ) %>%
  
  mutate(
    
    START = as.Date(START),
    
    BIRTHDATE = as.Date(BIRTHDATE),
    
    age_at_encounter = time_length(
      
      interval(BIRTHDATE, START),unit = "years"))

  
      # Acute Viral Pharyngitis ----
temp <- filter(dat$conditions.csv,DESCRIPTION=="Acute viral pharyngitis (disorder)") %>% 
  mutate(month=floor_date(START, unit = "month")) %>% 
  group_by(month) %>%  summarise(count=n())   

lm(count~month,temp)

# topconditions----
conditionslope <-  mutate(dat$conditions.csv,month=floor_date(START, unit = "month")) %>% 
  group_by(month,CODE,DESCRIPTION) %>%  summarise(count=n())   %>% 
  group_by(CODE,DESCRIPTION) %>% filter(year(month)>=2023 & length(unique(month))>10) %>% 
  summarise(events=lm(count~month)$coefficients[2]) %>% arrange(desc(events)) 
 
plot(conditionslope$events,type="l")                                                                 
abline(v=25,col="green")
#25 seems like a reasonable cutoff for the conditions
topconditionslope<- head(conditionslope,25)$CODE #to pull out column from data frame square brackets, to pull from vector put comma before


# code/name mapping
codemap<- dat$conditions.csv [c("CODE","DESCRIPTION")]%>%unique()%>% {setNames(.$DESCRIPTION,.$CODE)} 

codemap <- dat$conditions.csv[c("CODE","DESCRIPTION")] %>% unique() %>%  # Removes all the duplicate rows 
  with(setNames(DESCRIPTION,CODE)) # no longer a dataframe, now a vector with names # curly brackets protect the period
  

#code co occurance----
#pa
paitentcodes <- filter(dat$conditions.csv,CODE %in% topconditionslope)%>% #%in% filters for a value within a vector, unique keeps all different rows,
  distinct(PATIENT, CODE) %>% 
  mutate(present=1) %>% 
  pivot_wider(names_from = CODE, values_from= present, values_fill = 0 )#piviot wider chnages the code entrys to the title and lets you know if the code is present in a patient sidenote: Change n/a to zero because you going to do math
encountercodes <- filter(dat$conditions.csv,CODE %in% topconditionslope)%>% #%in% filters for a value within a vector, unique keeps all different rows,
  distinct(ENCOUNTER, CODE) %>% 
  mutate(present=1) %>% pivot_wider(names_from = CODE, values_from= present, values_fill = 0 )
#see unique patient ids
npatients <- nrow(dat$patients.csv) #specifc to getting number of patients
nencounters <- nrow(dat$encounters.csv)
codecombos <- combn(topconditionslope,2,simplify = FALSE) #2 gives you combinations of pairs, then FALSE gives you all those possible combinations of the pairs

fn_lift <- function(xx,codesource=paitentcodes,denom=npatients){ #creates a custom in R,ac
  counta <- sum(codesource[[xx[1]]])
  countb <- sum(codesource[[xx[2]]])
  expected <- counta*countb/denom
  observed <- sum(codesource[[ xx[1]]]*codesource[[ xx[2]]])
  out <- if(expected == 0){1} else{observed/expected} 
  data.frame(cnd_a=xx[],cnd_b=rev(xx[]),lift=out) #
}
paitentliftmatrix <- map(codecombos,fn_lift) %>% list_rbind() %>% 
   mutate(cnd_a=codemap[cnd_a],cnd_b=codemap[cnd_b]) %>% 
  xtabs(lift~cnd_a+cnd_b,data=.)
heatmap(log1p(paitentliftmatrix),symm = T,scale = "none",col=hcl.colors(50, "RdBu", rev=TRUE))



encounterliftmatrix <- map(codecombos,fn_lift,codesource=encountercodes,denom=nencounters) %>% list_rbind() %>% 
  mutate(cnd_a=codemap[cnd_a],cnd_b=codemap[cnd_b]) %>% 
  xtabs(lift~cnd_a+cnd_b,data=.)
e_heat <- heatmap(log1p(encounterliftmatrix),symm = T,scale = "none",col=hcl.colors(50, "RdBu", rev=TRUE))
colnames(encounterliftmatrix)[e_heat$colInd]

rownames(paitentliftmatrix)
codemap[colnames(paitentliftmatrix)]


fn_lift(codecombos[[12]])
#we are assigning a function and naming it with fn, curly brakctes bundle a bunch of code and allow you to use it in a function, brackets pull something from a dataframe
Map() #Map takes a list and has it perform a function, actual count divided by the expected count