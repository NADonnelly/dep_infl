
#Introduction --------

#I would like to see if we can do a decent job of imputing some of the missing 
#data from the ALSPAC dataset. We have - 

# - Much missing data from longitudinal dropout which causes the sample to 
#become rather biased over time. Lets say this is the missing data from CISR

# - Missing blood data relative to CISR data at F24

# We could try and impute either the blood data, or the CISR data, or both, but
#we would need to think carefully about what size data we are going to try to 
#impute up to. It might be worthwhile considering imputing the blood data
#at F24 using bloods from the F9 Olink panel - this means we would impute 
#everyone who has data from either F9 or F24

pacman::p_load(tidyverse,
               tidymodels,
               patchwork,
               crayon,
               visdat,
               mice,
               missRanger,
               plotly)

tidymodels_prefer()

source('./Final/0-Common_data.R')

# Load pre-split data -----

#Load up our pre-split data so we can focus on the training data only
d_split <- 
  read_rds(paste(data_dir,"alspac_data_split.rds",sep = '//'))


d_train <- 
  d_split |>
  training()

# MI with mice ------ 

# WE will use MICE to get multiple imputed datasets

## MI prep  ------

d_mi <- 
  d_train |>
  select(-c(id,multi_demo,eth,smfq_25,
            starts_with("smfq_c"),
            starts_with("ad_"),
            atypical,psychological,somatic,anxiety))

#We are removing ethnicity since we are not using it as a variable and nor do 
#we want to either impute it or using it in imputation.

#We are going to make sure we remove all the post-F24 variables too, and we are 
#going to remove the derived symptom scores as they are not going to be useful
#because they are made from variables already in the data

d_mi|>
  visdat::vis_miss(warn_large_data = FALSE) +
  theme(axis.text.x.top = element_text(angle = 90,size = 6))

#28.5% missingness overall

## MI with RFs ------

#Run mice with RFs - there are some theoretical reasons to prefer this, plus
#it can cope with categorical data, which we need. The combination of 50 
#imputations and 20 iterations e.g. 1000 runs, means this takes a VERY long time
#like a day. Adjust your expectations according. We use 50 imputations on the 
#heuristic than you need the same number of imputations as you have % missing 
#data, which for us varies up to 50% (that was the cut off we used to include
#variables), The 20 iterations is to aim for models that converge well.

imputed.datasets.rf <- mice(d_mi,
                            method = 'rf',
                            m = 50, maxit = 20,seed = 210823)

#Trace plots
#plot(imputed.datasets.rf)

#Inspect Density plots.
densityplot(imputed.datasets.rf, ~ hb_f24| .imp,groups = sex)
densityplot(imputed.datasets.rf, ~ il6_f24| .imp,groups = cisr_dep)

#This appears to have proudced reasoanble results

## Save MICED RF data -----

## Save miced-data using RFs
write_rds(imputed.datasets.rf, paste(data_dir,"alspac_data_split_imputed.rds",sep = '/'))


# Use the MI data =======

## Use the MICE output to determine number of imputations -----

#See the 3c MI Modeling script for full details....

#Simple model applied to the f24 blood variables
mice.models.u <- 
  with(imputed.datasets.rf, glm(formula = cisr_dep ~ alt_f24 *(sex + bmi24 + audit_c + smk24), 
                 family = binomial(link = 'logit')))

howManyImputations::how_many_imputations(mice.models.u)

#Apparently we need in the region of 40-50 imputations, which fits with our heuristic
