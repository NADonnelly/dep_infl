
#Introduction --------

#I would like to see if we can do a decent job of imputing some of the missing 
#data from the ALSPAC dataset. This could potentially go in multiple directions
#because we have

# - Much missing data from longitudinal dropout which causes the sample to 
#become rather biased over time. Lets say this is the missing data from CISR

# - Missing blood data relative to CISR data at F24

# We could try and impute either the blood data, or the CISR data, or both, but
#we would need to think carefully about what size data we are going to try to 
#impute up to. It might be worthwhile considering imputing the blood data
#up to the full CISR set at least

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

# Load data -----

#Load up our data
d <- 
  read_rds(paste(data_dir,"alspac_data_final.rds",sep = '//'))

# F9 - F24 correlations -----

#How much do the variables correlate? Useful to know for multiple imputation
d_corr_pairs <-
    bind_rows(
      d_9 |>
        mutate(dataset = "f9"),
      d_24 |>
        mutate(dataset = "f24"))  |>
    pivot_longer(-c(id,dataset)) |>
    drop_na(value) |>
    nest_by(name) |>
    filter(dim(data)[1] > 4000) |>
    mutate(
      m1 =
        list(
          data |>
            pivot_wider(names_from = dataset,values_from = value) |>
            select(-id) |>
            lm(data = _,
               formula = f24 ~ f9) |>
            broom::glance()
        )
    )

  d_corr_pairs |>
    select(name,m1) |>
    unnest(cols = m1) |>
    ungroup() |>
    ggplot(aes(adj.r.squared)) +
    geom_histogram(binwidth = 0.01)


#So most correlations are small, but all but 4/103 measures available at both timepoints have significant correlations
d_corr_pairs |>
  select(name,m1) |>
  unnest(cols = m1) |>
  ungroup() |>
  mutate(p.adj = p.adjust(p.value,method = "fdr")) |>
  arrange(-p.adj)

#The 4 being tf4ebp1, stambp, axin1 and il1alpha

#Some have rather decent correlations; 15 have adj r2 > 0.2
d_corr_pairs |>
  select(name,m1) |>
  unnest(cols = m1) |>
  ungroup() |>
  filter(adj.r.squared > 0.2)

#What is the mean and sd of the r2?
d_corr_pairs |>
  select(name,m1) |>
  unnest(cols = m1) |>
  ungroup() |>
  summarise(mu = mean(adj.r.squared),
            sig = sd(adj.r.squared))

#Mean is 9.7% and SD is 9.1%

d_0 |>
  select(IL6_F9,IL6_F24) |>
  mutate(across(IL6_F9:IL6_F24,
                ~case_when( . < 0 ~ NA_integer_ ,
                            TRUE ~ .))) |>
  drop_na(IL6_F24)

#1829 had IL6 measured at both time-points
#2999 measurements at F9
#3019 measurements at F24



# MI with mice ------ 

# Mice gives you multiple datasets to apply Rubin's rules to

## MI prep  ------

d_mi <- 
  d |>
  select(-c(id,eth,smfq_25,ph_type,
            starts_with("smfq_c"),
            starts_with("ad_"),
            atypical,psychological,somatic,anxiety))

#We are removing ethnicity since we are not using it as a variable and nor do 
#we want to either impute it or using it in imputation.

#We are going to remove the ph_type variable as it is so closely linked to the ph
#variable

#We are going to make sure we remove all the future variables too, and we are 
#going to remove the derived symptom scores as they are not going to be useful
#because they are made from variables already in the data

d_mi|>
  visdat::vis_miss(warn_large_data = FALSE) +
  theme(axis.text.x.top = element_text(angle = 90,size = 6))

#28.3% missingness overall

## MI with RFs ------

#Run mice with RFs - there are some theoretical reasons to prefer this, plus
#it can cope with categorical data, which we need. The combination of 50 
#imputations and 20 iterations e.g. 1000 runs, means this takes a VERY long time -
#like a day. Adjust your expectations according. We use 50 imputations on the 
#heuristic than you need the same number of imputations as you have % missing 
#data, which for us varies up to 50% (that was the cut off we used to include
#variables), The 20 iterations is to aim for models that converge well.

imputed.datasets.rf <- mice(d_mi,
                            method = 'rf',
                            m = 50, maxit = 20,seed = 210823)

#Trace plots
#plot(imputed.datasets.rf)

#Density plots.
densityplot(imputed.datasets.rf, ~ hb_f24| .imp,groups = sex)
densityplot(imputed.datasets.rf, ~ il6_f24| .imp,groups = cisr_dep)

#This appears promising

## Save MICED RF data -----

## Save miced-data using RFs
#write_rds(imputed.datasets.rf, paste(data_dir,"alspac_data_imputed.rds",sep = '/'))
write_rds(imputed.datasets.rf, "./Models/alspac_data_imputed.rds")

# Use the MI data =======

## Use the MICE output to determine number of imputations -----

#See the 3c MI Modeling script for full details....

#Simple model applied to the f24 blood variables
mice.models.u <- 
  with(imputed.datasets.rf, glm(formula = cisr_dep ~ il6_f24 *(sex + bmi24 + audit_c + smk24 + ph), 
                 family = binomial(link = 'logit')))

howManyImputations::how_many_imputations(mice.models.u)

#Apparently we need in the region of 11 imputations?? 
