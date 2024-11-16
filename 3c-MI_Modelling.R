
# Introduction ------

#load packages
pacman::p_load(tidyverse,
               survey,
               mice,
               broom,
               marginaleffects,
               ggdist)

source('./Final/0-Common_data.R')

#The F9 variables have a different set of covariates to adjust for compared to the F24 data
# For F9 - bmi9, ses and sex
# For F24 - sex + bmi24 + audit_c + smk24

covars = c("sex","bmi24","smk24","audit_c")

dv = 'cisr_dep'

# Use the imputated datasets to explore the cross-sectional associations at F24
imputed.datasets.rf = read_rds(paste(data_dir,"alspac_data_imputed.rds",sep = '/'))

#Reshape each imputation into a long stack
long1 <- 
  complete(imputed.datasets.rf, action ='long', include=TRUE)


midsobj2 <-  
  long1 |> 
  as_tibble() |> 
  
  #while we are here, remove the variables we aren't interested in (auxillaries)
  select(c(.imp,.id,all_of(dv), all_of(covars),contains("_f24"))) |>
  
  mutate(across(.cols = all_of(dv),.fns = ~as.integer(.x)-1)) |>
  
  #Zscore each imputation
  group_by(.imp) |>
  mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  ungroup() |>
  
  #put it back into a mids
  as.mids()



# Depression average difference -----

#We are going to apply a model where we are comparing the level of the blood
#measure between participants with ICD10 depression and those without


#Define a function and then apply to our F24 immuno-metabolic variables
fit_dep_model <- function(data,iv,dv,covars){
  

  #Try to sort out environments
  form_obj_1 <- 
    paste0(dv, ' ~ ', iv) |> 
    as.formula()
  
  mice.models.u <- 
    with(data, {
    # set form_obj's environment to the current one
    environment(form_obj_1) = environment()
    lm(form_obj_1)
  })
  
  
  form_obj_2 <- 
    paste0(dv, ' ~ ', iv, ' * (', paste0(covars[!str_detect(covars,"bmi")],collapse = ' + '),')') |> 
    as.formula()
  
  mice.models.n <- 
    with(data, {
      # set form_obj's environment to the current one
      environment(form_obj_2) = environment()
      lm(form_obj_2)
    })
  
  
  form_obj_3 <- 
    paste0(dv, ' ~ ',iv,' * (', paste0(covars,collapse = ' + '),')') |> 
    as.formula()
  
  mice.models.a <- 
    with(data, {
      # set form_obj's environment to the current one
      environment(form_obj_3) = environment()
      lm(form_obj_3)
    })
  
  
  
  
  #M1: Unadjusted Model
  # mice.models.u <- 
  #   with(data, lm(formula = paste0(dv, ' ~ ', iv) |> as.formula()))
  
  mice.results.u.tab  <-
    mice.models.u |>
    avg_comparisons(variables = iv) |> 
    tidy() |>
    mutate(model = 'M1: unadjusted') |> 
    select(-contrast) |> 
    mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 8))) |>
    as.matrix() |>
    as_tibble() |>
    mutate(across(estimate:conf.high,~as.double(.x))) 
  
  #M2: Adjusted Model (covariates minus BMI)
  # mice.models.n <- 
  #   with(data, lm(formula = paste0(dv, ' ~ ', iv, ' * (', paste0(covars[!str_detect(covars,"bmi")],collapse = ' + '),')') |> as.formula()))
  # 
  mice.results.n.tab <- 
    mice.models.n |>
    marginaleffects::avg_comparisons(variables = iv) |> 
    tidy() |>
    mutate(model = 'M2: adjusted no bmi') |> 
    select(-contrast) |> 
    mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 8))) |>
    as.matrix() |>
    as_tibble() |>
    mutate(across(estimate:conf.high,~as.double(.x))) 
  
  #M3: Fully Adjusted Model
  # mice.models.a <- 
  #   with(data, lm(formula = paste0(dv, ' ~ ',iv,' * (', paste0(covars,collapse = ' + '),')') |> as.formula()))

   mice.results.a.tab <- 
    marginaleffects::avg_comparisons(mice.models.a,variables = iv) |> 
    tidy() |>
    mutate(model = 'M3: adjusted + bmi') |> 
    select(-contrast) |> 
    mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 8))) |>
    as.matrix() |>
    as_tibble() |>
    mutate(across(estimate:conf.high,~as.double(.x))) 
   
   
   return(bind_rows( mice.results.u.tab,
                     mice.results.n.tab,
                     mice.results.a.tab))
  
  
}


#Now apply this function to our f24 variables
d_dep_avg_results = 
  tibble(vars =
           midsobj2$data |> 
           colnames() ) |> 
  filter((vars %nin% covars) & (vars %nin% dv) & !str_detect(vars,"id")) |>
  mutate(im_model = vector("list", 1)) |>
  filter(str_detect(vars,'_f24')) 

#Check we have all the variables we expect we should have
setdiff(d |> select(contains("_f24")) |> colnames(),
        d_dep_avg_results$vars)

#We have to do this in a wacky way because svyglm/mice/with seem to get very
#unhappy if you try to build a model formula using variables inside a function,
#they seem to only be able to take variables from the global environment and I 
#don't have the r skills to even formulate how to fix this, apart from changing 
#the value of iv in the global environment in each loop iteration. Sorry. 

for(i in 1:dim(d_dep_avg_results)[1]){
  
  tictoc::tic()
  dv_i = d_dep_avg_results$vars[[i]]
  
  print(paste0(i,' : ', dv_i))
  
  d_dep_avg_results$im_model[[i]] <- fit_dep_model(data   = midsobj2,
                                                   iv     = 'cisr_dep',
                                                   dv     = dv_i,
                                                   covars = covars) 
  tictoc::toc()
}



#So now we can interrogate the model output and ask what, for each variable, 
#what the average effect of going from undepressed to depressed is on the 
#blood variable of interest

d_dep_tab <- 
  d_dep_avg_results|>
  unnest(im_model) |>
  select(-term) |>
  janitor::clean_names() |>
  mutate(across(estimate:conf_high,as.double)) |> 
  relocate(vars,model)
  
#write_csv(d_dep_tab,paste(data_dir,"alspac_imputed_results.csv",sep = '/'))
write_csv(d_dep_tab,"./Models/alspac_imputed_results.csv")

d_dep_tab <-
  read_csv("./Models/alspac_imputed_results.csv")


#What about a single plot before we go back to the 3a script to assemble
#these results alongside the complete case data?
d_dep_tab |>
  mutate(vars = str_remove(vars,"_f24")) |>
  # filter(vars %in% c("alt","cdcp1","hb","il6","insulin","neutrophils","wbc")) |>
  mutate(across(estimate:conf_high,as.double)) |>
  ggplot(aes(x = estimate, y = vars,color = model,xmin = conf_low,xmax = conf_high)) +
  geom_pointrange(position = position_dodge(width = 0.2)) +
  geom_vline(xintercept = 0,lty = 2)

# Experimental: Matching ------

# Single Examples

# #Do a quick example
# d <-  
#   midsobj2 |> 
#   complete(0) |> 
#   as_tibble() |> 
#   select(cisr_dep,il6_f24,bmi24,sex,smk24,audit_c) |>
#   drop_na()
# 
# d |>
#   lm(il6_f24 ~ cisr_dep * (bmi24 + sex + smk24 + audit_c),data = _) |>
#   avg_comparisons(variables = list(cisr_dep = 0:1)) 
# 
# #And do the same thing using matching
# d_match <- 
#   MatchIt::matchit(cisr_dep ~ bmi24 + sex + smk24 + audit_c,
#                    data = d,
#                    method = "full",
#                    estimand = "ATE")
# 
# #plot(d_match,type = "jitter",interactive = FALSE)
# #Black is treated, grey is control
# plot(d_match,type = "density",interactive = FALSE,
#      which.xs = ~sex + smk24 + audit_c + bmi24)
# 
# #So things get better after matching which is nice
# 
# #Now we fit our model
# d_match_fit <- 
#   d_match |>
#   match.data() |>
#   lm(il6_f24 ~ cisr_dep * (bmi24 + sex + smk24 + audit_c),
#      data = _,
#      weights = weights) 
# 
# 
# d_match_fit |>
#   avg_comparisons(variables = list(cisr_dep = 0:1),
#                   vcov = ~subclass,
#                   #newdata = subset(d_match |> match.data(),cisr_dep == 1), #ATT
#                   newdata = d_match |> match.data(),                      #ATE
#                   wts = "weights") 

#I think we will do matching, but perhaps as a sensitivity analysis, alongside
#the multiple imputation, not necessarily here

# ,
#          
#          matched = map(data,~MatchIt::matchit(cisr_dep ~ bmi24 + sex + smk24 + audit_c,
#                                       data = .x,
#                                       method = "full",
#                                       distance = "glm",
#                                       link = "logit",
#                                       estimand = "ATE")),
#          m4 = map(matched, ~.x |>
#                     match.data() |>
#                     lm(infl_value ~ cisr_dep * (bmi24 + sex + smk24 + audit_c),
#                        data = _,
#                        weights = weights)  |>
#                     avg_comparisons(variables = list(cisr_dep = 0:1),
#                                     vcov = ~subclass,
#                                     newdata = d_match |> match.data(),                      #ATE
#                                     wts = "weights") |>
#                     as_tibble()))
# 
# 
# 
# #Now we do the same thing but with the imputed dataset
# mice.models.a <- 
#   with(midsobj2, lm(formula = il6_f24 ~ cisr_dep * (bmi24 + sex + smk24 + audit_c)))
# 
# #Now we pool our imputed datasets to generate pooled estimates
# mice.results.a.tab = 
#   avg_comparisons(mice.models.a,variables = list(cisr_dep = 0:1)) |> 
#   tidy() |>
#   mutate(model = 'adjusted + bmi') |> 
#   select(-contrast) |> 
#   mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 8))) |>
#   as.matrix() |>
#   as_tibble() |>
#   mutate(across(estimate:conf.high,~as.double(.x)))  
# 
# #Combine with matching
# d_mi_match <- 
#   matchthem(cisr_dep ~ bmi24 + sex + smk24 + audit_c,
#             datasets = midsobj2,
#             method = "nearest",
#             approach = "within",
#             estimand = "ATT")
# 
# cobalt::love.plot(d_mi_match,stars = "std")
# 
# #So matching does help with our covariate balances, although the situation is
# #not perfect
# d_mi_results <- 
#   with(d_mi_match,
#        lm(il6_f24 ~ cisr_dep * (bmi24 + sex + smk24 + audit_c)))
# 
# 
# mice.results.m.tab = 
#   avg_comparisons(d_mi_results,
#                   variables = list(cisr_dep = 0:1),
#                   vcov = 'HC3') |> 
#   tidy() |>
#   mutate(model = 'matched inc bmi') |> 
#   select(-contrast) |> 
#   mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 8))) |>
#   as.matrix() |>
#   as_tibble() |>
#   mutate(across(estimate:conf.high,~as.double(.x)))  



# Experimental: Infl as DV -----

#This is not the approach we go on to use in the paper, but it is used in the 
#Symptom models, so here it is for reference

#So now we would ideally use a function to package all this up and then apply
#it to a list of our F24 variables. 

#We will fit the same models to our imputed data as we did for the complete-case
#data

#Make a function that takes an individual variable and returns the fit model
fit_mi_model <- function(data,dv,iv,covars){
  
  
  mice.models.u <- 
    with(data, glm(formula = paste0(dv, ' ~ ',iv) |>    as.formula(), 
                   family = binomial(link = 'logit')))
  
  mice.models.a <- 
    with(data, glm(formula = paste0(dv, ' ~ ',iv, ' * (', paste0(covars[!str_detect(covars,'bmi24')],collapse = ' + '),')') |> as.formula(), 
                   family = binomial(link = 'logit')))
  
  mice.models.f <- 
    with(data, glm(formula = paste0(dv, ' ~ ',iv, ' * (', paste0(covars,collapse = ' + '),')') |> as.formula(), 
                   family = binomial(link = 'logit')))
  
  
  #Now we pool our imputed datasets to generate pooled estimates
  mice.results.u.tab  <-
    marginaleffects::avg_comparisons(mice.models.u,variables = iv) |> 
    tidy() |>
    mutate(model = 'M1: unadjusted') |> 
    select(-contrast) |> 
    mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 8))) |>
    
    #This does seem insane, but for some reason the tibble that gets output by
    #this function is coming out as 300 MB, whereas it should be 3 kb. I have no
    #idea why either
    as.matrix() |>
    as_tibble() |>
    mutate(across(estimate:conf.high,~as.double(.x))) 
  
  
  mice.results.a.tab = 
    marginaleffects::avg_comparisons(mice.models.a,variables = iv) |> 
    tidy() |>
    mutate(model = 'M2: adjusted no bmi') |> 
    select(-contrast) |> 
    mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 8))) |>
    as.matrix() |>
    as_tibble() |>
    mutate(across(estimate:conf.high,~as.double(.x))) 
  
  
  mice.results.f.tab = 
    marginaleffects::avg_comparisons(mice.models.f,variables = iv) |> 
    tidy() |>
    mutate(model = 'M3: adjusted + bmi') |> 
    select(-contrast) |> 
    mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 8))) |>
    as.matrix() |>
    as_tibble() |>
    mutate(across(estimate:conf.high,~as.double(.x))) 
  
  
  
  #Do weighting
  weighted.datasets <- weightthem(paste0(iv , ' ~ ', paste0(covars,collapse = ' + ')) |> as.formula(),
                                  data = data, 
                                  approach = 'within',
                                  estimand = 'ATE',
                                  method = "optweight")
  
  #Inspect the results of the matching process
  # cobalt::love.plot(weighted.datasets)
  
  #Now fit the main model to the weighted, multiply imputed datasets
  weighted.models <- with(weighted.datasets,
                          glm(formula = paste0(dv, ' ~ ',iv, ' * (', paste0(covars,collapse = ' + '),')') |> as.formula(),
                              family = stats::quasibinomial(link = "probit")))
  
  #We use the quasibinomial distribution because the weights cause errors from fit.glm
  
  #Now we pool our imputed datasets to generate pooled estimates
  mice.results.w.tab  <- 
    weighted.models |> 
    marginaleffects::avg_comparisons(variables = iv,
                                     vcov = 'HC3') |> 
    tidy()  |>
    mutate(model = 'weighted') |> 
    select(-contrast) |> 
    mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 8))) |>
    as.matrix() |>
    as_tibble() |>
    mutate(across(estimate:conf.high,~as.double(.x))) 
  
  #Get all outcomes together to return
  
  return(bind_rows( mice.results.u.tab,
                    mice.results.a.tab,
                    mice.results.f.tab,
                    mice.results.w.tab))
  
}

#Now the problem with the energy balancing method, despite it being good - the 
#results make more sense than the weird ebal results and appear to be more well
#balanced than other methods - is that its really really slow. I think we might 
#have to use a different method for rolling this out to all variables.What about
#using optweight? It is faster
d_mi_results <- 
  tibble(vars =
           midsobj2$data |> 
           colnames() ) |> 
  filter((vars %nin% covars) & (vars %nin% dv) & !str_detect(vars,"id")) |>
  mutate(im_model = vector("list", 1)) |>
  filter(str_detect(vars,'_f24'))

#We have to do this in a whacky way becuase svyglm/mice/with seem to get very
#unhappy if you try to build a model formula using variables inside a function,
#they seem to only be able to take variables from the global environment and I 
#dont have the r skills to even formulate how to fix this, apart from changing 
#the value of iv in the global environment in each loop iteration. Sorry. 

for(i in 1:dim(d_mi_results)[1]){
  
  
  print(i)
  
  iv = d_mi_results$vars[[i]]
  d_mi_results$im_model[[i]] <- fit_mi_model(midsobj2,
                                             dv = 'cisr_dep', 
                                             iv = d_mi_results$vars[[i]],
                                             covars = covars) 
}




d_mi_results |>
  unnest(im_model) |>
  select(-term) |>
  janitor::clean_names() |>
  mutate(across(estimate:conf_high,as.double)) |>
  filter(p_value < 0.05 & model == "adjusted") |>
  #mutate(model = factor(model , levels = c("unadjusted","adjusted","weighted"))) |>
  mutate(vars = fct_reorder(vars,estimate,max)) |>
  ggplot(aes(y = vars,x = estimate,xmin = conf_low,xmax = conf_high,
             colour = model)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0,lty =2) +
  theme(panel.grid = element_blank())


#Lets save the mi model results for combining with the straightforward cross-
#sectional data
write_csv(d_mi_results |>
            unnest(im_model) |>
            select(-term) |>
            janitor::clean_names() |>
            mutate(across(estimate:conf_high,as.double)) ,
          paste(data_dir,"alspac_split_imputed_results_rev.csv",sep = '/'))


d_mi_results <- 
  read_csv(paste(data_dir,"alspac_split_imputed_results_rev.csv",sep = '/'))

#Uh oh
top_vars = 
  d_mi_results |>
  filter(model == "adjusted") |>
  mutate(p_adj = p.adjust(p_value,method = 'fdr')) |> 
  arrange(p_adj) |>
  filter(p_value < 0.05) |>
  pull(vars)


