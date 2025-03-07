
# Introduction ------

#load packages
pacman::p_load(tidyverse,
               survey,
               mice,
               broom,
               marginaleffects,
               ggdist,
               brglm2,
               MatchThem)

source('./Final/0-Common_data.R')

#The F9 variables have a different set of covariates to adjust for compared to the F24 data
# For F9 - bmi9, ses and sex
# For F24 - sex + bmi24 + audit_c + smk24

covars = c("sex","bmi24","smk24","audit_c","ph")

dv = 'cisr_dep'

# Use the imputated datasets to explore the cross-sectional associations at F24
# imputed.datasets.rf <- read_rds(paste(data_dir,"alspac_data_imputed.rds",sep = '/'))
imputed.datasets.rf <- read_rds("./Models/alspac_data_imputed.rds")

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



# Depression as DV -----

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
    glm(form_obj_1,
        family = binomial(link = "logit"),
        method = "brglmFit")
  })
  
  
  form_obj_3 <- 
    paste0(dv, ' ~ ',iv,' * (', paste0(covars,collapse = ' + '),')') |> 
    as.formula()
  
  mice.models.a <- 
    with(data, {
      # set form_obj's environment to the current one
      environment(form_obj_3) = environment()
      glm(form_obj_3,
          family = binomial(link = "logit"),
          method = "brglmFit")
      
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
  iv_i = d_dep_avg_results$vars[[i]]
  
  print(paste0(i,' : ', iv_i))
  
  d_dep_avg_results$im_model[[i]] <- fit_dep_model(data   = midsobj2,
                                                   dv     = 'cisr_dep',
                                                   iv     = iv_i,
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
