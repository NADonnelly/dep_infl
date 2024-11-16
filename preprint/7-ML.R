
# Introduction -----

#We are asking two questions with machine learning:

# - Can we predict ICD-10 depression status in the whole population, using a 
# combination of immuno-metabolic variables and demograhpic/past MH data

# - Within cases of depression, can we differential between profiles of symptoms 
# in individuals with depression using the same sets of variables

# Loading -----

pacman::p_load(tidyverse,
               tidymodels,
               themis,
               embed,
               finetune,
               bestNormalize,
               doFuture,
               tidyposterior,
               patchwork,
               crayon,
               probably,
               UpSetR,
               vip)

tidymodels_prefer()

conflicted::conflict_prefer('select','dplyr')
conflicted::conflict_prefer('filter','dplyr')
conflicted::conflict_prefer('map'   ,'purrr')
conflicted::conflicts_prefer(crayon::`%+%`)

#Set up some preferences for the crayon package so we can print progress 
#information to the command line as our ML process can take a few hours and its 
#helpful to have updates on progress

fancy  <- combine_styles(make_style("ivory"), 
                         make_style("grey20", bg = TRUE))

#Load our standard variables
source('./Final/0-Common_data.R')

#Load up our split training/test data
d_split <-
   read_rds(paste(data_dir,"alspac_data_split.rds",sep = '//'))


#Set our covariates
covars = c("sex","bmi24","smk24","audit_c")
 

#Do some clearing up
rm(list = c("wdn","valid_vars_f24","valid_vars_f9","v_0"))

#And get the training data for analysis
d_train <- 
 d_split |>
 training()

# Symptom scales ------

#At this point we need to simplify how we break up depression if we are going 
#investigate subgroups

#We can do something similar to covid paper - 
#https://www.nature.com/articles/s41590-024-01778-0#Abs1

#We define a  set of outcomes to try and predict. These can be based on the 
#symptoms we find most associated with immuno-metabolic variables in the PLS, 
#which are also similar to the existing literature e.g. in the Milaneschi
#paper

#-Combination 1: FTG, SLP, and SOM (aches and pains)
#-Combination 2: ANX and WOR
#-Combination 3: DEP and DID

d_symptom_score <- 
  d_train |>
  select(id,cisr_dep,ftg,slp,som,anx,wor,dep,did) |>
  mutate(across(.cols = ftg:did,.fns = ~as.numeric(.x)-1)) |>
  transmute(id,cisr_dep,
            c1 = ftg+slp+som,
            c2 = anx+wor,
            c3 = dep+did)

#Look at how these synthetic scores are distributed
p_cisr_scores <- 
  d_symptom_score |>
  pivot_longer(-c(cisr_dep,id)) |>
  drop_na(cisr_dep) |>
  mutate(name = case_when(name == "c1" ~ "Somatic \nSymptoms",
                          name == "c2" ~ "Anxiety \nSymptoms",
                          name == "c3" ~ "Depressive \nSymptoms")) |>
  mutate(cisr_dep = case_when(cisr_dep == 0 ~ "No Depression",
                              cisr_dep == 1 ~ "ICD10 Depression")) |>
  ggplot(aes(value,fill = cisr_dep)) +
  geom_histogram(binwidth = 1,colour = "black",linewidth = 0.25) +
  facet_wrap(~cisr_dep+name,ncol = 3, scales = "free") +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  scale_fill_manual(values = c("#32324B","#E1AF64"))+
  theme(panel.grid = element_blank(),
        axis.title  = element_text(size = 8),
        axis.text   = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(x = "Score", 
       y = "Count") +
  geom_vline(xintercept = c(5.5),lty = 2,linewidth = 0.25)

p_cisr_scores

ggsave('./Figures/depression_subtype_scales.pdf',plot = p_cisr_scores,width = 6,height = 4)


#So we would be looking at the anxiety and somatic score domains within people with depression - 
d_symptom_score |>
  filter(cisr_dep == 1) |>
  drop_na()

#Counts for depression

d_symptom_score |>
  filter(cisr_dep == 1) |>
  pivot_longer(-c(id,cisr_dep)) |>
  group_by(name) |>
  mutate(value = if_else(value > 0,1,0)) |>
  count(value)


#Do an upset plot to see overlap
p_overlap <- 
  d_symptom_score|>
  filter(cisr_dep == 1) |>
  drop_na() |>
  mutate(dep_split      = if_else(c3 >= 6,1,0),
         som_split      = if_else(c1 >= 5,1,0),
         anxiety_split  = if_else(c2 >= 6,1,0))|>
  select(-c(id,cisr_dep,c1:c3)) |>
  as.data.frame() |>
  upset(order.by = "freq")

p_overlap



#What is the relationship between our symptom scores and plain old depression sev
d_sev <- 
  d_train |>
  select(id,cisr_dep,cisr_dep_mod,cisr_dep_sev) |>
  mutate(across(.cols = cisr_dep,.fns = ~as.numeric(.x)-1)) |>
  transmute(id,ds = cisr_dep+cisr_dep_mod+cisr_dep_sev) |>
  mutate(dep_sev = case_when(ds == 0 ~ "no_dep",
                             ds == 1 ~ "mild_dep",
                             ds == 2 ~ "mod_dep",
                             ds == 3 ~ "sev_dep")) |>
  mutate(dep_sev = factor(dep_sev, levels = c("no_dep","mild_dep","mod_dep","sev_dep")))

#How do these guys add up?
d_sev |>
  drop_na(ds) |>
  count(dep_sev) |>
  ggplot(aes(x = dep_sev,y = n)) +
  geom_col()

#Compare our scores to depression severity
d_symptom_score |>
  left_join(d_sev,by = "id") |>
  pivot_longer(c1:c3) |>
  mutate(name = case_when(name == "c1" ~ "Somatic \nSymptoms",
                          name == "c2" ~ "Anxiety \nSymptoms",
                          name == "c3" ~ "Depressive \nSymptoms")) |>
  drop_na() |>
  ggplot(aes(x = dep_sev,y = value)) +
  geom_boxplot() +
  facet_wrap(~name,ncol = 3)

#So all three subscales have higher scores in more severe depression, which you would expect.
#Its not simply a case that one is more related to severity than the others


# Prepare Data for ML ------

## Make the new outcomes -----

d_train_outcomes <- 
  d_split |>
  training() |>
  select(id,cisr_dep,cisr_dep_mod,cisr_dep_sev,ftg,slp,som,anx,wor,dep,did) |>
  mutate(across(.cols = c(cisr_dep,ftg:did),.fns = ~as.numeric(.x)-1)) |>
  transmute(id,cisr_dep,
            c1 = ftg+slp+som,
            c2 = anx+wor,
            c3 = dep+did) |>
  drop_na() |>
  mutate(som_plus   = if_else((c1 >= 6) & (cisr_dep == 1),1,0),
         anx_plus   = if_else((c2 >= 6) & (cisr_dep == 1),1,0),
         dep_plus   = if_else((c3 >= 6) & (cisr_dep == 1),1,0)) |>
  select(-c(c1:c3)) |>
  mutate(across(.cols = contains("_split")|contains("_any")|contains("_plus"),.fns = ~factor(.x)),
         dep_class = factor(dep_class))



d_test_outcomes <- 
  d_split |>
  testing() |>
  select(id,cisr_dep,cisr_dep_mod,cisr_dep_sev,ftg,slp,som,anx,wor,dep,did) |>
  mutate(across(.cols = c(cisr_dep,ftg:did),.fns = ~as.numeric(.x)-1)) |>
  transmute(id,cisr_dep,
            c1 = ftg+slp+som,
            c2 = anx+wor,
            c3 = dep+did) |>
  drop_na() |>
  mutate(som_plus   = if_else((c1 >= 6) & (cisr_dep == 1),1,0),
         anx_plus   = if_else((c2 >= 6) & (cisr_dep == 1),1,0),
         dep_plus   = if_else((c3 >= 6) & (cisr_dep == 1),1,0)) |>
  select(-c(c1:c3)) |>
  mutate(across(.cols = contains("_split")|contains("_any")|contains("_plus"),.fns = ~factor(.x)),
         dep_class = factor(dep_class))

#How common are our outcomes in the test data? This will be a potential issue


#Diagnosis - related

#CISR depression - 264 cases in the full training data, 73 in the test data
d_train_outcomes |> count(cisr_dep)
d_test_outcomes  |> count(cisr_dep)


# Depressive symptoms

# Any depressive symptoms - 166 cases, 485 controls 
d_train_outcomes  |> filter(cisr_dep == 1) |> count(dep_plus)
d_test_outcomes  |> filter(cisr_dep == 1) |> count(dep_plus)


#Somatic symptom - related

# Depression and somatic symptoms - 29 = 0, 44 = 1
d_train_outcomes  |> filter(cisr_dep == 1) |> count(som_plus)
d_test_outcomes  |> filter(cisr_dep == 1) |> count(som_plus)


# Anxiety - related

#Anxiety symptoms + depression diagnosis - 30 cases, 621 controls
d_train_outcomes   |> filter(cisr_dep == 1)|> count(anx_plus)
d_test_outcomes   |> filter(cisr_dep == 1)|> count(anx_plus)


## Do the standard Data Prep ------

#Load the cluster data 
infl_clus <- 
  read_csv('./Models/f24_infl_variable_clusters.csv')

#Prepare the training data
d_train <- 
  d_split |>
  training() |>
  
  #Remove some variables we do not want including in any model
  select(-c(cisr_dep_mod,cisr_dep_sev,eth, 
            cisr_dep:pan,
            atypical:cape_total,
            starts_with("ad_"), starts_with("smfq_c"))) |>
  
  #Convert sex and smoking to being dummy coded
  mutate(sex = if_else(sex == "Male",0,1)) |>
  mutate(across(.cols = c(sex,smk24),.fns = as.integer)) |>
  
  #Some other manipulations. We are going to take a limited subset of clinical
  #relevant variables. Hopefully although this is inevitably subjective, we can 
  #make sensible choices
  relocate(c(sex,bmi24,smk24,audit_c,sc,me,md),.after = id) |>
  
  #Now select the final sets of variables for our analysis
  select(
    
    #ID, outcomes, covariates
    id:md,
    
    #Our F24 bloods
    contains("_f24"),
    
    
    #SMFQ Data
    
    #The previous year's SMFQ data on the basis that you could potentially
    #ask how the person was the previous year?
    smfq_10,smfq_12.5, smfq_13.5,smfq_16.5, smfq_17.5,smfq_23,
    
    #Early life difficulties
    sdq7,sdq9,
    
    #Age 17 CISR depression but we will turn it into a binary
    cisr_mild_17,cisr_mod_17,cisr_sev_17,
    
    #Trauma/ACEs (although the aces variable is nearly 50% missing)
    ace_total_classic,trauma_total,
    
    #Keep this for stratification later
    multi_demo) |>
  
  mutate(dep_17 = if_else((cisr_mild_17 + cisr_mod_17 + cisr_sev_17) > 0,1,0 )) |>
  mutate(dep_17 = as.integer(dep_17))|>
  select(-c(cisr_mild_17,cisr_mod_17,cisr_sev_17)) |>
  
  #Now you need to graft on your outcomes
  left_join(d_train_outcomes,by = "id") |>
  relocate(cisr_dep:dep_sev,.after = id)


## Extreme values ------

#We need to reproduce our extreme value analysis for cluster 3
#Make a set of normative values for the mu and sd for each blood measure
d_norm_val <- 
  d_train |>
  select(id,cisr_dep,contains('_f24')) |>
  drop_na(cisr_dep) |>
  pivot_longer(-c(id,cisr_dep)) |>
  nest_by(name,cisr_dep) |>
  pivot_wider(names_from = cisr_dep,names_prefix = "dep_",values_from = data) |>
  mutate(dep_0_mu = map_dbl(dep_0, ~mean(.x$value,na.rm = T)),
         dep_0_sd = map_dbl(dep_0, ~sd(.x$value,na.rm = T))) |>
  select(-c(dep_0,dep_1)) |>
  
  #Filter for our cluster 3 variables - 
  filter(name %in% (infl_clus |> mutate(infl = paste0(infl,'_f24')) |> filter(clus == 3) |> pull(infl)))

#Now apply to the training data
d_train_ev <- 
  d_train |>
  select(id,cisr_dep,contains('_f24')) |>
  drop_na(cisr_dep) |>
  pivot_longer(-c(id,cisr_dep)) |>
  right_join(d_norm_val,by = "name") |>
  mutate(norm_value = (value - dep_0_mu)/dep_0_sd) |>
  mutate(ex_dev = if_else(abs(norm_value) > 2.6,1,0)) |>
  select(id,cisr_dep,name,ex_dev) |>
  group_by(id) |>
  reframe(cluster_3_ev = sum(ex_dev,na.rm = T)) |>
  arrange(id) 

d_train <- 
  d_train |>
  left_join(d_train_ev,by = "id")


## Apply to test data ------

d_test <- 
  d_split |>
  testing() |>
  
  #Remove some variables we do not want including in any model
  select(-c(cisr_dep_mod,cisr_dep_sev,eth, 
            cisr_dep:pan,
            atypical:cape_total,
            starts_with("ad_"), starts_with("smfq_c"))) |>
  
  #Dummy coding
  mutate(sex = if_else(sex == "Male",0,1)) |>
  mutate(across(.cols = c(sex,smk24),.fns = as.integer)) |>
  relocate(c(sex,bmi24,smk24,audit_c,sc,me,md),.after = id) |>
  
  mutate(dep_17 = if_else((cisr_mild_17 + cisr_mod_17 + cisr_sev_17) > 0,1,0 )) |>
  mutate(dep_17 = as.integer(dep_17)) |>
  
  left_join(d_test_outcomes,by = "id") |>
  relocate(cisr_dep:dep_sev,.after = id)



#Get the extreme value data
d_test_ev <- 
  d_test |>
  select(id,cisr_dep,contains('_f24')) |>
  drop_na(cisr_dep) |>
  pivot_longer(-c(id,cisr_dep)) |>
  
  #Use the normative values from the training data
  right_join(d_norm_val,by = "name") |>
  mutate(norm_value = (value - dep_0_mu)/dep_0_sd) |>
  mutate(ex_dev = if_else(abs(norm_value) > 2.6,1,0)) |>
  select(id,cisr_dep,name,ex_dev) |>
  group_by(id) |>
  reframe(cluster_3_ev = sum(ex_dev,na.rm = T)) |>
  arrange(id) 

d_test <- 
  d_test |>
  left_join(d_test_ev,by = "id")


d_test <- 
  d_test |>
  
  #Remove same variables
  select(c(all_of(d_train |> colnames()))) |>
  drop_na(cisr_dep)


#Prove we have the same columns in both datasets
setdiff(d_train |> colnames(),d_test |> colnames())

#Now, the size of our final training dataset is
dim(d_train)



## Variables and recipes -----

#Get variable types
d_var_type = 
  d_train |>
  
  #Remove our categorical variables
  select(-c(id:dep_sev,sc, me,sex,smk24,multi_demo)) |>
  pivot_longer(everything()) |>
  drop_na() |>
  group_by(name) |>
  mutate(value = factor(value)) |>
  summarise(nl = length(unique(value))) 

#Lets say that above 30 you count as a continuous variable
cont_vars = 
  d_var_type |>
  filter(nl > 30) |>
  pull(name)

#One exception is the basophils, which are normalish 
#distributed but appear to be discretised somehow. We add this to the continuous
#variable list
cont_vars <- c(cont_vars,"basophils_f24")


#Set 1: Bloods + Covariates

#Take the cluster 1 and cluster 3 variables for consistency with other 
#aspects of the paper

f24_blood_names <- 
  infl_clus |> 
  filter(clus %in% c(1,3)) |> 
  mutate(infl = paste0(infl,"_f24")) |>
  pull(infl)

iv_f24_blood <- 
  c("sex","audit_c","smk24","bmi24",
    f24_blood_names,
    "cluster_3_ev")

#Note we include our covariates in the model as we can see that IM-variables are
#pretty good at classifying on our covariates e.g. you can get sex from the Hb
#so we want the model the get information about depression from the IM variables
#and it can get sex data from the sex variable

#Update the continuous variables to include only those variables we are 
#including
cont_vars <- intersect(iv_f24_blood,cont_vars)

#Set 2: Social Variables

#Make a set of social variables we know are likely to be related to
#depression but are not blood-related, to act as a positive control
iv_soc <- 
  c("sex","audit_c","smk24","bmi24", # Demographic variables
    "sc",                            # Maternal Social class- correlates 
                                     # substantially with maternal education so 
                                     # we will only pick one
    "ace_total_classic",             # ACES
    "trauma_total"                   # Trauma
  )


#Set 3: Past MH variables

#Make a set of variables related to personal history of depression as another 
#positive control

iv_dep <- 
  c("md",   # Maternal EPDS
    d_train |> select(starts_with("smfq")) |> colnames(),  #Earlier smfq data
    d_train |> select(starts_with("sdq")) |> colnames(),   #Earlier sdq data - use age 7 and 9
    "dep_17")

#Set 4: Full dataset

#Make a full set of variables combining the blood and clinical variables
iv_full <- 
  c(iv_soc,iv_dep,iv_f24_blood[str_detect(iv_f24_blood,"_f24")],"cluster_3_ev")

#Make recipes 
f_blood_24  <- as.formula(paste0('class ~ ', paste0(iv_f24_blood , collapse = ' + '), collapse = ""))
f_soc       <- as.formula(paste0('class ~ ', paste0(iv_soc       , collapse = ' + '), collapse = ""))
f_dep       <- as.formula(paste0('class ~ ', paste0(iv_dep       , collapse = ' + '), collapse = ""))
f_full      <- as.formula(paste0('class ~ ', paste0(iv_full      , collapse = ' + '), collapse = ""))


## Models -----

#Model specs

#Elastic net
en_mod <- 
  logistic_reg(mode = "classification",
               penalty = tune(), 
               mixture = tune()) %>% 
  set_engine("glmnet")


## Splits -----

#Set a random number seed for reproducibility
set.seed(1)

#Set our data splitting ratios and number of folds
split_prop = 4/5

n_outer    = 20
n_inner    = 20
n_fold     = 5

#Prepare for (nested) cross validation
d_folds = nested_cv(d_train, 
                    outside = vfold_cv(v = n_fold, repeats = n_outer/n_fold, strata = multi_demo), 
                    inside  = bootstraps(times = n_inner, strata = multi_demo)) 


# Double Loooooop -------

#Fit our models using nested cross validation, nesting over each symptom
dv_list <- c("cisr_dep",
             "som_plus", 
             "anx_plus", 
             "dep_plus")

#Make output table
ml_result = 
  tibble(symptom = dv_list,
         ml      = vector(mode = 'list',length = 1))


for(j in 1:nrow(ml_result)){
  
  dv = ml_result$symptom[[j]]
  
  cat(magenta(
    'This is outer loop ' %+%
      blue$underline$bold(j) %+%
      ' of ' %+%
      blue$underline$bold(nrow(ml_result)) %+% 
      '; Current symptom = ' %+% blue$underline$bold(dv) %+% '\n'
  ))
  
  #Preallocate an output tibble
  d_fold_results = 
    d_folds |>
    dplyr::select(id,id2) |>
    mutate(.results       = vector(mode = "list",length = nrow(d_folds)),
           .mod_posterior = vector(mode = "list",length = nrow(d_folds)))
  
  
  
  #Loop through each outer fold
  for(i in 1:nrow(d_folds)){
    
    #Tell the user where we are in our looping
    cat(green(
      'This is inner loop ' %+%
        blue$underline$bold(i) %+%
        ' of ' %+%
        blue$underline$bold(nrow(d_folds)) %+% '\n'
    ))
    
    
    #Get the outer data for this iteration of the loop
    d_outer <-  
      d_folds$splits[[i]]
    
    
    #Make pre-processing recipes to prepare the training data for this loop
    cat(fancy("Pre-processing..."), "\n")
    
    
    #We do some prep to rename the dependent variable and make sure its a factor
    if(str_detect(dv,'_plus')){
      
      d_outer_train <- 
        d_outer |>
        analysis() |>
        filter(cisr_dep == 1) |>
        dplyr::select(-c(id,multi_demo,cisr_dep)) |>
        rename(class = !!dv) |>
        mutate(class = factor(class)) |>
        dplyr::select(class, all_of(iv_full)) |> 
        drop_na(class)
      
      d_inner <- 
        d_outer_train |>
        bootstraps(times = n_inner)
      
      d_outer_test <- 
        d_outer |>
        assessment()|>
        filter(cisr_dep == 1) |>
        dplyr::select(-c(id,multi_demo,cisr_dep)) |>
        rename(class = !!dv) |>
        mutate(class = factor(class)) |>
        dplyr::select(class, all_of(iv_full))|> 
        drop_na(class)
      
      
    }else{
      
      d_outer_train <- 
        d_outer |>
        analysis() |>
        dplyr::select(-c(id,multi_demo)) |>
        rename(class = !!dv) |>
        mutate(class = factor(class)) |>
        dplyr::select(class, all_of(iv_full)) |> 
        drop_na(class)
      
      d_inner <- 
        d_outer_train |>
        bootstraps(times = n_inner)
      
      d_outer_test <- 
        d_outer |>
        assessment()|>
        dplyr::select(-c(id,multi_demo)) |>
        rename(class = !!dv) |>
        mutate(class = factor(class)) |>
        dplyr::select(class, all_of(iv_full))|> 
        drop_na(class)
      
    }
    
    #We now rebuild a split object with the prepared data
    d_outer_split =
      make_splits(d_outer_train,d_outer_test)
    
    
    
    #Recipe with blood data
    rec_blood_24 =
      d_outer_train |>
      recipe(formula = f_blood_24) |>
      
      #Filter variables with near-zero variance (should do nothing)
      step_nzv(all_predictors()) |>
      
      #Filter variables with >50% missingness (again should do nothing)
      step_filter_missing(all_predictors(), threshold = 0.5) |>
      
      #Scale continuous variables to mean zero and SD 1
      step_orderNorm(!!cont_vars) |>
      
      #Impute missing predictors using KNN
      step_impute_knn(all_predictors(),
                      impute_with = imp_vars(all_predictors()),
                      neighbors   = 5) |>
      
      #Remove any rows that retain NAs (should not be in issue after imputation
      #but we include for quality assurance)
      step_naomit(all_predictors(),
                  skip = TRUE) |>
      
      #Use adasyn to generate new synthetic data to account for class imbalance
      step_adasyn(class,
                  over_ratio = 1,
                  neighbors  = 5,
                  seed       = 1,
                  skip       = TRUE)
    
    
    #The sociodemographic data
    rec_soc =
      d_outer_train |>
      recipe(formula = f_soc) |>
      step_nzv(all_predictors()) |>
      step_filter_missing(all_predictors(), threshold = 0.5) |>
      step_lencode_glm(sc, outcome = vars(class)) |>
      step_impute_knn(all_predictors(),
                      impute_with = imp_vars(all_predictors()),
                      neighbors   = 5) |>
      step_naomit(all_predictors(),
                  skip = TRUE) |>
      step_adasyn(class,
                  over_ratio = 1,
                  neighbors  = 5,
                  seed       = 1,
                  skip       = TRUE)
    
    
    #Past MH Data
    rec_dep =
      d_outer_train |>
      recipe(formula = f_dep) |>
      step_nzv(all_predictors()) |>
      step_filter_missing(all_predictors(), threshold = 0.5) |>
      step_impute_knn(all_predictors(),
                      impute_with = imp_vars(all_predictors()),
                      neighbors   = 7) |>
      step_naomit(all_predictors(),
                  skip = TRUE) |>
      step_adasyn(class,
                  over_ratio = 1,
                  neighbors  = 5,
                  seed       = 1,
                  skip       = TRUE)
    
    
    
    #All variables together in a big soup
    rec_full =
      d_outer_train |>
      recipe(formula = f_full) |>
      step_nzv(all_predictors()) |>
      step_filter_missing(all_predictors(), threshold = 0.5) |>
      step_orderNorm(!!cont_vars) |>
      step_lencode_glm(sc, outcome = vars(class)) |>
      step_impute_knn(all_predictors(),
                      impute_with = imp_vars(all_predictors()),
                      neighbors   = 5)|>
      step_naomit(all_predictors(),
                  skip = TRUE) |>
      step_adasyn(class,
                  over_ratio = 1,
                  neighbors  = 5,
                  seed       = 1,
                  skip       = TRUE)
    
    
    #Now we can move on to model fitting
    cat(fancy("Fitting Models..."), "\n")
    
    
    #Make workflow sets so we can run all out models together using the workflowsets package
    wf_par <- 
      workflow_set(
        preproc = list(
          
          blood = rec_blood_24,
          soc   = rec_soc,
          dep   = rec_dep,
          full  = rec_full

        ),
        models  = list(
          en    = en_mod))
    
    
    #Prepare parameters
    for(k in 1:dim(wf_par)[1]){
      
      t_params <- extract_parameter_set_dials(wf_par,id = wf_par$wflow_id[k])
      
      t_params <- 
        t_params |>
        finalize(d_train)
      
      wf_par <- 
        wf_par |>
        option_add(param_info = t_params , id = wf_par$wflow_id[k])  
      
    }
    
    
    #Set a timer
    tictoc::tic()
    
    #Set up controls for the finetune anova race method
    race_ctrl <-
      control_race(
        parallel_over = "everything",
        save_pred     = TRUE,
        save_workflow = TRUE
      )
    
    
    #Set up to run these in parallel (I have a 16 core CPU so we use all 16 cores)
    registerDoFuture()
    plan(multisession, workers = 16)
    
    #Fit our workflows that work in parallel
    wf_results <-
      wf_par |>
      workflow_map(
        "tune_race_anova",
        seed      = 1,
        resamples = d_inner,
        control   = race_ctrl,
        grid      = 30,
        metrics   = metric_set(bal_accuracy,mcc)
      )
    
    
    #Turn off parallel mode
    plan(sequential)  
    
    
    #Only keep models that ran properly
    wf_results <- 
      wf_results |>
      mutate(is_ok = map_lgl(result,~is_tibble(.x))) |>
      filter(is_ok) |>
      select(-is_ok)
    
    tictoc::toc()
    
    # rankings <-
    #   wf_results |>
    #   mutate(is_ok = map_lgl(result,~is_tibble(.x))) |>
    #   filter(is_ok) |>
    #   select(-is_ok) |>
    #   rank_results(select_best = TRUE,rank_metric = "bal_accuracy")

    
    cat(fancy("Calculating Model Performance..."), "\n")

    

    #Evaluate models on the outer fold data
    tictoc::tic()
    
    #Get the settings for the  best version of every model, fit to
    #the outer fold training data, and then test on the outer fold test
    #data
    best_results_wf = 
      tibble(wflow_id = wf_results$wflow_id) |>

      #Here we extract the workflows and the best parameters in each workflow
      mutate(res         = purrr::map(wflow_id,~workflowsets::extract_workflow_set_result(wf_results,id = .x)))  |>
      mutate(best_params = purrr::map(res, tune::select_best,metric = "bal_accuracy")) |>
      mutate(best_wf     = purrr::map(wflow_id, ~workflowsets::extract_workflow(wf_results,id = .x))) 
    
    
    best_results_wf = 
      best_results_wf |>
      
      #We finalize the workflow by combining the best parameters with the workflow
      mutate(best_wf_final  = map2(best_wf,best_params,~finalize_workflow(.x,.y))) |>
      
      #Fit the best models to the full training data for this outer fold and apply to the test 
      #data for this outer fold (not previously seen by any part of the fitting process)
      mutate(outer_full_fit = purrr::map(best_wf_final, ~last_fit(.x,
                                                                  split = d_outer_split,
                                                                  metrics = metric_set(roc_auc,mcc,brier_class,bal_accuracy)))) 
    
    best_results_wf = 
      best_results_wf  |>
      
      #An error checking step
      mutate(has_metrics = map_lgl(outer_full_fit,~if_else(dim(.x |> 
                                                                 dplyr::select(.metrics) |> 
                                                                 unnest(.metrics))[1] > 0,
                                                           TRUE,FALSE))) |> 
      filter(has_metrics) |>
      
      
      #Put those metrics somewhere easy to grab
      mutate(best_metric = map_dbl(outer_full_fit,~.x |> 
                                     dplyr::select(.metrics) |> 
                                     unnest(.metrics) |> 
                                     filter(.metric == "bal_accuracy") |> 
                                     pull(.estimate))) |>
      arrange(best_metric)
    
    
    #Store the results for this fold
    d_fold_results$.results[[i]] =
      best_results_wf |>
      
      #Can we try to save memory by reducing this down?
      dplyr::select(-c(res,best_wf,has_metrics)) |>
      
      unnest(outer_full_fit) |> 
      
      #Remove the split data as we don't need it any more
      dplyr::select(c(wflow_id,best_params,best_wf_final,.metrics,best_metric)) |>
      arrange(-best_metric) 
    
    
    tictoc::toc()
    
    
    # #Quick plot of the outer fold test performance
    # best_results_wf |>
    #   separate_wider_delim(cols = wflow_id,names = c("vars","mod_type"),delim = "_") |>
    #   mutate(vars = factor(vars,levels = c("blood","soc","dep","full"))) |>
    #   ggplot(aes(x = vars,y = best_metric,fill = mod_type)) +
    #   geom_col(position = position_dodge(width = 0.5))

    
    #lets do some cleaning up
    rm(list = c("d_inner", "d_outer","d_outer_test","d_outer_train",
                "rec_impute","rec_full","rec_blood_24","rec_soc","rec_dep",
                "wf_par","wf_results",
                "best_results_wf",
                "roc_mod","roc_mod_post"))  
    
  }
  
  
  ml_result$ml[[j]] = d_fold_results
  

  rm(list = c("d_fold_results"))
  
}

write_rds(ml_result,paste("./Models/depression_subtype_results.rds"))
#ml_result <- read_rds("./Models/depression_subtype_results.rds")

## Prep results -----

#Make a reduced set of results to make a less-memory intense object
ml_result_compact = 
  ml_result |>
  unnest(ml) |>
  #select(-.mod_posterior) |>
  unnest(.results) |>
  select(-c(best_wf_final)) |> 
  rename(id1 = id) |>
  select(c(symptom,id1,id2,wflow_id,best_params,.metrics))

lobstr::obj_size(ml_result_compact)


#Do some clearing up
rm(list = c('d_folds','logistic_reg_spec','pca_train',
            'cor_vars','d_test_id','cor_train','infl_clus',
            'iv_dep','iv_f24_blood','iv_full','iv_onlyBlood','iv_soc',
            'n_fold','n_inner','n_outer','j','k','i','mr','ms',
            'priority_bloods','split_prop','d_final_split'))

gc()

#We need to reduce the size of the ml_result object in order to get anything done
ml_result <- 
  ml_result |>
  select(symptom,ml) |>
  unnest(ml) |>
  mutate(id_loop = interaction(id,id2,sep = "_")) |>
  select(-c(id,id2)) |>
  unnest(.results)|>
  #select(-c(id,.notes,.predictions,splits,has_metrics)) |>
  nest_by(symptom,.key = 'ml')

gc()

lobstr::obj_size(ml_result)

#OK this does free up a decent amount of memory


## Basic Plotting -----

#Really quick and dirty plot


#Plot the balanced accuracy
ml_result_compact |>
  select(symptom:wflow_id,.metrics) |>
  unnest(.metrics) |>
  filter(.metric == "bal_accuracy") |>
  separate_wider_delim(wflow_id,delim = '_',names = c("vars","mod_type")) |>
  ggplot(aes(x = symptom,y = .estimate,shape = vars,colour = vars,
             group = vars)) +
  geom_point(position = position_dodge(width = 0.3)) +
  facet_wrap(~mod_type,ncol = 1) +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_text(angle = -90)) +
  geom_hline(yintercept = .50,lty = 2) +
  labs(x = "Model",y = "MCC")

## Performance model ------

#Which metric should we use for comparison purposes?
cv_results <- 
  ml_result_compact|>
  select(symptom:wflow_id,.metrics) |>
  unnest(.metrics) |>
  select(-c(.estimator)) |>
  drop_na(.estimate) |>
  mutate(id = interaction(id1,id2,sep = "_")) |>
  select(-c(id1,id2)) |>
  nest_by(symptom,.metric)

#We can directly compare with our own rstanarm model. The ~ 0 + formula means 
#the model is estimating the posterior distribution for the performance of each
#model individually without an intercept
cv_mod <- 
  cv_results |>
  mutate(model = list(
    rstanarm::stan_glmer(.estimate ~ 0 + wflow_id + (1|id),
                         data = data,
                         seed = 1, refresh = 0,
                         chains = 4, cores = 8,
                         iter = 10000, warmup = 2000,
                         prior = rstanarm::normal(0,1),
                         prior_intercept = rstanarm::normal(0,1.5),
                         prior_aux = rstanarm::exponential(rate = 1))
    
  ))


cv_mod <- 
  cv_mod |>
  mutate(s_ord = 
           list(
             model |> 
               tidybayes::tidy_draws() |> 
               select(.chain:.draw,starts_with('wflow_id')) |>
               pivot_longer(starts_with("wflow_id"))  |> 
               mutate(name = str_remove(name,"wflow_id")) |>
               group_by(name) |>
               tidybayes::mean_hdci(value) |>
               arrange(value)|>
               pull(name) 
           )
  )


## Plots ------

wflow_name <- 
  c( "Full"                  = "full",
     "Sociodemographic"      = "soc",
     "Mental Health History" = "dep",
     "Immuno-metabolic"      = "blood")


metric_name <- 
  c( "Balanced Accuracy" = "bal_accuracy",
     "Brier Score"       = "brier_class",
     "MCC"               = "mcc",
     "AUROC"             = "roc_auc")


#Sort out prep for plotting
cv_mod <- 
  cv_mod |>
  mutate(disp_line = case_when(.metric == "AUROC" ~ 0.5,
                               .metric == "Balanced Accuracy" ~ 0.5,
                               .metric == "Brier Score" ~ 0.25,
                               .metric == "MCC" ~ 0),
         disp_lims = case_when(.metric == "AUROC" ~ c(0.45,0.85) |> list(),
                               .metric == "Balanced Accuracy" ~ c(0.45,0.75) |> list(),
                               .metric == "Brier Score" ~ c(0,0.3) |> list(),
                               .metric == "MCC" ~ c(-0.05,0.4) |> list()))  |>
  mutate(.metric = fct_recode(.metric, !!!metric_name)) |>
  mutate(.metric = fct_relevel(.metric, c("Brier Score","Balanced Accuracy","AUROC","MCC"))) |>
  arrange(.metric) |>
  rowwise(.metric) 


d_cv_mod_plot <- 
  cv_mod |>
  # filter(.metric %in% c("MCC")) |>
  filter(.metric %in% c("Balanced Accuracy")) |>
  mutate(symptom_group = case_when(str_detect(symptom,"anx")   ~ "anxiety symptoms",
                                   str_detect(symptom,"dep")   & !str_detect(symptom,"cisr") ~ "depressive symptoms",
                                   str_detect(symptom,"som")   ~ "somatic symptoms")) |>
  mutate(symptom = factor(symptom, levels = c("cisr_dep",
                                              "anx_plus",
                                              "dep_plus",
                                              "som_plus"))) |>
  arrange(symptom) |>
  mutate(draws = 
           list(model |> 
                  tidybayes::tidy_draws()  |> 
                  select(.chain:.draw,starts_with('wflow_id')) |>
                  pivot_longer(starts_with("wflow_id"))  |> 
                  mutate(name = str_remove(name,"wflow_id")))) |>
  ungroup() |>
  select(symptom_group,symptom,draws) |>
  unnest(cols = draws) |>
  separate_wider_delim(name, delim = "_",names = c("var_set","ml_method")) |>
  mutate(var_set = fct_recode(var_set, !!!wflow_name)) |>
  mutate(var_set = fct_relevel(var_set, c('Immuno-metabolic','Sociodemographic','Mental Health History','Full')))


p_all_wflow <- 
  d_cv_mod_plot |>
  ggplot(aes(y = symptom_group, x = value,fill = var_set)) +
  ggdist::stat_halfeye()+
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey80"),
        legend.position = "none") +
  labs(x = "Balanced Accuracy", 
       y = "Variable Set") +
  # coord_cartesian(xlim = c(-0.1,0.5)) +
  # geom_vline(xintercept = c(0,0.1,0.2,0.3,0.4),lty = 2,linewidth = 0.5) +
  coord_cartesian(xlim = c(0.4,0.8)) +
  geom_vline(xintercept = c(0.5,0.6,0.7,0.8),lty = 2,linewidth = 0.5) +  
  scico::scale_fill_scico_d(palette = "oleron") + 
  facet_wrap(~var_set,nrow = 1)

p_all_wflow

#Lets boil this down into two plots which I think will illustrate the strengths of this model

#First - is the model better able to predict whether there are any depressive symptoms, compared 
#to anxiety or somatic symptoms?
d_cv_mod_plot |>
  filter(symptom %in% c("anx_plus","dep_plus","som_plus")) |>
  ggplot(aes(y = symptom, x = value,fill = var_set)) +
  ggdist::stat_halfeye(alpha = 0.4)+
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey80")) +
  labs(x = "Balanced Accuracy", 
       y = "Variable Set") +
  coord_cartesian(xlim = c(0.4,0.8)) +
  geom_vline(xintercept = c(0.5,0.6,0.7,0.8),lty = 2,linewidth = 0.5) +  
  scico::scale_fill_scico_d(palette = "oleron")

#Within the people with depression the IM model does better for the somatic symptoms 
#than depression and anxiety, MH history is much better for depressive symptoms


# Final performance + VI -------

## Bootstrap performance ------

#Fit the best model to the full training set and measure performance in the testing set

n_boot <- 1000
n_core <- 12

boot_out <- 
  tibble(symptom     = dv_list,
         test_result = vector(mode = 'list',length = 1),
         test_vi     = vector(mode = 'list',length = 1),
         boot_result = vector(mode = 'list',length = 1),
         boot_sens   = vector(mode = 'list',length = 1))

for(j in 1:nrow(ml_result)){
  
  #Get the DV for this loop
  dv = dv_list[[j]]
  
  
  if(str_detect(dv,'_plus')){
    
      f_train <- 
        
        #Prepare the full training data
        d_train |>
        dplyr::select(-c(multi_demo)) |>
        filter(cisr_dep == 1) |>
        rename(class = !!dv) |>
        drop_na(class) |>
        dplyr::select(-c(any_of(dv_list),ends_with("_any"),ends_with("_class"),ends_with("_sev")) )|>
        mutate(class = factor(class))
      
      f_test <- 
        
        #Prepare the testing data
        d_test |>
        dplyr::select(-c(multi_demo)) |>
        filter(cisr_dep == 1) |>
        rename(class = !!dv) |>
        drop_na(class) |>
        dplyr::select(-c(any_of(dv_list),ends_with("_any"),ends_with("_class"),ends_with("_sev")) )|>
        mutate(class = factor(class)) 
      
  } else{
    
    f_train <- 
      
      #Prepare the full training data
      d_train |>
      dplyr::select(-c(multi_demo)) |>
      rename(class = !!dv) |>
      drop_na(class) |>
      dplyr::select(-c(any_of(dv_list),ends_with("_any"),ends_with("_class"),ends_with("_sev")) )|>
      mutate(class = factor(class))
    
    f_test <- 
      
      #Prepare the testing data
      d_test |>
      
      dplyr::select(-c(multi_demo)) |>
      rename(class = !!dv) |>
      drop_na(class) |>
      dplyr::select(-c(any_of(dv_list),ends_with("_any"),ends_with("_class"),ends_with("_sev")) )|>
      mutate(class = factor(class)) 
    
  }

  ff <- 
    make_splits(f_train,f_test)
  
  #Get the best performing model
  final_performance <- 
    ml_result |>
    filter(symptom == dv) |>
    pluck('ml',1) |>
    group_by(wflow_id) |>
    slice_max(best_metric) |>
    mutate(data = list(ff)) |>
    mutate(final_fit = map2(best_wf_final,data,~last_fit(.x,split = .y,
                                                         metrics = metric_set(roc_auc,mcc,bal_accuracy))))
  
  final_performance <- 
    final_performance |>
    select(wflow_id,final_fit,best_params) |>
    unnest(cols = final_fit) |>
    mutate(bal_acc_mw = map(.predictions,~bal_accuracy(.x,truth = class,estimate = .pred_class,estimator = "macro_weighted")),
           cm         = map(.predictions,~conf_mat(.x,class,.pred_class))) |>
    unnest(cols = bal_acc_mw) |>
    select(-c(id,.notes,.metric))

  #This is the raw performance
  # final_performance |>
  #   select(wflow_id,.metrics) |>
  #   unnest(cols = .metrics) |>
  #   filter(.metric == "bal_accuracy")

  boot_out$test_result[[j]] <- final_performance
  
  
  #VI ? Can do this for the elastic net models
  final_vi <- 
    final_performance |> 
    ungroup() |>
    filter(str_detect(wflow_id,"_en")) |>
    mutate(.penalty = map_dbl(best_params,~.x$penalty),
           mp       = map(.workflow,~extract_fit_parsnip(.x))) |>
    mutate(vi = map2(mp,.penalty,~vip::vi(.x,penalty =  .y))) |>
    select(wflow_id,.estimate,.penalty,vi)
  
  boot_out$test_vi[[j]] <- final_vi
  

  #Now make the bootstraps
  f_boot <-
    f_test |>
    select(id,class,sex:cluster_3_ev) |>

    # #Glue in the sensitivity analysis data for later
    # left_join(d_sens,by = "id") |>

    #Draw bootstraps
    bootstraps(times = n_boot,strata = class,apparent = FALSE) |>
    mutate(testing = purrr::map(splits,~analysis(.x)))|>
    dplyr::select(-splits)

  final_boot <-
    final_performance |>
    select(wflow_id,.workflow) |>
    mutate(boot = list(f_boot)) |>
    unnest(cols = boot) |>
    mutate(test_results = vector(mode = "list", length = 1))


  #Get the size of our new table (number of workflows x number of bootstraps)
  bl <- dim(final_boot)[1]



  #Run on all bootstraps
  tictoc::tic()


  
  for(i in 1:nrow(final_boot)){
    
    cat(magenta(
      'This is bootstrap ' %+%
        blue$underline$bold(i) %+%
        ' of ' %+%
        blue$underline$bold(nrow(final_boot)) %+%  '\n'
    ))
    
    
    final_boot$test_results[[i]] <- predict(final_boot$.workflow[[i]],final_boot$testing[[i]])
    
  }

  tictoc::toc()


  final_boot <-
    final_boot |>
    rowwise() |>
    mutate(
      .mcc = list(
        tibble(truth    = testing$class,
               estimate = test_results$.pred_class) |>
          yardstick::mcc(truth,estimate)),
      .ba = list(
        tibble(truth    = testing$class,
               estimate = test_results$.pred_class) |>
          yardstick::bal_accuracy(truth,estimate)      
        ))


  boot_metrics <-
    final_boot |>
    select(wflow_id,.ba,.mcc) |>
    pivot_longer(-wflow_id,names_to = '.metric',values_to = '.estimate') |>
    select(-.metric) |>
    unnest(cols = .estimate) |>
    select(-.estimator) |>
    drop_na(.estimate) |>
    reframe(conf_low  = coxed::bca(.estimate,conf.level = 0.95)[1],
            conf_high = coxed::bca(.estimate,conf.level = 0.95)[2],
            .estimate = mean(.estimate),
            .by = c(wflow_id,.metric)) |>
    relocate(.estimate,.before = conf_low)
  
  boot_metrics |>
    ggplot(aes(x = wflow_id,y = .estimate,ymin = conf_low,ymax = conf_high)) +
    geom_pointrange() +
    geom_hline(yintercept = 0,lty = 2) + 
    facet_wrap(~.metric,ncol = 2,scales = 'free_y')


  boot_out$boot_result[[j]] <- boot_metrics
 
}

#Save the bootstrapped data
write_rds(boot_out,paste("./Models/depression_subtype_test_performance.rds"))
#boot_out <- read_rds("./Models/depression_subtype_test_performance.rds")

## Plot results ------

wflow_name <- 
  c( "Full"                  = "full",
     "Sociodemographic"      = "soc",
     "Mental Health History" = "dep",
     "Immuno-metabolic"      = "blood")


metric_name <- 
  c( "Balanced Accuracy" = "bal_accuracy",
     "Brier Score"       = "brier_class",
     "MCC"               = "mcc",
     "AUROC"             = "roc_auc")

diagnosis_names <- 
  c( "ICD-10 Depression"   = "cisr_dep",
     "Depressive Symptoms" = "dep_plus",
     "Anxiety Symptoms"    = "anx_plus",
     "Somatic Symptoms"    = "som_plus")


#Inspect results - first raw performance
p_final_performance <- 
  boot_out |>  
  select(-boot_sens) |>
  unnest(cols = boot_result) |>
  select(symptom,wflow_id,.metric,.estimate,conf_low,conf_high) |> 
  filter(.metric == "bal_accuracy" ) |>
  # filter(.metric == "mcc" ) |>
  mutate(symptom_group = case_when(str_detect(symptom,"anx")   ~ "anxiety symptoms",
                                   str_detect(symptom,"dep")   & !str_detect(symptom,"cisr") ~ "depressive symptoms",
                                   str_detect(symptom,"som")   ~ "somatic symptoms",
                                   str_detect(symptom,"cisr")  ~ "depression diagnosis")) |>
  mutate(symptom = factor(symptom, levels = c("cisr_dep",
                                              "dep_plus",
                                              "anx_plus",
                                              "som_plus"))) |>
  filter(symptom %in% c("cisr_dep","dep_plus","anx_plus","som_plus")) |>
  mutate(symptom = fct_recode(symptom, !!!diagnosis_names)) |>
  separate_wider_delim(wflow_id, delim = "_",names = c("var_set","ml_method")) |>
  mutate(var_set = fct_recode(var_set, !!!wflow_name)) |>
  mutate(var_set = fct_relevel(var_set, c('Immuno-metabolic','Sociodemographic','Mental Health History','Full'))) |>
  
  ggplot(aes(y = var_set,x = .estimate,colour = var_set,xmin = conf_low,xmax = conf_high)) +
  geom_pointrange(position = position_dodge(width = 0.4)) +
  geom_vline(xintercept = c(0.5, 0.6,0.7,0.8),lty = 2) +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~symptom,ncol = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        #panel.grid = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "grey95"),
        legend.position = "none") +
  # coord_cartesian(xlim = c(0.47,0.75)) +
  labs(x = "Balanced Accuracy", y  = NULL)

#Performance is pretty much what you would expect
p_final_performance

ggsave('./Figures/depression_ml_test_performance.pdf',plot = p_final_performance,width = 3,height = 6)
  

#Tabulate outputs
boot_out |>  
  select(-boot_sens) |>
  unnest(cols = boot_result) |>
  select(symptom,wflow_id,.metric,.estimate,conf_low,conf_high) |> 
  mutate(symptom_group = case_when(str_detect(symptom,"anx")   ~ "anxiety symptoms",
                                   str_detect(symptom,"dep")   & !str_detect(symptom,"cisr") ~ "depressive symptoms",
                                   str_detect(symptom,"som")   ~ "somatic symptoms",
                                   str_detect(symptom,"cisr")  ~ "depression diagnosis")) |>
  mutate(symptom = factor(symptom, levels = c("cisr_dep",
                                              "dep_plus",
                                              "anx_plus",
                                              "som_plus"))) |>
  filter(symptom %in% c("cisr_dep","dep_plus","anx_plus","som_plus")) |>
  mutate(symptom = fct_recode(symptom, !!!diagnosis_names)) |>
  separate_wider_delim(wflow_id, delim = "_",names = c("var_set","ml_method")) |>
  mutate(var_set = fct_recode(var_set, !!!wflow_name)) |>
  mutate(var_set = fct_relevel(var_set, c('Immuno-metabolic','Sociodemographic','Mental Health History','Full'))) |>

  mutate(.metric = case_when(.metric == "bal_accuracy" ~ "Balanced Accuracy",
                             .metric == "mcc" ~ "MCC",
                             .metric == "roc_auc" ~ "AUROC"))  |>
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(estimate = paste0(.estimate," (",conf_low,", ",conf_high,")")) |>
  select(-c(conf_low,conf_high,.estimate,ml_method,symptom_group)) |>
  pivot_wider(names_from = .metric,values_from = estimate) |> 
  arrange(symptom,var_set) |>
  rename('Variable Set' = var_set,
         Symptom = symptom) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


## Plot final VI ######

#Let use the rank data as thats directly comparable and what Harrell advises in his textbook
p_final_vi <- 
  boot_out |>
  select(symptom,test_vi) |>
  unnest(test_vi) |>
  select(symptom,wflow_id,vi) |>
  unnest(vi) |>
  mutate(wflow_id = str_remove(wflow_id,'_en')) |>
  mutate(wflow_id = fct_recode(wflow_id,!!!wflow_name)) |>
  mutate(symptom = fct_recode(symptom, !!!diagnosis_names)) |>
  mutate(symptom = fct_relevel(symptom, c('ICD-10 Depression','Sociodemographic','Mental Health History','Full'))) |>
  ggplot(aes(x = Importance,
             colour = symptom,
             group = symptom,
             y = Variable)) +
  geom_col(position = position_dodge(width = 0.3)) +
  facet_grid(wflow_id~symptom) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        #panel.grid = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "grey80"),
        legend.position = "none") +
  labs(x = "Variable Importance", y = "Symptom") +
  scico::scale_color_scico_d(palette = "oleron")

#Performance is pretty much what you would expect
p_final_vi

#The problem is that the y axes are confusing and the plot is massive...
#gsave('./Figures/symptom_ml_test_vi.pdf',plot = p_final_vi,width = 18,height = 36)

#What are the top ten variables for ICD-10 depression diagnosis in the best performing model?
boot_out |>
  select(symptom,test_vi) |>
  unnest(test_vi) |>
  select(symptom,wflow_id,vi) |>
  unnest(vi) |>
  mutate(wflow_id = str_remove(wflow_id,'_en')) |>
  filter(symptom == "dep_plus" & wflow_id == "full") |>
  slice_max(Importance,n = 10)


#Just look at one variable set and just ICD-10 Depression e.g.
p_final_dep_vi <- 
  boot_out |>
  select(symptom,test_vi) |>
  unnest(test_vi) |>
  select(symptom,wflow_id,vi) |>
  unnest(vi) |>
  mutate(Variable = str_remove(Variable,'_f24')) |>
  mutate(Variable = str_remove(Variable,'_')) |>
  mutate(Variable = fct_reorder(Variable,Importance,max)) |>
  filter(wflow_id == "full_en") |>
  mutate(symptom = fct_recode(symptom, !!!diagnosis_names)) |>
  mutate(symptom = fct_relevel(symptom, c("ICD-10 Depression","Any Depressive Cognitions","Any Anxiety Symptoms","Any Somatic Symptoms"))) |>
  ggplot(aes(x    = Importance,
             y    = Variable)) +
  geom_col() +
  facet_wrap(~symptom,nrow = 1) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text.x.bottom = element_text(angle = 90,vjust=0.5),
        strip.text = element_text(size = 8),
        panel.background = element_rect(fill = "grey95")) +
  labs(x = "Variable Importance", y = "Variable") +
  scale_fill_brewer(palette = "Dark2")

p_final_dep_vi

ggsave('./Figures/diagnosis_ml_test_vi_combo.pdf',plot = p_final_dep_vi, width = 7,height = 7)


