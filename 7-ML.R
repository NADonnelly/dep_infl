
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

#Load up our pre-prepared data
d <- 
  read_rds(paste(data_dir,"alspac_data_final.rds",sep = '//'))



#Set our covariates
covars = c("sex","bmi24","smk24","audit_c","ph")
 

#Do some clearing up
rm(list = c("wdn","valid_vars_f24","valid_vars_f9","v_0"))



# Symptom scales ------

#We are interested in:

# - Are biomarkers any good for depression diagnosis?
# - Are biomarkers any good at diagnosing depression symptom profiles?

#We can do something similar to this covid paper - 
#https://www.nature.com/articles/s41590-024-01778-0#Abs1

#We define a  set of outcomes to try and predict. These are based on the 
#symptoms we find most associated with immuno-metabolic variables in the PLS, 
#which are also similar to the existing literature e.g. in the Milaneschi
#paper

#-Combination 1 (i.e. PLS component 1): FTG, SLP_DEC, MOT, and ACH
#-Combination 2 (i.e. PLS component 2): ANX and WOR

#Note that depressed mood (dep) is common to both and a requirement for a depression
#diagnosis itself so we do not include it in these subgroup profiles

d_symptom_score <- 
  d |>
  select(id,cisr_dep,ftg,slp_dec,mot,ach,anx,wor) |>
  mutate(across(.cols = ftg:wor,.fns = ~as.numeric(.x)-1)) |>
  transmute(id,cisr_dep,
            c1 = ftg+slp_dec+mot+ach,
            c2 = anx+wor)

#Look at how these synthetic scores are distributed
p_cisr_scores <- 
  d_symptom_score |>
  pivot_longer(-c(cisr_dep,id)) |>
  drop_na(cisr_dep) |>
  mutate(name = case_when(name == "c1" ~ "Somatic \nSymptoms",
                          name == "c2" ~ "Anxiety \nSymptoms")) |>
  mutate(cisr_dep = case_when(cisr_dep == 0 ~ "No Depression",
                              cisr_dep == 1 ~ "ICD10 Depression")) |>
  ggplot(aes(value,fill = cisr_dep)) +
  geom_histogram(binwidth = 1,colour = "black",linewidth = 0.25) +
  facet_wrap(~cisr_dep+name,ncol = 2, scales = "free") +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(0,16,2)) +
  scale_fill_manual(values = c("#32324B","#E1AF64"))+
  theme(panel.grid = element_blank(),
        axis.title  = element_text(size = 8),
        axis.text   = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(x = "Score", 
       y = "Count") +
  geom_vline(xintercept = c(4.5),lty = 2,linewidth = 0.25) +
  geom_vline(xintercept = c(8.5),lty = 2,linewidth = 0.25)

p_cisr_scores

#Look like a cut at 4 for anxiety and 8 for somatic would be reasonable for defining top/bottom halfs

#Save this for a supplementary figure
ggsave('./Figures/depression_subtype_scales.pdf',plot = p_cisr_scores,width = 6,height = 4)


#So we would be looking at the anxiety and somatic score domains within people with depression - 
#which is 337 individuals
d_symptom_score |>
  filter(cisr_dep == 1) |>
  drop_na()

#Counts by group
d_symptom_score|>
  filter(cisr_dep == 1) |>
  drop_na() |>
  mutate(anx_split      = if_else(c2 >= 5,1,0),
         som_split      = if_else(c1 >= 9,1,0))|>
  select(-c(id,cisr_dep,c1:c2)) |>
  pivot_longer(everything()) |>
  count(name,value) |>
  group_by(name) |>
  mutate(cc = cumsum(n))



#Do an upset plot to see overlap
p_overlap <- 
  d_symptom_score|>
  filter(cisr_dep == 1) |>
  drop_na() |>
  mutate(som_split  = if_else(c1 >= 9,1,0),
         anx_split  = if_else(c2 >= 5,1,0))|>
  select(-c(id,cisr_dep,c1:c2)) |>
  as.data.frame() |>
  upset(order.by = "freq")

p_overlap



#What is the relationship between our symptom scores and plain old depression sev
d_sev <- 
  d |>
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
  pivot_longer(c1:c2) |>
  mutate(name = case_when(name == "c1" ~ "Somatic \nSymptoms",
                          name == "c2" ~ "Anxiety \nSymptoms")) |>
  drop_na() |>
  ggplot(aes(x = dep_sev,y = value)) +
  geom_boxplot() +
  facet_wrap(~name,ncol = 3)

#So both subscales have higher scores in more severe depression, which you would expect.
#Its not simply a case that one is more related to severity than the others


# Prepare Data for ML ------

## Make the new outcomes -----

d_outcomes <- 
  d |>
  select(id,cisr_dep,cisr_dep_mod,cisr_dep_sev,ftg,slp_dec,mot,ach,anx,wor) |>
  mutate(across(.cols = c(cisr_dep,ftg:wor),.fns = ~as.numeric(.x)-1)) |>
  transmute(id,cisr_dep,
            c1 = ftg+slp_dec+mot+ach,
            c2 = anx+wor) |>
  drop_na() |>
  mutate(som_plus   = if_else((c1 >= 9) & (cisr_dep == 1),1,0),
         anx_plus   = if_else((c2 >= 5) & (cisr_dep == 1),1,0)) |>
  select(-c(c1:c2)) |>
  mutate(across(.cols = contains("_split")|contains("_any")|contains("_plus"),.fns = ~factor(.x)))


#Diagnosis - related

#CISR depression - 337 cases in the full data
d_outcomes |> count(cisr_dep)


# High Anxiety symptoms

# Among participants with ICD-10 depression: 140 no, 197 yes
d_outcomes  |> filter(cisr_dep == 1) |> count(anx_plus)

#High Somatic symptoms

# Among participants with ICD-10 depression: 179 no, 158 yes
d_outcomes  |> filter(cisr_dep == 1) |>  count(som_plus)



## Do the Data Prep ------

#Load the cluster data 
infl_clus <- 
  read_csv('./Models/f24_infl_variable_clusters.csv')

#Prepare the data for resampling
d_train <- 
  d |>
  mutate(bmi_strata = rsample::make_strata(bmi24,breaks = 4)) |>
  mutate(multi_demo = interaction(cisr_dep,sex,bmi_strata))  |>
  
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
  relocate(c(sex,bmi24,smk24,audit_c,ph,sc,me,md),.after = id) |>

  
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


    #Keep this for stratification later
    multi_demo) |>
  
  mutate(dep_17 = if_else((cisr_mild_17 + cisr_mod_17 + cisr_sev_17) > 0,1,0 )) |>
  mutate(dep_17 = as.integer(dep_17))|>
  select(-c(cisr_mild_17,cisr_mod_17,cisr_sev_17)) |>
  
  #get the ace varaibles
  left_join(d_0 |> 
              select(alnqlet,clon122,clon170) |>
              rename(id = alnqlet,
                     ace_total_classic  = clon122,
                     trauma_total       = clon170),
            by = "id") |>
  
  #Now you need to graft on your outcomes
  left_join(d_outcomes,by = "id") |>
  relocate(cisr_dep:anx_plus,.after = id)



#Now, the size of our final training dataset is
dim(d_train)



## Variables and recipes -----

#Get variable types
d_var_type = 
  d_train |>
  
  #Remove our categorical variables
  select(-c(id:anx_plus,sc, me,sex,smk24,ph,multi_demo)) |>
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



### Set 1: Bloods + Covariates -----

#Take the cluster 2 and cluster 3 variables for consistency with other 
#aspects of the paper

f24_blood_names <- 
  infl_clus |> 
  filter(clus %in% c(2,3)) |> 
  mutate(infl = paste0(infl,"_f24")) |>
  pull(infl)

iv_f24_blood <- 
  c("sex","audit_c","smk24","bmi24","ph",
    f24_blood_names,
    "cluster_3_ev")

#Note we include our covariates in the model as we can see that IM-variables are
#pretty good at classifying on our covariates e.g. you can get sex from the Hb
#so we want the model the get information about depression from the IM variables
#and it can get sex data from the sex variable

#Update the continuous variables to include only those variables we are 
#including
cont_vars <- intersect(iv_f24_blood,cont_vars)



### Set 2: Social Variables -----

#Make a set of social variables we know are likely to be related to
#depression but are not blood-related, to act as a positive control
iv_soc <- 
  c("sex","audit_c","smk24","bmi24", # Demographic variables
    "sc","ph",                       # Maternal Social class- correlates 
                                     # substantially with maternal education so 
                                     # we will only pick one
    "ace_total_classic",             # ACES
    "trauma_total"                   # Trauma
  )


### Set 3: Past MH variables -----

#Make a set of variables related to personal history of depression as another 
#positive control

iv_dep <- 
  c("md",   # Maternal EPDS
    d_train |> select(starts_with("smfq")) |> colnames(),  #Earlier smfq data
    d_train |> select(starts_with("sdq")) |> colnames(),   #Earlier sdq data - use age 7 and 9
    "dep_17")


### Set 4: Full dataset -----

#Make a full set of variables combining the blood and clinical variables
iv_full <- 
  c(iv_soc,iv_dep,iv_f24_blood[str_detect(iv_f24_blood,"_f24")],"cluster_3_ev")

### Formulas for the recipes --------

#Make formulas 
f_blood_24  <- as.formula(paste0('class ~ ', paste0(iv_f24_blood , collapse = ' + '), collapse = ""))
f_soc       <- as.formula(paste0('class ~ ', paste0(iv_soc       , collapse = ' + '), collapse = ""))
f_dep       <- as.formula(paste0('class ~ ', paste0(iv_dep       , collapse = ' + '), collapse = ""))
f_full      <- as.formula(paste0('class ~ ', paste0(iv_full      , collapse = ' + '), collapse = ""))


#Also make a list of variables excluding the evs, which we make inside loops
iv_pre_ev <- 
  c(iv_soc,iv_dep,iv_f24_blood[str_detect(iv_f24_blood,"_f24")])

## Models -----

#Model specs

### Elastic net -----
en_mod <- 
  logistic_reg(penalty = tune(), 
               mixture = tune()) %>% 
  set_engine("glmnet") %>% 
  set_mode("classification")


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


# Fit model -------

#Fit our models using nested cross validation, nesting over each symptom
dv_list <- c("cisr_dep",
             "som_plus", 
             "anx_plus")

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
    
    
    #First of all lets make that cluster 3 ev metric in the training and then then
    #test data
    d_outer_1 <- 
      d_outer |>
      analysis() 
    
    d_outer_2 <- 
      d_outer |>
      assessment()
    
    
    #EV model
    
    #First apply the model to the training data for this outer loop
    d_res_model <- 
      left_join(
        d_train |>
          select(id,cisr_dep,contains('_f24'),all_of(covars)) |>
          drop_na(cisr_dep) |>
          filter(cisr_dep == 0) |>
          select(-cisr_dep) |>
          pivot_longer(-c(id,all_of(covars)),names_to = 'infl',values_to = 'infl_value') |>
          drop_na() |>
          nest_by(infl,.key = "data_nd"),
        
        d_train |>
          select(id,cisr_dep,contains('_f24'),all_of(covars)) |>
          drop_na(cisr_dep) |>
          filter(cisr_dep == 1) |>
          select(-cisr_dep) |>
          pivot_longer(-c(id,all_of(covars)),names_to = 'infl',values_to = 'infl_value') |>
          drop_na() |>
          nest_by(infl,.key = "data_d"),
        by = "infl") |>
      mutate(infl = str_remove(infl,"_f24")) |>
      
      #Only for the cluster 3 variables
      filter(infl %in% c(infl_clus |> filter(clus == 3) |> pull(infl))) |>
      
      mutate(m1 = 
               list(
                 lm(infl_value ~ sex + bmi24 + smk24 + audit_c,data = data_nd)
               ),
             r1 = 
               list(rstudent(m1))  ) |>
      mutate(
        pred_nd = 
          list(
            predict.lm(m1,newdata = data_nd) |>
              as_tibble() |>
              rename(pred = value) |>
              bind_cols(data_nd) |>
              mutate(resid = infl_value-pred)
          ),
        pred_d = 
          list(
            predict.lm(m1,newdata = data_d) |>
              as_tibble() |>
              rename(pred = value)|>
              bind_cols(data_d) |>
              mutate(resid = infl_value-pred)
          ),
        mu_nd = list(mean(pred_nd$resid,na.rm = T)),
        sd_nd = list(sd(pred_nd$resid,na.rm = T))
      ) |>
      ungroup() |>
      unnest(cols = c(mu_nd,sd_nd))
    
    #Now we have the models and the values to use for standardisation
    
    #Next step is to get the extreme value counts in the training data
    d_outer_1_ev <- 
      d_res_model |>
      select(infl,pred_nd,pred_d,mu_nd,sd_nd) |>
      pivot_longer(starts_with("pred_")) |>
      unnest(value) |>
      
      #Standardise the residuals using the mean and sd of residuals from the non-depressed group
      mutate(resid_z = (resid-mu_nd)/sd_nd) |>
      
      mutate(ex_dev = if_else(abs(resid_z) > 2.6,1,0) |> as.integer())  |>
      select(infl,id,all_of(covars),resid_z,ex_dev) |>
      group_by(id) |>
      reframe(cluster_3_ev = sum(ex_dev))

    #Glue this onto the outer loop training data
    d_outer_1 <- 
      d_outer_1 |>
      left_join(d_outer_1_ev,by = "id")
    
    #Next we apply the same model to the outer loop testing data
    d_outer_2_ev <- 
      d_outer_2 |>
      select(id,cisr_dep,contains('_f24'),all_of(covars)) |>
      pivot_longer(-c(id,all_of(covars)),names_to = 'infl',values_to = 'infl_value') |>
      mutate(infl = str_remove(infl,"_f24")) |>
      filter(infl %in% c(infl_clus |> filter(clus == 3) |> pull(infl))) |>
      nest_by(infl,.key = "data_t") |>
      
      left_join(d_res_model |>
                  select(infl,m1,mu_nd,sd_nd),
                         by = "infl")|>
      mutate(
        pred_t = 
          list(
            predict.lm(m1,newdata = data_t) |>
              as_tibble() |>
              rename(pred = value) |>
              bind_cols(data_t) |>
              mutate(resid = infl_value-pred)
          )
      )|>
      select(infl,pred_t,mu_nd,sd_nd) |>
      pivot_longer(starts_with("pred_")) |>
      unnest(value) |>
      
      #Standardise the residuals using the mean and sd of residuals from the non-depressed group
      mutate(resid_z = (resid-mu_nd)/sd_nd) |>
      
      mutate(ex_dev = if_else(abs(resid_z) > 2.6,1,0) |> as.integer())  |>
      select(infl,id,all_of(covars),resid_z,ex_dev) |>
      group_by(id) |>
      reframe(cluster_3_ev = sum(ex_dev))
    
    #Glue this onto the outer loop training data
    d_outer_2 <- 
      d_outer_2 |>
      left_join(d_outer_2_ev,by = "id")
    
    
    #Dataset prep
    
    
    #Now we do some prep to rename the dependent variable and make sure its a factor
    #This is somewhat different depending on the data.
    
    if(str_detect(dv,'_plus')){
      
      d_outer_train <- 
        d_outer_1 |>
        filter(cisr_dep == 1) |>
        dplyr::select(-c(id,multi_demo,cisr_dep)) |>
        rename(class = !!dv) |>
        mutate(class = factor(class)) |>
        dplyr::select(class, all_of(iv_full)) |> 
        drop_na(class)
      
      d_outer_test <- 
        d_outer_2 |>
        filter(cisr_dep == 1) |>
        dplyr::select(-c(id,multi_demo,cisr_dep)) |>
        rename(class = !!dv) |>
        mutate(class = factor(class)) |>
        dplyr::select(class, all_of(iv_full))|> 
        drop_na(class)
      
    }else{
      
      d_outer_train <- 
        d_outer_1 |>
        dplyr::select(-c(id,multi_demo)) |>
        rename(class = !!dv) |>
        mutate(class = factor(class)) |>
        dplyr::select(class, all_of(iv_full)) |> 
        drop_na(class)
      
      d_outer_test <- 
        d_outer_2 |>
        dplyr::select(-c(id,multi_demo)) |>
        rename(class = !!dv) |>
        mutate(class = factor(class)) |>
        dplyr::select(class, all_of(iv_full))|> 
        drop_na(class)
      
    }
    

    #Finally make the inner bootstraps from the outer loop training data
    d_inner <- 
      d_outer_train |>
      bootstraps(times = n_inner)
    


    #We now rebuild a split object with the prepared data
    d_outer_split =
      make_splits(d_outer_train,d_outer_test)
    
    
    #Now do the recipe prep
    
    
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
                "d_outer_1","d_outer_1_ev","d_outer_2","d_outer_2_ev",
                "d_res","d_res_model",
                "rec_impute","rec_full","rec_blood_24","rec_soc","rec_dep",
                "wf_par","wf_results",
                "best_results_wf",
                "roc_mod","roc_mod_post"))  
    
  }
  
  
  ml_result$ml[[j]] = d_fold_results
  

  rm(list = c("d_fold_results"))
  
}

write_rds(ml_result,paste("./Models/depression_ml_results.rds"))
#ml_result <- read_rds("./Models/depression_ml_results.rds")

# Inspect output -----

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
            'cor_vars','d_test_id','cor_train',
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
  # filter(.metric == "bal_accuracy") |>
  separate_wider_delim(wflow_id,delim = '_',names = c("vars","mod_type")) |>
  mutate(symptom = factor(symptom,levels = c("cisr_dep","anx_plus","som_plus"))) |>
  ggplot(aes(x = symptom,y = .estimate,shape = vars,colour = vars,
             group = vars)) +
  geom_point(position = position_dodge(width = 0.3)) +
  facet_grid(mod_type~.metric) +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_text(angle = -90)) +
  geom_hline(yintercept = c(0,.5),lty = 2) +
  scale_y_continuous(breaks = seq(-0.2,1,0.2)) +
  labs(x = "Model",y = "Metric")

## Performance model ------

#Prepare results for modelling
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

symptom_name <- 
  c("ICD-10 Depression" = "cisr_dep",
    "Anxiety Symptoms"  = "anx_plus",
    "Somatic Symptoms"  = "som_plus")


#Sort out prep for plotting
cv_mod <- 
  cv_mod |>
  mutate(.metric = fct_recode(.metric, !!!metric_name)) |>
  mutate(.metric = fct_relevel(.metric, c("Brier Score","Balanced Accuracy","AUROC","MCC"))) |> 
  
  mutate(symptom = fct_recode(symptom, !!!symptom_name)) |>
  mutate(symptom = fct_relevel(symptom, c("ICD-10 Depression","Anxiety Symptoms","Somatic Symptoms"))) |>  
  
  mutate(disp_line = case_when(.metric == "AUROC" ~ 0.5,
                               .metric == "Balanced Accuracy" ~ 0.5,
                               .metric == "Brier Score" ~ 0.25,
                               .metric == "MCC" ~ 0),
         disp_lims = case_when(.metric == "AUROC" ~ c(0.45,0.85) |> list(),
                               .metric == "Balanced Accuracy" ~ c(0.45,0.75) |> list(),
                               .metric == "Brier Score" ~ c(0,0.3) |> list(),
                               .metric == "MCC" ~ c(-0.05,0.4) |> list()))  |>

  arrange(.metric) |>
  rowwise(.metric) 


d_cv_mod_plot <- 
  cv_mod |>
  # filter(.metric %in% c("MCC")) |>
  # filter(.metric %in% c("Balanced Accuracy")) |>
  arrange(symptom) |>
  mutate(draws = 
           list(model |> 
                  tidybayes::tidy_draws()  |> 
                  select(.chain:.draw,starts_with('wflow_id')) |>
                  pivot_longer(starts_with("wflow_id"))  |> 
                  mutate(name = str_remove(name,"wflow_id")))) |>
  ungroup() |>
  select(symptom,.metric,draws,disp_line,disp_lims) |>
  unnest(cols = draws) |>
  separate_wider_delim(name, delim = "_",names = c("var_set","ml_method")) |>
  mutate(var_set = fct_recode(var_set, !!!wflow_name)) |>
  mutate(var_set = fct_relevel(var_set, c('Immuno-metabolic','Sociodemographic','Mental Health History','Full')))

#Make a plot just of the balanced accuracy
p_all_wflow <- 
  d_cv_mod_plot |>
  filter(.metric == "Balanced Accuracy" & ml_method == "en") |>
  ggplot(aes(y = var_set, x = value,fill = var_set,colour = var_set)) +
  #ggdist::stat_halfeye()+
  ggdist::stat_pointinterval(point_interval = "mean_hdci",
                             .width = c(0.95))+
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        #panel.grid = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "grey95"),
        legend.position = "none") +
  labs(x = "Balanced Accuracy", 
       y = "Variable Set") +
  coord_cartesian(xlim = c(0.4,0.8)) +
  geom_vline(xintercept = c(0, 0.5),lty = 2) +
  #scico::scale_fill_scico_d(palette = "oleron") + 
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~symptom,ncol = 1)
  #facet_grid(symptom~ml_method)


p_all_wflow

#Save for plotting
ggsave('./Figures/depression_ml_performance.pdf',plot = p_all_wflow,width = 3,height = 4)


#Within the people with depression the IM model does better for discriminating high somatic symptoms 
#than high depression and anxiety symptoms, MH history is much better for depressive symptoms


## Tabulate output ======

d_cv_mod_plot |>  
  select(-c(disp_line,disp_lims)) |>
  group_by(.metric,symptom,var_set) |>
  ggdist::mean_hdci(value,.width = 0.95) |>
  select(-c(.width,.point,.interval)) |>

  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(estimate = paste0(value," (",.lower,", ",.upper,")")) |>
  select(-c(.lower,.upper,value)) |>
  arrange(symptom,var_set) |>
  rename('Variable Set' = var_set,
         Symptom = symptom) |>
  pivot_wider(names_from = .metric,values_from = estimate) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


## Directly compare variable sets ------

#We can use the posterior samples to make a direct comparison between the immuno-metabolic variable's predictive
#accuracy for the psychic vs somatic symptom profiles:
cv_mod |>
  filter(.metric %in% c("Balanced Accuracy") & symptom %in% c("Anxiety Symptoms","Somatic Symptoms")) |>
  mutate(draws = 
           list(model |> 
                  tidybayes::tidy_draws()  |> 
                  select(.chain:.draw,starts_with('wflow_id')) |>
                  pivot_longer(starts_with("wflow_id"))  |> 
                  mutate(name = str_remove(name,"wflow_id")))) |>
  ungroup() |>
  select(symptom,.metric,draws) |>
  unnest(cols = draws) |>
  separate_wider_delim(name, delim = "_",names = c("var_set","ml_method")) |>
  mutate(var_set = fct_recode(var_set, !!!wflow_name)) |>
  mutate(var_set = fct_relevel(var_set, c('Immuno-metabolic','Sociodemographic','Mental Health History','Full'))) |>
  filter(var_set == "Immuno-metabolic" & ml_method == "en" & .metric %in% c("Balanced Accuracy")) |>
  select(-c(.metric,var_set,ml_method)) |>
  pivot_wider(names_from = symptom,values_from = value) |>
  mutate(value_diff = `Somatic Symptoms` - `Anxiety Symptoms`) |>
  ggdist::mean_hdci(value_diff,width = 0.95) |>
  as_tibble() |>
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3)))
