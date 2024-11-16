# Introduction --------

#Sometimes its just not very nice having to make markdown documents. Lets 
#reproduce the alspac-olink-basics.qmd document here but we will save the figures
#in nice ways, and we might work out how to export the tables nicely too

pacman::p_load(broom,
               tidyverse,
               tidymodels,
               rsample,
               viridis,
               janitor,
               marginaleffects,
               mixOmics,
               tidybayes,
               patchwork,
               gtsummary)



tidymodels_prefer()

conflicted::conflict_prefer('select','dplyr')
conflicted::conflict_prefer('filter','dplyr')
conflicted::conflict_prefer('map'   ,'purrr')
conflicted::conflicts_prefer(crayon::`%+%`)

# Load our datasets - invoke a script that loads from the secure RDSF folder
source('./Final/0-Common_data.R')

#Load up our pre-split data so we can focus on the training data only
d_split <- 
  read_rds(paste(data_dir,"alspac_data_split.rds",sep = '//'))

d_train <- 
  d_split |>
  training()

#Set our covariates
covars = c("sex","bmi24","smk24","audit_c")

## Data Dimensions ---------

# The total size of the alspac dataset is:
dim(d_0)[1]

#Of this total cohort, this many have Olink F24 panel bloods:
d_24 |> select(il6) |> drop_na() %>% dim(.) %>% .[1]

#Of the total ALSPAC cohort, this many have have complete CISR data
d_24 |> select(cisr_dep) |> drop_na() %>% dim(.) %>% .[1]

#Of whom, this many have an ICD-10 diagnosis of depression
d_24 |> select(cisr_dep) |> drop_na() |> count(cisr_dep) |> filter(cisr_dep ==1) |> pull(n)

#A total of this many participants have data for both age 24 Olink panel bloods 
#and CISR data
d_24 |> select(cisr_dep,il6) |> drop_na()%>% dim(.) %>% .[1]

#And this many in that sample have depression
d_24 |> select(cisr_dep,il6)  |> drop_na() |> count(cisr_dep) |> filter(cisr_dep ==1) |> pull(n)

#Including Olink panel bloods and the existing F24 blood test data, we have a 
#total of this many blood test variables.
d_24 |> select(-c(id,cisr_dep)) |> colnames() |> length()



#Now, we have our training data - what about its dimensions?
dim(d_train)[1] #3331

#Of the training set, this many have Olink F24 panel bloods:
d_train |> select(il6_f24) |> drop_na() %>% dim(.) %>% .[1]

#And this many had a CISR dataset
d_train |> select(cisr_dep) |> drop_na() %>% dim(.) %>% .[1]

#Depression count - 
d_train |> select(cisr_dep) |> drop_na() |> count(cisr_dep) |> filter(cisr_dep ==1) |> pull(n)

d_train |> select(cisr_dep) |> drop_na() |> count(cisr_dep) |> mutate(c = cumsum(n),p = n/c)

#Make a pseudo "complete case" n for tabulating based on having both a CISR value and at 
#least one immunometabolic variable
d_cc <- 
  d_train |> 
  select(id,c(contains("_f24"),cisr_dep)) |> 
  mutate(across(.cols = -c(id,cisr_dep),.fns = ~if_else(is.na(.x),1,0))) |>
  mutate(sum = rowSums(across(contains("_f24")))) |>
  select(id,cisr_dep,sum) |>
  filter(sum < 93) |>
  drop_na(cisr_dep) 



# F24 Covariates ------

d_24_cv = 
  d_cv |>
  left_join(d_train |> select(id,cisr_dep),
            by = "id") |>
  select(id,cisr_dep,all_of(covars)) |>
  filter(id %in% d_cc$id)

dim(d_24_cv)[1]

## Covariate Data ------ 

#From the set of participants with CISR data, this many were male:
d_24_cv |> count(sex) |> filter(sex == "Male") |> pull(n)

#And this many were female:
d_24_cv |> count(sex) |> filter(sex == "Female") |> pull(n)

#So rather female

#Median BMI at age 24 by sex was this for males:
d_24_cv |> drop_na(bmi24) |> group_by(sex) |>  summarise(bmi_m = median(bmi24) |> janitor::round_half_up(digits = 2)) |>  filter(sex == "Male") |>  pull(bmi_m)

#And this for females:
d_24_cv |> drop_na(bmi24) |>  group_by(sex) |>  summarise(bmi_m = median(bmi24) |> janitor::round_half_up(digits = 2)) |>  filter(sex == "Female") |>  pull(bmi_m)

#A total of this many males smoked daily
d_24_cv |> drop_na(smk24) |> count(sex,smk24) |>  filter(sex == "Male") |> arrange(-smk24) |> mutate(n = cumsum(n))|> pull(n) |> paste0(collapse = '/')

#And this many females:
d_24_cv |> drop_na(smk24) |> count(sex,smk24) |>  filter(sex == "Female") |> arrange(-smk24) |> mutate(n = cumsum(n))|> pull(n) |> paste0(collapse = '/')

#Median AUDIT-C by sex was this for males:
d_24_cv |> drop_na(audit_c) |> group_by(sex) |>  summarise(au_m = median(audit_c)) |>  filter(sex == "Male") |>  pull(au_m)

#And this for females: 
d_24_cv |>  drop_na(audit_c) |>  group_by(sex) |>  summarise(au_m = median(audit_c)) |>  filter(sex == "Female") |>  pull(au_m)


## Covariate Table -----

#We can make a nice looking table 1
d_24_cv %>% 
  select(-id) %>%
  tbl_summary(by = c(sex)) %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**sex**") %>%
  bold_labels()

#There you go

#And we have our model DAG too, which is found in 1-Model_DAG.R


## Covariate Plot ------

#This is going to go into the supplement with the DAG

#Make a nice plot
p_sex = 
  d_24_cv |>
  count(sex) |>
  ggplot(aes(x = sex, y = n,fill = sex)) +
  geom_col(colour = "black") +
  theme(panel.grid = element_blank(),
        title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_fill_brewer(palette = 2) +
  labs(x = "Sex", y = "n @ Age 24",title = "Sex of Participants")

p_bmi24 = 
  d_24_cv |>
  drop_na(bmi24) |>
  ggplot(aes(x = bmi24, fill = sex)) +
  geom_histogram(binwidth = 2,colour = "black") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_fill_brewer(palette = 2) +
  facet_wrap(~sex,ncol = 1,scales = "free_y") +
  labs(x = "BMI @ Age 24", y = "n @ Age 24", title = "BMI of Participants")

p_smoking = 
  d_24_cv |>
  drop_na(smk24) |>
  mutate(smk24 = case_when(smk24 == 0 ~ "No",
                           smk24 == 1 ~ "Yes")) |>
  count(sex,smk24)  |>
  ggplot(aes(x = smk24, y = n,fill = sex)) +
  geom_col(colour = "black",position = position_dodge()) +
  theme(panel.grid = element_blank(),
        title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_fill_brewer(palette = 2) +
  labs(x = "Smoking Status", y = "n @ Age 24",title = "Daily Smoking in Participants")

p_drinking = 
  d_24_cv |>
  drop_na(audit_c) |>
  ggplot(aes(x = audit_c, fill = sex)) +
  geom_histogram(binwidth = 1,colour = "black") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_fill_brewer(palette = 2) +
  facet_wrap(~sex,ncol = 1) +
  labs(x = "AUDIT-C Score", y = "n @ Age 24",title = "Drinking in Participants")

#Lets also plot our depression variable
p_cisr = 
  d_24_cv |>
  drop_na(cisr_dep) |>
  mutate(cisr_dep = case_when(cisr_dep == 0 ~ "No",
                              cisr_dep == 1 ~ "Yes")) |>
  count(sex,cisr_dep)  |>
  ggplot(aes(x = cisr_dep, y = n,fill = sex)) +
  geom_col(colour = "black",position = position_dodge()) +
  theme(panel.grid = element_blank(),
        title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_fill_brewer(palette = 2) +
  labs(x = "ICD-10 Depression", y = "n @ Age 24",title = "ICD-10 Depression in Participants")


p_cisr_prop <- 
  d_24_cv |>
  drop_na(cisr_dep) |>
  mutate(cisr_dep = case_when(cisr_dep == 0 ~ "No",
                              cisr_dep == 1 ~ "Yes")) |>
  count(sex,cisr_dep) |>
  group_by(sex) |>
  mutate(cs = sum(n),
         prop = n/cs)|>
  ggplot(aes(x = cisr_dep, y = prop,fill = sex)) +
  geom_col(colour = "black",position = position_dodge()) +
  theme(panel.grid = element_blank(),
        title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_fill_brewer(palette = 2) +
  labs(x = "ICD-10 Depression", y = "Proportion of Participants",title = "ICD-10 Depression as Proportion")


#Assemble 
p_f24_covars = 
  ((p_sex/p_smoking) | p_bmi24 |p_drinking|(p_cisr/p_cisr_prop)) + 
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

#And save
ggsave(filename = "./Figures/f24_covars.pdf",
       height = 3,
       width  = 8,
       plot = p_f24_covars)

#This will need a little wrangling in inkscape to get everything nicely aligned


# F24 Correlation Plot -------

#Make a correlation plot of the F24 bloods to replicate the Wellcome data note
cor_f24 = 
  d_train |>
  select(contains("_f24")) |>
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  mutate(across(.cols = everything(),.fns = zscore)) |>
  cor(use = "pairwise.complete.obs",
      method = "pearson") 

#Do our own clustering so we can save the variable names for later - I have 
#taken this code directly from the ggcorrplot github so it should produce 
#identical results to using the hc.ord method
clust_cor_f24 <- 
  hclust(as.dist(1 - cor_f24),method = "complete")

#Extract the ordering based on the clustering
ord_f24 <- clust_cor_f24$order

#Save variable names for later
vars_ordered <- clust_cor_f24$labels[clust_cor_f24$order]

#Sort the variables
cor_f24 <- cor_f24[ord_f24, ord_f24]

#Now make the plot
p_f24_cor <- 
  cor_f24 %>% 
  ggcorrplot::ggcorrplot(method = "square",type = "full",
                         show.diag = TRUE,
                         #colors = scico::scico(n = 3,palette = 'roma',direction = -1),
                         hc.order = F,
                         tl.cex = 6) +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_text(angle = 90),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  scico::scale_fill_scico(palette = 'vik',name = "Correlation",direction = 1,
                          limits = c(-1,1))

p_f24_cor


#F24 Models --------

## Depression Group Difference -----

#Now we fit our models, first as a complete case analysis with the inflammatory
#variable as the DV, so we are calculating group differences between depressed
#and non-depressed

d_reg <- 
  d_train |>
  select(id,cisr_dep,all_of(covars),contains("_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  
  #Sort out variable names
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  
  #Put each blood variable into its own dataset
  pivot_longer(-c(id,cisr_dep,all_of(covars)), names_to = "infl", values_to = "infl_value") |> 
  drop_na() |>
  nest_by(infl)

#So the n for these models varies between
d_reg |> mutate(nd = nrow(data) |> list()) |> ungroup() |> reframe(r = range(nd))



## Frequentist Models ------

#Fit our models
d_reg = 
  d_reg |>
  mutate(m1 = lm(formula = infl_value ~ cisr_dep, 
                 data = data) |>
           list(),
         m2 = lm(formula = infl_value ~ cisr_dep * (sex + audit_c + smk24), 
                 data = data) |>
           list(),
         m3 = lm(formula = infl_value ~ cisr_dep * (sex + audit_c + smk24 + bmi24), 
                 data = data) |>
           list()) 

d_reg = 
  d_reg |>
  select(-c(data)) |>
  pivot_longer(-infl,names_to = "model_name",values_to = "model") |>
  mutate(me = map(model,~avg_comparisons(model = .x,
                                         variables = list(cisr_dep = 0:1),conf_level = 0.95) |>
                    
                    #This stuff is to get marginaleffects to behave itself - at one point 
                    #its raw output used huge amountso f memeory for some reason
                    as_tibble() |>
                    as.matrix() |>
                    as_tibble() |>
                    mutate(across(.cols = -c(term:contrast), .fns = ~as.double(.)))
                  )
         )

#Quickly tabulate our model results
d_reg_result <- 
  d_reg |>
  mutate(me_tab = map(me, ~as_tibble(.x) |> 
                        select(estimate:conf.high) |>
                        mutate(across(.cols = everything(),.fns = ~as.double(.x))))) |>
  select(-c(model:me)) |>
  unnest(me_tab) |>
  arrange(p.value) |>
  group_by(model_name)|>
  mutate(p.adj = p.adjust(p.value,method = 'fdr')) |>
  ungroup()



### Table of results -----

#Get the top markers for the base model
sig_markers = 
  d_reg_result |>
  filter(model_name == "m1") |>
  filter(p.adj < 0.05) |>
  pull(infl)

#Get the variable descriptions
sig_desc <- 
  v_0 |> 
  as_tibble() |> 
  select(name,lab) |> 
  filter(str_detect(name,"_F24")) |>
  mutate(name = tolower(str_remove(name,"_F24"))) |> 
  filter(name %in% sig_markers) |>
  rename(infl = name)


#Tabulate  
d_reg_result |>
  ungroup() |>
  filter(infl %in% sig_markers) |>
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(estimate = paste0(estimate," [",conf.low,", ",conf.high,"], P FDR = ",p.adj)) |>
  arrange(infl) |>
  select(-c(std.error:p.adj)) |>
  mutate(model_name = case_when(model_name == "m1" ~ "Model 1: Unadjusted",
                                model_name == "m2" ~ "Model 2: Adjusted (no BMI)",
                                model_name == "m3" ~ "Model 3: Adjusted (including BMI)")) |>
  pivot_wider(names_from = model_name,values_from = estimate) |>
  relocate(`Model 3: Adjusted (including BMI)`,.after = `Model 2: Adjusted (no BMI)`) |>
  #left_join(sig_desc,by = "infl") |>
  rename(`Immunometabolic Variable` = infl) |> 
 # relocate(lab,.after = `Immunometabolic Variable`) |>
  #rename(`Description` = lab) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)

#Note that the model coefficients represent the marginal change in the 
#inflammatory variable of interest between non-depressed and depressed 
#participants (in SD units)

#Tabulate all associations
d_reg_result |>
  ungroup() |>
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(estimate = paste0(estimate," [",conf.low,", ",conf.high,"], PFDR = ",p.adj)) |>
  arrange(infl) |>
  select(-c(std.error:p.adj)) |>
  mutate(model_name = case_when(model_name == "m1" ~ "Model 1: Unadjusted",
                                model_name == "m2" ~ "Model 2: Adjusted (no BMI)",
                                model_name == "m3" ~ "Model 3: Adjusted (including BMI)")) |>
  pivot_wider(names_from = model_name,values_from = estimate) |>
  relocate(`Model 3: Adjusted (including BMI)`,.after = `Model 2: Adjusted (no BMI)`) |>
  #left_join(sig_desc,by = "infl") |>
  rename(`Immunometabolic Variable` = infl) |> 
  # relocate(lab,.after = `Immunometabolic Variable`) |>
  #rename(`Description` = lab) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


### Plot ----

p_reg <- 
  d_reg_result |>
  filter(infl %in% sig_markers) |>
  mutate(model_name = case_when(model_name == "m1" ~ "Unadjusted",
                                model_name == "m2" ~ "Adjusted - BMI",
                                model_name == "m3" ~ "Adjusted + BMI")) |>
  mutate(model_name = factor(model_name,levels = c("Unadjusted","Adjusted - BMI","Adjusted + BMI"))) |>

  ggplot(aes(y = infl, x = estimate, xmin = conf.low,xmax = conf.high,
             colour = model_name))+
  geom_pointrange(position = position_dodge(width = 0.5))+
  scale_color_brewer(palette = "Set1")+
  geom_vline(xintercept = 0,lty = 2,linewidth = 0.25) +
  theme(panel.grid = element_blank(),
        axis.title  = element_text(size = 8),
        axis.text   = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_blank()) +
  labs(y = "Blood Variable", 
       x = "Change in standardised measure")


ggsave(filename = "./Figures/f24_lineplots.pdf",
       height = 5,
       width  = 5,
       plot = p_reg)

### Hb -----

# A quick note on haemoglobin; this comes out as significant in unadjusted 
#models, however, we can see that Hb is strongly related to sex at birth:
d_hb <- 
  d_train |>
  select(id,all_of(covars),cisr_dep,hb_f24) |>
  drop_na()

#Plot Hb Distribution by sex
d_hb |>
  ggplot(aes(hb_f24,fill = sex)) +
  geom_density(alpha = 0.4,colour = "black")  +
  theme(panel.grid = element_blank(),
        title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_fill_brewer(palette = 2) +
  labs(x = "Hb at age 24", y = "Density",title = "Hb by Sex")

#Fit a linear model on Hb
m_hb <- 
  lm(hb_f24 ~ sex, data = d_hb)

#So what are the average Hb values for men and women?
avg_predictions(m_hb, by = "sex") |> as_tibble()
avg_comparisons(m_hb) |> as_tibble()

#On average women have Hb 17.6 mg/dL lower than men

#What about depression?
m_hb2 <- 
  glm(cisr_dep ~ sex,data = d_hb,
      family = binomial(link = "logit"))
  
avg_predictions(m_hb2, by = "sex") |> as_tibble()
avg_comparisons(m_hb2,by = TRUE, type = "response") |> as_tibble()

#On average women have a 6.2% higher rate of depression than men

#Finally, lets stratify our full model by sex
m_hb_strat <- 
  d_hb |>
  nest_by(sex) |>
  mutate(
    m3 = list(
      lm(hb_f24 ~ cisr_dep * (audit_c + smk24 + bmi24),
         data = data)
      ),
    me = list(
      avg_comparisons(model = m3,
                      variables = list(cisr_dep = 0:1),conf_level = 0.95) |>
        as_tibble()
    ))

m_hb_strat |>
  select(sex,me) |>
  unnest(cols = me) |>
  ungroup()

## Multiple Imputation Results -------


#See 3b and 3c-MI_Modelling.R for the code that produces these results; it takes
#about 2 days of computation on an OK laptop to run so we don't run it here.

d_dep_tab <- 
  read_csv(paste(data_dir,"alspac_split_imputed_results.csv",sep = '/'))

d_mi <- 
  d_dep_tab |>
  rename(infl = vars) |>
  mutate(infl = map_chr(infl,~str_remove(.x,"_f24"))) |>
  select(-c(s_value)) |>
  mutate(model = factor(model,
                        levels = c("M1: unadjusted","M2: adjusted no bmi","M3: adjusted + bmi"))) |>
  arrange(infl,model)


## Table -----

#Results Table
d_mi|>
  group_by(model) |>
  mutate(p.adj = p.adjust(p_value,method = "fdr")) |>
  # filter(infl %in% sig_markers) |>
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(estimate = paste0(estimate," (",conf_low,", ",conf_high,"), PFDR = ",p.adj)) |>
  select(-c(starts_with("predicted"),conf_low,conf_high,std_error,statistic,p_value,df,p.adj)) |>
  arrange(infl) |>
  pivot_wider(names_from = model,values_from = estimate) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


## Plot ----

#Comparison Results Plot
p_mi <- 
  d_mi |>
  filter(infl %in% sig_markers) |>
  mutate(mi = "Multiple Imputation") |>
  rename(std.error = std_error,
         p.value = p_value,
         conf.low = conf_low,
         conf.high = conf_high) |>
  bind_rows(
    
    d_reg_result |>
      filter(infl %in% sig_markers) |>
      mutate(infl = fct_reorder(.f = infl,.x = estimate,.fun = max)) |>
      mutate(model_name = case_when(model_name == "m1" ~ "M1: unadjusted",
                                    model_name == "m2" ~ "M2: adjusted no bmi",
                                    model_name == "m3" ~ "M3: adjusted + bmi")) |>
      mutate(model_name = factor(model_name,levels = c("M1: unadjusted",
                                                       "M2: adjusted no bmi",
                                                       "M3: adjusted + bmi"))) |>
    select(-c(s.value)) |>
    rename(model = model_name) |>
      mutate(mi = "Complete Case")
    ) |>
  
  mutate(mi = factor(mi,levels = c("Complete Case","Multiple Imputation"))) |>
  arrange(infl) |>
  select(-df) |>
  mutate(mm = interaction(mi,model,sep = '_')) |>
  
  ggplot(aes(y = infl, x = estimate, xmin = conf.low, xmax = conf.high,
             colour = model,
             group = mm))+
  geom_point(position = position_dodge(width = 0.35),size = 2) +
  geom_linerange(position = position_dodge(width = 0.35), linewidth = 0.8) +
  geom_vline(xintercept = 0,lty = 2) +
  scale_color_brewer(palette = "Set1")+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "gray95"),
        axis.text.x.bottom = element_text(angle = 0),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) +
  labs(y = "Immunometabolic Variable", 
       x = "Change in standardised measure") +
  facet_wrap(~mi,ncol = 2)


ggsave(filename = "./Figures/f24_lineplots_mi.pdf",
       height = 3,
       width  = 7,
       plot = p_mi)
