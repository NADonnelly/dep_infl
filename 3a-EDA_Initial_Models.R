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

#Load up our pre-prepared data
d <- 
  read_rds(paste(data_dir,"alspac_data_final.rds",sep = '//'))

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

##



#Make a pseudo "complete case" n for tabulating based on having both a CISR value and at 
#least one immunometabolic variable
d_cc <- 
  d |> 
  select(id,c(contains("_f24"),cisr_dep)) |> 
  mutate(across(.cols = -c(id,cisr_dep),.fns = ~if_else(is.na(.x),1,0))) |>
  mutate(sum = rowSums(across(contains("_f24")))) |>
  select(id,cisr_dep,sum) |>
  filter(sum < 118) |>
  drop_na(cisr_dep) 



# F24 Covariates ------

d_24_cv = 
  d_cv |>
  left_join(d |> select(id,cisr_dep),
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

#How many have depression?
d_24_cv |> count(cisr_dep)

#And express that as a percentage
d_24_cv |> count(cisr_dep) |> mutate(t = cumsum(n),perc_d = (n/t)*100)

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

#We should probably make some kind of nice looking table 1
d_24_cv %>% 
  select(-id) %>%
  tbl_summary(by = c(sex)) %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**sex**") %>%
  bold_labels()

#There you go

#And we have our model DAG too, which is found in 0-Model_DAG.R

## Blood marker table ------

wdn |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


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

#Make a correlation plot of the F24 bloods
cor_f24 = 
  d |>
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

#We are no longer including this as a supplementary figure but please enjoy it as a little treat

#F24 Models --------

## Depression Group Difference -----

#Now we fit our models, first as a complete case analysis with the inflammatory
#variable as the DV, so we are calculating group differences between depressed
#and non-depressed
d_reg <- 
  d |>
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


### Fit models ------
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
                    as_tibble() |>
                    as.matrix() |>
                    as_tibble() |>
                    mutate(across(.cols = -c(term:contrast), .fns = ~as.double(.)))
                  )
         )

### Multiple-comparison adjustment -----

#Have a look at our p value distributions
d_reg_p <- 
  d_reg |> 
  select(-c(data,model)) |> 
  unnest(cols = me) |> 
  filter(model_name == "m1") |>
  mutate(p.adj = p.adjust(p.value,method = "BH"))

d_reg_p |>
  ggplot(aes(p.value)) +
  geom_histogram(binwidth = 0.05)


#Make q values
D0q <- qvalue::qvalue(d_reg_p$p.value)


#So what does q values think? Half a chance that the null hypothesis is true, which is fair
D0q$pi0

summary(D0q)


#Glue on the q values to our dataset
d_reg_p <-
  d_reg_p |>
  mutate(qv = D0q$qvalues) 

d_reg_p |>
  filter(qv < 0.05)

#Do the swfdr method including N as covaraite given it varies between dataset
olink_lm_qvalue <- lm_qvalue(d_reg_p$p.value, X=d_reg_p[, c("n")])

#Glue on these new values
d_reg_p <-
  d_reg_p |>
  mutate(sw_qv = olink_lm_qvalue$qvalues) 

d_reg_p |>
  filter(sw_qv < 0.05)


### Table of results -----

#Quickly tabulate all our model 1 results
d_reg_result <- 
  d_reg_p |>
  arrange(p.value) 



#Get the top markers for the base model
sig_markers = 
  d_reg_result |>
  filter(model_name == "m1") |>
  filter(sw_qv < 0.05) |>
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
d_reg |>
  ungroup() |>
  filter(infl %in% sig_markers) |>
  select(-c(data,n,model)) |>
  unnest(cols = me) |>
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(estimate = paste0(estimate," [",conf.low,", ",conf.high,"]")) |>
  arrange(infl) |>
  select(-c(std.error:predicted)) |>
  mutate(model_name = case_when(model_name == "m1" ~ "Model 1: Unadjusted",
                                model_name == "m2" ~ "Model 2: Adjusted (no BMI)",
                                model_name == "m3" ~ "Model 3: Adjusted (including BMI)")) |>
  pivot_wider(names_from = model_name,values_from = estimate) |>
  relocate(`Model 3: Adjusted (including BMI)`,.after = `Model 2: Adjusted (no BMI)`) |>
  left_join(sig_desc,by = "infl") |>
  rename(`Immunometabolic Variable` = infl) |> 
  relocate(lab,.before =  `Immunometabolic Variable`) |>
  rename(`Description` = lab) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)

#Note that the model coefficients represent the marginal change in the 
#inflammatory variable of interest between non-depressed and depressed 
#participants (in SD units)

#Tabulate all associations
d_reg |>
  ungroup() |>
  select(-c(data,n,model)) |>
  unnest(cols = me) |>
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(estimate = paste0(estimate," [",conf.low,", ",conf.high,"]")) |>
  arrange(infl) |>
  select(-c(std.error:predicted)) |>
  mutate(model_name = case_when(model_name == "m1" ~ "Model 1: Unadjusted",
                                model_name == "m2" ~ "Model 2: Adjusted (no BMI)",
                                model_name == "m3" ~ "Model 3: Adjusted (including BMI)")) |>
  pivot_wider(names_from = model_name,values_from = estimate) |>
  relocate(`Model 3: Adjusted (including BMI)`,.after = `Model 2: Adjusted (no BMI)`) |>
  left_join(  v_0 |> 
                as_tibble() |> 
                select(name,lab) |> 
                filter(str_detect(name,"_F24")) |>
                mutate(name = tolower(str_remove(name,"_F24"))) |> 

                rename(infl = name),
              by = "infl") |>
  rename(`Immunometabolic Variable` = infl) |> 
  relocate(lab,.before =  `Immunometabolic Variable`) |>
  rename(`Description` = lab) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


### Plot ----

p_reg <- 
  d_reg |>
  ungroup() |>
  select(-c(data,n,model)) |>
  unnest(cols = me) |>
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
  d |>
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

#On average women have a 5.9% higher rate of depression than men

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

#So within sex groups you see no association between Hb and depression

## Multiple Imputation Results -------


#See 3b and 3c-MI_Modelling.R for the code that produces these results; it takes
#about 2 days of computation to run so we don't run it here.
d_dep_tab <-
  read_csv("./Models/alspac_imputed_results.csv")


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
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(estimate = paste0(estimate," (",conf_low,", ",conf_high,")")) |>
  select(-c(starts_with("predicted"),conf_low,conf_high,std_error,statistic,p_value,df)) |>
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
    
    d_reg |>
      ungroup() |>
      select(-c(data,n,model)) |>
      unnest(cols = me) |>
      filter(infl %in% sig_markers) |>
      mutate(infl = fct_reorder(.f = infl,.x = estimate,.fun = max)) |>
      mutate(model_name = case_when(model_name == "m1" ~ "M1: unadjusted",
                                    model_name == "m2" ~ "M2: adjusted no bmi",
                                    model_name == "m3" ~ "M3: adjusted + bmi")) |>
      mutate(model_name = factor(model_name,levels = c("M1: unadjusted",
                                                       "M2: adjusted no bmi",
                                                       "M3: adjusted + bmi"))) |>
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
       height = 5,
       width  = 7,
       plot = p_mi)



# Experimental ------

## PCA ------

### Fit PCA Model ------

#Ask how many components to extract
n_comp <- 
  d |>
  select(contains("_f24")) |>
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  mutate(across(.cols = everything(),.fns = zscore)) |>
  drop_na() |>
  parameters::n_components(type = "PCA")

#Well apparently 28 according to parallel analysis?
as_tibble(n_comp)

#Do PCA on the variables
X = 
  d |>
  select(contains("_f24")) |>
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  mutate(across(.cols = everything(),.fns = zscore)) |>
  drop_na() |> 
  as.matrix()


# PCA on the independent variables - asking for 10 components
pca.f24 = pca(X, ncomp = 10, center = TRUE, scale = TRUE) 

#Print explained variance
pca.f24$prop_expl_var$X |> 
  as_tibble(rownames = "component") |>
  mutate(v_p = janitor::round_half_up(value*100,digits = 2))

#Plot explained variance per component
p_f24_pca_var <- 
  pca.f24$prop_expl_var$X |> 
  as_tibble(rownames = "PC") |>
  mutate(PC = str_remove(PC,"PC") |> as.integer()) |>
  ggplot(aes(x = PC,y = value)) +
  geom_col() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) +
  labs(x = "Principal Component", y = "Explained Variance") +
  scale_x_continuous(breaks = c(2,4,6,8,10))

#So it looks a lot like there is one big component and maybe another 1-3 which explain
#most of the variance in the data

# #Plot where the variables sit in PC space
# p_var_12 <-
#   plotVar(pca.f24, comp = c(1, 2), var.names = TRUE,
#         title = 'PC1 & 2')
# 
# circles <-
#   data.frame(
#   x0 = c(0,0),
#   y0 = c(0,0),
#   r  = c(0.5,1)
# )
# 
# # Behold some circles
# p_var_12_repel <-
#   p_var_12 |>
#   as_tibble() |>
#   ggplot(aes(x=x,y=y,label = names)) +
#   geom_point(colour = "red") +
#   ggrepel::geom_text_repel(box.padding = 0.5, max.overlaps = Inf)+
#   theme(panel.background = element_rect(fill = "grey95")) +
#   coord_equal(xlim = c(-1,1), ylim = c(-1,1)) +
#   ggforce::geom_circle(aes(x0 = x0,y0 = y0,r = r,
#                            x = NULL,y=NULL,label = NULL),data = circles) +
#   geom_vline(xintercept = 0,lty = 2) +
#   geom_hline(yintercept = 0, lty = 2) +
#   labs(x = "PC1", y ="PC2")
# 
# 
# p_var_23 <-
#   plotVar(pca.f24, comp = c(2, 3), var.names = TRUE,
#           title = 'PC2 & 3')
# 
# 
# p_var_23_repel <-
#   p_var_23 |>
#   as_tibble() |>
#   ggplot(aes(x=x,y=y,label = names)) +
#   geom_point(colour = "red") +
#   ggrepel::geom_text_repel(box.padding = 0.5, max.overlaps = Inf)+
#   theme(panel.background = element_rect(fill = "grey95")) +
#   coord_equal(xlim = c(-1,1), ylim = c(-1,1)) +
#   ggforce::geom_circle(aes(x0 = x0,y0 = y0,r = r,
#                            x = NULL,y=NULL,label = NULL),data = circles) +
#   geom_vline(xintercept = 0,lty = 2) +
#   geom_hline(yintercept = 0, lty = 2)+
#   labs(x = "PC2", y ="PC3")
# 
# p_var_12_repel|p_var_23_repel


#This looks like a total mess sadly and there are no really clear clusters
#within variables in  PC space

### Plot =------

#Plot scatter plots of the PCA results for the first 2 components, coloured by 
#phenotype (ICD10 depression) using the mixedomics package
p_pca_pairs_12 <- 
  plotIndiv(pca.f24,  
            comp = c(1,2),
            ind.names = FALSE, # plot the samples projected
            legend = TRUE, 
            pch=20,
            cex = 0.3,
            alpha = 0.3,
            col = "grey30") # onto the PCA subspace

p_pca_pairs_12 <- 
  p_pca_pairs_12$graph +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        panel.background = element_rect(fill = "grey95"),
        legend.position = "none",
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank())


p_pca_pairs_13 <- 
  plotIndiv(pca.f24,  
            comp = c(1,3),
            ind.names = FALSE, # plot the samples projected
            legend = TRUE, 
            pch=20,
            cex = 0.3,
            alpha = 0.3,
            col = "grey30") # onto the PCA subspace

p_pca_pairs_13 <- 
  p_pca_pairs_13$graph +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        panel.background = element_rect(fill = "grey95"),
        legend.position = "none",
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank())


#Loadings
p_pca_loading <- 
  pca.f24$loadings$X |>
  as_tibble(rownames = "variable") |>
  janitor::clean_names() |>
  dplyr::select(variable:pc4) |>
  
  mutate(variable = factor(variable,levels = vars_ordered)) |>
  arrange(variable) |>
  
  #What would be nifty would then be do order them in the same way as the 
  #correlation matrix above...
  pivot_longer(-variable) %>%
  mutate(pc = map_dbl(name,~str_remove(.x,"pc") |> as.double())) |>
  ggplot(aes(x = pc,y = variable,fill = value)) +
  geom_tile() +
  scico::scale_fill_scico(palette = 'vik') +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title.y.left = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.title.x.bottom = element_text(size = 8) ,
        axis.text = element_text(size = 6),
        axis.text.x.bottom = element_text(angle = 90)) +
  coord_fixed() +
  scale_y_discrete(position = "right") +
  labs(x = "Principal Component", y = NULL)

p_corr_loading <- 
  (p_f24_cor+p_pca_loading) + plot_layout(widths = c(15,1))


p_pca_combo = 
  (p_corr_loading/(p_f24_pca_var|p_pca_pairs_12|p_pca_pairs_13|plot_spacer())) + 
  plot_layout(heights = c(4,1)) +
  plot_annotation(tag_levels = "A")

#Save our plots
ggsave(filename = "./Figures/f24_pca_combined.pdf",
       height = 8,
       width  = 10,
       plot = p_pca_combo)

### PCA Regression ------

#What if we took our PCA scores and used them in regression, as a little test
#of what they relate to

pc_scores <- 
  pca.f24$variates |> 
  as.data.frame() |> 
  as_tibble() |> 
  rename_with(.cols = everything(),.fn = ~str_remove(.x,'X.')) |> 
  janitor::clean_names()

pc_iv <- 
  d |>
  select(id,sex,cisr_dep,smk24,bmi24,audit_c,contains("_f24")) |>
  mutate(across(.cols = contains("_f24"),.fns = zscore)) |>
  drop_na(contains("_f24")) |>
  select(-contains("_f24"))

pc_reg <- 
  bind_cols(pc_iv,pc_scores) |>
  select(-c(pc5:pc10)) |>
  pivot_longer(cols = starts_with("pc"),
               names_to = "pc",
               values_to = "pc_score") |>
  
  #Fit as one big model??
  nest_by(pc) |>
  mutate(
    m1 = list(
      lm(pc_score ~ cisr_dep * (sex + bmi24 + smk24 + audit_c),data = data)
    ),
    ac = list(
      avg_comparisons(m1) |> as_tibble()
    )
  )

#Plot all the marginal effects

p_pca_loading_wide <- 
  pca.f24$loadings$X |>
  as_tibble(rownames = "variable") |>
  janitor::clean_names() |>
  dplyr::select(variable:pc4) |>
  
  mutate(variable = factor(variable,levels = vars_ordered)) |>
  arrange(variable) |>
  
  #What would be nifty would then be do order them in the same way as the 
  #correlation matrix above...
  pivot_longer(-variable) %>%
  mutate(pc = map_dbl(name,~str_remove(.x,"pc") |> as.double())) |>
  ggplot(aes(y = pc,x = variable,fill = value)) +
  geom_tile() +
  scico::scale_fill_scico(palette = 'vik') +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title.y.left = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.title.x.bottom = element_text(size = 8) ,
        axis.text = element_text(size = 8),
        axis.text.x.bottom = element_text(angle = 90)) 


p_pc_reg <- 
  pc_reg |> 
  select(ac) |>
  unnest(cols = c('ac')) |> 
  mutate(term = fct_relevel(term,"cisr_dep")) |>
  arrange(term) |>
  
  ggplot(aes(x = term,y = estimate,ymin = conf.low,ymax = conf.high)) + 
  geom_pointrange() +
  facet_wrap(~pc,ncol = 4) +
  geom_hline(yintercept = 0,lty = 2) +
  theme(panel.grid = element_blank(),
        # panel.background = element_blank(),
        axis.title.y.left = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.title.x.bottom = element_text(size = 8) ,
        axis.text = element_text(size = 6),
        axis.text.x.bottom = element_text(angle = 90))

p_pca_loading_wide/p_pc_reg

p_pca_reg = 
  (p_pca_loading_wide/p_pc_reg) + 
  plot_layout(heights = c(1,1)) +
  plot_annotation(tag_levels = "A")

p_pca_reg

#Save our plots
ggsave(filename = "./Figures/f24_pca_reg.pdf",
       height = 4,
       width  = 8,
       plot = p_pca_reg)

#Or as a heat plot I guess
pc_reg |> 
  select(ac) |>
  unnest(cols = c('ac')) |> 
  mutate(term = fct_relevel(term,"cisr_dep")) |>
  ggplot(aes(x = pc,y = term,fill = -log10(p.value))) +
  geom_tile() +
  scico::scale_fill_scico(palette = 'devon')

#Look at the regression output, adjusting for multiple comparisons
pc_reg |>
  mutate(dep_me = list(
    avg_comparisons(model = m1,
                    variables = list(cisr_dep = 0:1),conf_level = 0.95) |>
      as_tibble() 
  )) |>
  select(pc,dep_me) |>
  unnest(cols = 'dep_me') |>
  ungroup() |>
  mutate(p_adj = p.adjust(p.value,method = 'fdr'))

#So none of the top 4 PCs are significantly associated with depression. By 
#contrast, we get lots of associations with covariates, especially sex, bmi and
#smoking status. 

## Split correlation plot ------

#Alternative version - top half people without depression, bottom half
#people with depression. Will need to do the clustering slightly differently
#to ensure consistency though.....

cor_no_dep <- 
  d |>
  filter(cisr_dep == 0) |> #gives you an n of 2287
  select(contains("_f24")) |>
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  mutate(across(.cols = everything(),.fns = zscore)) |>
  cor(use = "pairwise.complete.obs",
      method = "pearson")

#Do clustering on these guys
clust_no_dep <- 
  hclust(as.dist(1 - cor_no_dep),method = "complete")

ord <- clust_no_dep$order

cor_no_dep <- cor_no_dep[ord, ord]

p_f24_cor_nd = 
  ggcorrplot::ggcorrplot(cor_no_dep,
                         method = "square",
                         type = "upper",
                         show.diag = FALSE,
                         hc.order = F,
                         tl.cex = 6) +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_text(angle = 90))

#Now get the correlation matrix for the depressed people
cor_dep <- 
  d |>
  filter(cisr_dep == 1) |> #gives you an n of 2287
  select(contains("_f24")) |>
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  mutate(across(.cols = everything(),.fns = zscore)) |>
  cor(use = "pairwise.complete.obs",
      method = "pearson")

#Apply same ordering
cor_dep <- cor_dep[ord, ord]

#Plot
p_f24_cor_d = 
  ggcorrplot::ggcorrplot(cor_dep,
                         method = "square",
                         type = "upper",
                         show.diag = FALSE,
                         hc.order = F,
                         tl.cex = 6) +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_text(angle = 90))

#Assembly
p_corr_combined = p_f24_cor_nd/p_f24_cor_d


#Save our plots
ggsave(filename = "./Figures/f24_corr_combined.pdf",
       height = 16,
       width  = 8,
       plot = p_corr_combined)


#Looks very similar to be honest, and I think makes the PCA more complicated

## Other adjustments -----

#So in the CIS-R there are some self-report variables about physical health. Inflammation
#can be increased by physical health problems after all, so should we be looking at these
#variables?

#We have:
# - FKDQ2500 - Number of times visited the GP (0 - 4 scale), presumably could include visits for MH
# - FKDQ2510 - List of conditions (but note MH included)
# - FKDQ2520 - longstanding illness,disability or infirmity (but presumably could include MH)

#So lets focus on FKDQ2510 and select the non-MH physical conditions (codes 1 - 7)

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

#Take CCL25 as an example
d <- d_reg$data[[11]]

#Lets grab FKDQ2510
d_ph <- 
  d |>
  left_join(d_0 |> select(alnqlet,FKDQ2510) |> rename(id = alnqlet),
            by = 'id') |>
  mutate(ph = case_when(FKDQ2510 %in% c(1:7) ~ 1,
                        FKDQ2510 == 8|FKDQ2510 == 9 ~ 0,
                        FKDQ2510 < 1 ~ NA_integer_)) |>
  select(-FKDQ2510)

d_ph |>
  lm(formula = infl_value ~ cisr_dep * (sex + audit_c + smk24 + bmi24 + ph), 
     data = _) |>
  avg_comparisons(model = _,
                  variables = list(cisr_dep = 0:1),conf_level = 0.95) |>
  tidy()


d_ph_adj <- 
  d_train |>
  select(id,cisr_dep,all_of(covars),contains("_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  
  #Sort out variable names
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  
  left_join(d_0 |> select(alnqlet,FKDQ2510) |> rename(id = alnqlet),
            by = 'id') |>
  mutate(ph = case_when(FKDQ2510 %in% c(1:7) ~ 1,
                        FKDQ2510 == 8|FKDQ2510 == 9 ~ 0,
                        FKDQ2510 < 1 ~ NA_integer_)) |>
  select(-FKDQ2510) |>
  
  #Put each blood variable into its own dataset
  pivot_longer(-c(id,cisr_dep,all_of(covars),ph), names_to = "infl", values_to = "infl_value") |> 
  drop_na() |>
  nest_by(infl) |>
  
  mutate(m1 = lm(formula = infl_value ~ cisr_dep, 
                 data = data) |>
           list(),
         m2 = lm(formula = infl_value ~ cisr_dep * (sex + audit_c + smk24), 
                 data = data) |>
           list(),
         m3 = lm(formula = infl_value ~ cisr_dep * (sex + audit_c + smk24 + bmi24), 
                 data = data) |>
           list(),
         m4 = lm(formula = infl_value ~ cisr_dep * (sex + audit_c + smk24 + bmi24 + ph), 
                 data = data) |>
           list())  |>
  
  select(-c(data)) |>
  pivot_longer(-infl,names_to = "model_name",values_to = "model") |>
  mutate(me = map(model,~avg_comparisons(model = .x,
                                         variables = list(cisr_dep = 0:1),conf_level = 0.95) |>
                    as_tibble() |>
                    as.matrix()))

#Plot?
d_ph_adj |>
  mutate(me_tab = map(me, ~as_tibble(.x) |> 
                        select(estimate:conf.high) |>
                        mutate(across(.cols = everything(),.fns = ~as.double(.x))))) |>
  select(-c(model:me)) |>
  unnest(me_tab) |>
  filter(infl %in% c("il6","wbc","cubcp","alt","ccl25","hgf","fgf19","fgf21")) |>
  ggplot(aes(x = estimate,xmin = conf.low,xmax = conf.high,
             fill = model_name,colour = model_name,y = infl)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbarh(position = position_dodge(width = 0.2)) +
  geom_vline(xintercept = 0,lty = 2)

#Quickly tabulate our model results
d_ph_adj |>
  mutate(me_tab = map(me, ~as_tibble(.x) |> 
                        select(estimate:conf.high) |>
                        mutate(across(.cols = everything(),.fns = ~as.double(.x))))) |>
  select(-c(model:me)) |>
  unnest(me_tab) |>
  arrange(p.value) |>
  group_by(model_name)|>
  filter(infl %nin% c("rbc","hct")) |>
  mutate(p.adj = p.adjust(p.value,method = 'fdr')) |>
  ungroup() |>
  filter(p.adj < 0.05) |>
  arrange(model_name)

#So adjusting for the self-reported physical health conditions does not seem to make a big difference.

#Another way to adjust might be to just exclude everyone with a raw CRP > 10

d_crp_adj <- 
  d_train |>
  select(id,cisr_dep,all_of(covars),contains("_f24")) |>
  
  #Sort out variable names
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  
  filter(crp < 10) |>
  
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  
  left_join(d_0 |> select(alnqlet,FKDQ2510) |> rename(id = alnqlet),
            by = 'id') |>
  mutate(ph = case_when(FKDQ2510 %in% c(1:7) ~ 1,
                        FKDQ2510 == 8|FKDQ2510 == 9 ~ 0,
                        FKDQ2510 < 1 ~ NA_integer_)) |>
  select(-FKDQ2510) |>
  
  #Put each blood variable into its own dataset
  pivot_longer(-c(id,cisr_dep,all_of(covars),ph), names_to = "infl", values_to = "infl_value") |> 
  drop_na() |>
  nest_by(infl) |>
  
  mutate(m1 = lm(formula = infl_value ~ cisr_dep, 
                 data = data) |>
           list(),
         m2 = lm(formula = infl_value ~ cisr_dep * (sex + audit_c + smk24), 
                 data = data) |>
           list(),
         m3 = lm(formula = infl_value ~ cisr_dep * (sex + audit_c + smk24 + bmi24), 
                 data = data) |>
           list(),
         m4 = lm(formula = infl_value ~ cisr_dep * (sex + audit_c + smk24 + bmi24 + ph), 
                 data = data) |>
           list())  |>
  
  select(-c(data)) |>
  pivot_longer(-infl,names_to = "model_name",values_to = "model") |>
  mutate(me = map(model,~avg_comparisons(model = .x,
                                         variables = list(cisr_dep = 0:1),conf_level = 0.95) |>
                    as_tibble() |>
                    as.matrix()))

#Plot?
d_crp_adj |>
  mutate(me_tab = map(me, ~as_tibble(.x) |> 
                        select(estimate:conf.high) |>
                        mutate(across(.cols = everything(),.fns = ~as.double(.x))))) |>
  select(-c(model:me)) |>
  unnest(me_tab) |>
  filter(infl %in% c("il6","wbc","neutrophils","cdcp1","alt","ccl25","hgf","fgf19","fgf21","insulin")) |>
  ggplot(aes(x = estimate,xmin = conf.low,xmax = conf.high,
             fill = model_name,colour = model_name,y = infl)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbarh(position = position_dodge(width = 0.2)) +
  geom_vline(xintercept = 0,lty = 2)

#Table?
d_crp_adj |>
  mutate(me_tab = map(me, ~as_tibble(.x) |> 
                        select(estimate:conf.high) |>
                        mutate(across(.cols = everything(),.fns = ~as.double(.x))))) |>
  select(-c(model:me)) |>
  unnest(me_tab) |>
  arrange(p.value) |>
  group_by(model_name)|>
  filter(infl %nin% c("rbc","hct")) |>
  mutate(p.adj = p.adjust(p.value,method = 'fdr')) |>
  ungroup() |>
  filter(p.adj < 0.05) |>
  arrange(model_name)

#The results are broadly the same, after correction for multiple comparisons we are left with the same set of 
#variables

## Batch ID ------

#Should we include Olink batch (as a random intercept?)
library(lme4)


d_batch <- 
  d_train |>
  select(id,cisr_dep,all_of(covars),contains("_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  
  #Sort out variable names
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  
  #Add batch ID
  left_join(d_0 |> select(alnqlet,Batch_F24) |> rename(id = alnqlet),
            by = "id") |>
  
  #Put each blood variable into its own dataset
  pivot_longer(-c(id,cisr_dep,all_of(covars),Batch_F24), names_to = "infl", values_to = "infl_value") |> 
  drop_na() |>
  mutate(Batch_F24 = factor(Batch_F24)) |>
  nest_by(infl)

#So the n for these models varies between
d_batch |> mutate(nd = nrow(data) |> list()) |> ungroup() |> reframe(r = range(nd))

d_batch = 
  d_batch |>
  mutate(m1 = lmer(formula = infl_value ~ cisr_dep + (1|Batch_F24), 
                   data = data) |>
           list(),
         m2 = lmer(formula = infl_value ~ cisr_dep * (sex + audit_c + smk24) + (1|Batch_F24), 
                   data = data) |>
           list(),
         m3 = lmer(formula = infl_value ~ cisr_dep * (sex + audit_c + smk24 + bmi24) + (1|Batch_F24), 
                   data = data) |>
           list()) 

d_batch = 
  d_batch |>
  select(-c(data)) |>
  pivot_longer(-infl,names_to = "model_name",values_to = "model") |>
  mutate(me = map(model,~avg_comparisons(model = .x,
                                         variables = list(cisr_dep = 0:1),
                                         conf_level = 0.95,
                                         re.form=NA) |>
                    as_tibble() |>
                    as.matrix() |>
                    as_tibble() |>
                    mutate(across(.cols = -c(term:contrast), .fns = ~as.double(.)))
  )
  )

#Quickly tabulate our model results
d_batch_result <- 
  d_batch |>
  mutate(me_tab = map(me, ~as_tibble(.x) |> 
                        select(estimate:conf.high) |>
                        mutate(across(.cols = everything(),.fns = ~as.double(.x))))) |>
  select(-c(model:me)) |>
  unnest(me_tab) |>
  arrange(p.value) |>
  group_by(model_name)|>
  filter(infl %nin% c("rbc","hct")) |>
  mutate(p.adj = p.adjust(p.value,method = 'fdr')) |>
  ungroup()

d_batch_result |>
  filter(p.adj<0.05)

#So the results are the same as with the lm first time around

## Frequentist Matching ------

#Out of interest, again for a supplement, using matching
d_reg <- 
  d |>
  select(id,cisr_dep,all_of(covars),contains("_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  
  #Sort out variable names
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  
  #Put each blood variable into its own dataset
  pivot_longer(-c(id,cisr_dep,all_of(covars)), names_to = "infl", values_to = "infl_value") |> 
  drop_na() |>
  nest_by(infl)

#Take CCL25 as an example
d1 <- d_reg$data[[13]]

d_match1 <-
  MatchIt::matchit(cisr_dep ~ bmi24 + sex + smk24 + audit_c,
                   data = d1,
                   method = "full",
                   estimand = "ATE")


plot(d_match1,type = "jitter",interactive = FALSE)

#Black is treated, grey is control
plot(d_match1,type = "density",interactive = FALSE,
     which.xs = ~sex + smk24 + audit_c + bmi24)

#Do loveplot here
cobalt::love.plot(d_match1,stars = "std")



#So things get better after matching which is nice

#Now we fit our model
d_match_fit <-
  d_match1 |>
  MatchIt::match.data() |>
  lm(infl_value ~ cisr_dep * (bmi24 + sex + smk24 + audit_c),
     data = _,
     weights = weights)

#calculate marginal effects
d_match_fit |>
  avg_comparisons(variables = list(cisr_dep = 0:1),
                  vcov = "HC3",
                  #newdata = subset(d_match1 |> match.data(),cisr_dep == 1), #ATT
                  newdata = d_match1 |> MatchIt::match.data(),                      #ATE
                  wts = "weights") |>
  as_tibble()


#Apply to each blood variable
d_match <- 
  d_reg |>
  mutate(match_data = list(
    MatchIt::matchit(cisr_dep ~ bmi24 + sex + smk24 + audit_c,
                     data = data,
                     method = "full",
                     estimand = "ATE") |>
      MatchIt::match.data()
  ))|>
  mutate(model =  list(
    lm(infl_value ~ cisr_dep * (bmi24 + sex + smk24 + audit_c),
       data = match_data,
       weights = weights))
  ) |>
  mutate(me = list(
    avg_comparisons(model = model,
                    variables = list(cisr_dep = 0:1),
                    vcov = "HC3",
                    #newdata = subset(match_data,cisr_dep == 1), #ATT
                    newdata = match_data,                      #ATE
                    wts = "weights") |>
      as_tibble()
  ))

#### Tabulate Results -------

#Inspect output
d_match |>
  select(infl,me) |>
  unnest(me) |>
  arrange(p.value)


## GAMLSS ------

#What if we model whether depression is also associated with increased variablity
#in our immunometabolic variables?
library(gamlss)

d_gamlss <- 
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

t <-  d_gamlss$data[[11]]

g_1 <- gamlss(data   = t,
              formula = infl_value ~ cisr_dep * (sex + bmi24 + smk24 + audit_c),
              sigma.formula = ~1,
              family  = NO(mu.link = "identity", sigma.link = "log")) 

marginaleffects::avg_comparisons(g_1,variables = list(cisr_dep = 0:1),
                                 what = 'mu')

marginaleffects::avg_predictions(g_1,by = "cisr_dep",
                                 what = 'mu')

#Now do the same thing but also with a model on sigma
g_2 <- gamlss(data   = t,
              formula = infl_value ~ cisr_dep * (sex + bmi24 + smk24 + audit_c),
              sigma.formula = ~ cisr_dep * (sex + bmi24 + smk24 + audit_c),
              family  = NO(mu.link = "identity", sigma.link = "log")) 

#Compare means
marginaleffects::avg_comparisons(g_2,variables = list(cisr_dep = 0:1),
                                 what = 'mu')

marginaleffects::avg_predictions(g_2,by = "cisr_dep",
                                 what = 'mu')

#Compare sigmas
marginaleffects::avg_comparisons(g_2,variables = list(cisr_dep = 0:1),
                                 what = 'sigma')

marginaleffects::avg_predictions(g_2,by = "cisr_dep",
                                 what = 'sigma')

#So sigmas don't differ, but mus do in this example.

#Now we can deploy this to all the immunometabolic variables. I recall we have to do this as
#a loop as gamlss gets upset inside map

d_gamlss <- 
  d_gamlss |>
  mutate(model = vector(mode = "list",length = 1),
         me    = vector(mode = "list",length = 1))

for(i in 1:dim(d_gamlss)[1]){
  
  t <-  d_gamlss$data[[i]]
  
  #Now do the same thing but also with a model on sigma
  g_t <- gamlss(data   = t,
                formula = infl_value ~ cisr_dep * (sex + bmi24 + smk24 + audit_c),
                sigma.formula = ~ cisr_dep * (sex + bmi24 + smk24 + audit_c),
                family  = NO(mu.link = "identity", sigma.link = "log")) 
  
  
  # marginaleffects::avg_predictions(g_t,variables = list(cisr_dep = 0:1),
  #                                  what = 'sigma',
  #                                  type = "response") |>
  #   as_tibble() |>
  #   mutate(parameter = "sigma")|>
  #   ggplot(aes(x = cisr_dep,y = estimate,ymin = conf.low,ymax = conf.high)) +
  #   geom_pointrange()
  # 
  # 
  # marginaleffects::avg_predictions(g_t,variables = list(cisr_dep = 0:1),
  #                                  what = 'mu',
  #                                  type = "response") |>
  #   as_tibble() |>
  #   mutate(parameter = "mu") |>
  #   ggplot(aes(x = cisr_dep,y = estimate,ymin = conf.low,ymax = conf.high)) +
  #   geom_pointrange()
  
  g_comp <- 
    bind_rows(
      #Compare means
      marginaleffects::avg_comparisons(g_t,variables = list(cisr_dep = 0:1),
                                       newdata = t,
                                       what = 'mu',
                                       type = "response") |>
        as_tibble() |>
        mutate(parameter = "mu")
      ,
      #Compare sigmas
      marginaleffects::avg_comparisons(g_t,variables = list(cisr_dep = 0:1),
                                       what = 'sigma',
                                       type = "response") |>
        as_tibble() |>
        mutate(parameter = "sigma")
    ) |>
    as.matrix() |>
    as_tibble() |>
    mutate(across(.cols = estimate:conf.high,.fns = ~janitor::round_half_up(as.double(.x),digits = 10)))
  
  
  d_gamlss$model[[i]] <- g_t
  d_gamlss$me[[i]] <- g_comp
  
  # marginaleffects::avg_predictions(g_t,by = "cisr_dep",
  #                                  what = 'mu')
  # 
  # marginaleffects::avg_predictions(g_t,by = "cisr_dep",
  #                                  what = 'sigma')
  
  
}

#See what we got - for the mu values
d_gamlss |>
  select(infl,me) |>
  unnest(me) |>
  ungroup() |>
  filter(parameter == "mu") |>
  arrange(p.value) |>
  mutate(p.adj = p.adjust(p.value,method = 'fdr')) |> 
  filter(p.value < 0.05)

#Nothing withstands correction for multiple comparisons here, but its the same 
#set of values you would expect

d_gamlss |>
  select(infl,me) |>
  unnest(me) |>
  ungroup() |>
  filter(parameter == "sigma") |>
  arrange(p.value) |>
  mutate(p.adj = p.adjust(p.value,method = 'fdr'))|> 
  # filter(p.value < 0.05) |>
  filter(p.adj < 0.05) |>
  print(n = 93)

#So apparently there are 11 variables with difference in sigma; they are 
#different to the variables with different means



## Depression as DV ------

library(brglm2)

# If one goal in immunopsychiatry research is to delineate the nature of the 
#relationships between immune mediators and mental health outcomes, we ought
#to flip our models around and ask whether the probability of depression
#variables across the distribution of our immune mediators, and whether this is
#linear. We may have to do logistic regression here, for which I apologise

d_slope <- 
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

m_t <- 
  d_slope |>
  filter(infl == "il6") |>
  pluck('data',1) |>
  glm(cisr_dep ~ infl_value * (sex + bmi24 + smk24 + audit_c),
      data = _,
      family = binomial(link = "probit"),
      method = "brglm.fit")

marginaleffects::plot_predictions(m_t,
                                  condition = list("infl_value" = c(-3:3),
                                                   "sex"),
                                  vcov = "HC3",
                                  type = "response")


m_t2 <- 
  d_slope |>
  filter(infl == "il6") |>
  pluck('data',1) |>
  glm(cisr_dep ~ ns(infl_value,3) * (sex + bmi24 + smk24 + audit_c),
      data = _,
      family = binomial(link = "probit"))


marginaleffects::plot_predictions(m_t2,
                                  condition = list("infl_value" = c(-3:3),
                                                   "sex"),
                                  vcov = "HC3",
                                  type = "response")


## Non-linear relationships -----

library(mgcv)

g1 <- 
  d_slope |>
  filter(infl == "il6") |>
  pluck('data',1) |>
  gam(formula = cisr_dep ~ infl_value + sex + bmi24 + smk24 + audit_c,
      data = _,
      family = binomial(link = "logit"))

g2 <- 
  d_slope |>
  filter(infl == "il6") |>
  pluck('data',1) |>
  gam(formula = cisr_dep ~ s(infl_value) + sex + bmi24 + smk24 + audit_c,
      data = _,
      family = binomial(link = "logit"),
      method = "REML")

g3 <- 
  d_slope |>
  filter(infl == "il6") |>
  pluck('data',1) |>
  gam(formula = cisr_dep ~ s(infl_value,by = sex) + sex + bmi24 + smk24 + audit_c,
      data = _,
      family = binomial(link = "logit"),
      method = "REML")

g4 <- 
  d_slope |>
  filter(infl == "il6") |>
  pluck('data',1) |>
  gam(formula = cisr_dep ~ s(infl_value,by = sex) + sex + s(bmi24) + smk24 + audit_c,
      data = _,
      family = binomial(link = "logit"),
      method = "REML")

#Try to make a fully-specified model
g5 <- 
  d_slope |>
  filter(infl == "il6") |>
  pluck('data',1) |>
  gam(formula = cisr_dep ~ s(infl_value,by = sex) + sex + s(bmi24) + te(infl_value,bmi24) + s(infl_value,by = smk24) + smk24 + s(audit_c) + te(infl_value,audit_c),
      data = _,
      family = binomial(link = "logit"),
      method = "REML")

AIC(g1,g2,g3,g4,g5)

#Notably the massive model does not fit better, although its slower

p_g1 <- marginaleffects::plot_predictions(g1,
                                          condition = list("infl_value" = seq(-3,3,0.2),
                                                           "sex"),
                                          type = "response")

p_g2 <- marginaleffects::plot_predictions(g2,
                                          condition = list("infl_value" = seq(-3,3,0.2),
                                                           "sex"),
                                          type = "response")

p_g3 <- marginaleffects::plot_predictions(g3,
                                          condition = list("infl_value" = seq(-3,3,0.2),
                                                           "sex"),
                                          type = "response")

p_g5 <- marginaleffects::plot_predictions(g5,
                                          condition = list("infl_value" = seq(-3,3,0.2),
                                                           "sex"),
                                          type = "response")


p_g1|p_g2|p_g3|p_g5

#You can also plot the slope of that line, although I do find it harder to interpret
marginaleffects::plot_slopes(g2,
                             variables = "infl_value",
                             condition = list("infl_value" = seq(-3,3,0.2),
                                              "sex"),
                             type = "response")

#So by AIC the spline model is not a much better fit than the standard model

#But this is just one blood measure; what about the other measures. We fit only
#our simple model with linear trends for the blood measure, and one with a 
#spline fit for males and female separately
d_slope <- 
  d_slope |>
  mutate(
    m_linear = 
      list(
        gam(formula = cisr_dep ~ infl_value * sex + bmi24 + smk24 + audit_c,
            data = data,
            family = binomial(link = "logit"))
      ),
    m_spline =
      list(
        gam(formula = cisr_dep ~ s(infl_value,by = sex) + sex + bmi24 + smk24 + audit_c,
            data = data,
            family = binomial(link = "logit"),
            method = "REML")
        
      )
  )

#Now lets look at some AICs
d_slope_aic <- 
  d_slope |>
  select(-data) |>
  pivot_longer(-infl) |>
  mutate(aic = map_dbl(value,AIC)) |>
  select(-value) |>
  pivot_wider(names_from = name,values_from = aic) |>
  mutate(aic_diff = m_spline - m_linear)

d_slope_aic|>
  ggplot(aes(aic_diff)) +
  geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 0,lty = 2) +
  theme(panel.grid = element_blank())

d_slope_aic |>
  summarise(mu = mean(aic_diff),
            sd = sd(aic_diff),
            med = median(aic_diff),
            mad = mad(aic_diff),
            min = min(aic_diff),
            max = max(aic_diff))

#So the AIC is basically similar by median, but with some exceptions in the 
#negative direction, which might be interesting to pursue - 

#Lets look at some extreme examples
d_slope_aic |>
  slice_min(aic_diff,n = 1)

#cst5
d_slope |>
  filter(infl == "cst5") |>
  pluck('m_linear',1) |>
  marginaleffects::plot_predictions(condition = list("infl_value" = seq(-3,3,0.2),
                                                     "sex"),
                                    type = "response")

d_slope |>
  filter(infl == "cst5") |>
  pluck('m_spline',1) |>
  marginaleffects::plot_predictions(condition = list("infl_value" = seq(-3,3,0.2),
                                                     "sex"),
                                    type = "response")

#Weird U shaped thing going on in spline world

d_slope_aic |>
  slice_max(aic_diff,n = 1)

#ggt
d_slope |>
  filter(infl == "ggt") |>
  pluck('m_linear',1) |>
  marginaleffects::plot_predictions(condition = list("infl_value" = seq(-3,3,0.2),
                                                     "sex"),
                                    type = "response")

d_slope |>
  filter(infl == "ggt") |>
  pluck('m_spline',1) |>
  marginaleffects::plot_predictions(condition = list("infl_value" = seq(-3,3,0.2),
                                                     "sex"),
                                    type = "response")

#Not much difference


#Lets focus on the values associated with differences in depression

d_slope_plots <- 
  d_slope |>
  filter(infl %in% sig_markers) |>
  select(-data) |>
  pivot_longer(-infl) |>
  mutate(plot = map2(value,infl,~marginaleffects::plot_predictions(.x,
                                                                   condition = list("infl_value" = seq(-3,3,0.2),
                                                                                    "sex"),
                                                                   type = "response") +
                       labs(x = "Inflammatory Value",title = .y)))


wrap_plots(d_slope_plots$plot,ncol = 2)

#So there are a few things going on which are curious

d_slope |>
  select(-data) |>
  filter(infl %in% sig_markers) |>
  pivot_longer(-infl) |>
  mutate(aic = map_dbl(value,AIC)) |>
  select(-value) |>
  pivot_wider(names_from = name,values_from = aic) |>
  mutate(aic_diff = m_spline - m_linear)

#The AIC differences are only very small. The curios are-

#IL6 might have a sigmoid relationship in depression in females?
#The ALT relationship appears to only really be in males - not sure why low
#ALT would be associated with more depression in males particularly?

#Is the ALT distribution really whack? No it looks like a pretty standard
#lognormal - although levels tend to be higher in males so is there something
#particular about males with lower levels?
d_0 |> 
  select(ALT_F24,kz021) |> 
  filter(ALT_F24 > -1) |> 
  ggplot(aes(log10(ALT_F24))) + 
  geom_density() + 
  facet_wrap(~kz021,ncol = 1)


#Something is a bit odd about neutrophils and WBC count - this middle peak thing
#with a lot of increased uncertainty at larger values. Could be another sigmoid
#type relationship in males where above average neutrohpils increases risk of 
#depression, but within that top half more is not necessarily more, a bit like
#the IL6 relationship.
d_slope_plots$plot[[12]]

#Neutrophil count tends to be lower in males
d_0 |> 
  select(Neutrophils_F24,kz021) |> 
  filter(Neutrophils_F24 > -1) |> 
  ggplot(aes(log10(Neutrophils_F24))) + 
  geom_density() + 
  facet_wrap(~kz021,ncol = 1)

# I think on the basis that the AIC differences are only small, and further that 
#our ML methods should be able to capture any non-linear relationships when
#we fit them, that we can leave this here



## Frequentist Weighting ======

library(WeightIt)
library(cobalt)

weight_test = 
  d_reg |> 
  filter(infl == "il6") |> 
  pluck('data',1) |>
  mutate(sex = if_else(sex == "Male",1,0)) 


bal.tab(infl_value ~ sex + bmi24 + audit_c + smk24,
        data = weight_test, 
        thresholds = c(m = .05)) |>
  love.plot()

#So we are not balanced with respect to bmi particularly

#Start with standard weights
W.glm <- weightit(infl_value ~ sex + bmi24 + audit_c + smk24,
                  data = weight_test, 
                  method = "glm")

summary(W.glm)

love.plot(W.glm)

#This does not work at all. Try entropy balancing

W.out <- weightit(infl_value ~ sex + bmi24 + audit_c + smk24,
                  data = weight_test, 
                  moments = 2,
                  int = TRUE,
                  method = "ebal")


love.plot(W.out)
summary(W.out)

#This does work which is nice 
bal.tab(W.out, stats = c("m", "v"), thresholds = c(m = .05))

#OK so entropy balancing worked nicely...although perhaps too nicely?

#Next we can fit our model with a builtin
fit_w <- glm_weightit(cisr_dep ~ splines::ns(infl_value,df = 4) * (sex + bmi24 + audit_c + smk24),
                      family = binomial(link = 'logit'),
                      data = weight_test,
                      weightit = W.out)

avg_comparisons(fit_w, variables = "infl_value")



#Representative values of infl_value:
values <- with(weight_test, 
               seq(quantile(infl_value, .1),
                   quantile(infl_value, .9),
                   length.out = 31))

#G-computation
p <- avg_predictions(fit_w,
                     variables = list(infl_value = values))

ggplot(p, aes(x = infl_value)) +
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3) +
  labs(x = "infl_value", y = "E[Y|A]") +
  theme_bw()


# Estimate the pointwise derivatives at representative
# values of Ac
s <- avg_slopes(fit_w,
                variables = "infl_value",
                newdata = datagrid(infl_value = values,
                                   grid_type = "counterfactual"),
                by = "infl_value")

# Plot the AMEF
ggplot(s, aes(x = infl_value)) +
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "infl_value", y = "dE[Y|A]/dA") +
  theme_bw()



#Directly compare to the plain model
fit0 <- glm(cisr_dep ~ splines::ns(infl_value,df = 4) * (sex + bmi24 + audit_c + smk24),
            family = binomial(link = 'logit'),
            data = weight_test)

p_0 <- avg_predictions(fit0,
                     variables = list(infl_value = values))

ggplot(p_0, aes(x = infl_value)) +
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3) +
  labs(x = "infl_value", y = "E[Y|A]") +
  theme_bw()


#Add in brglm2 for comparison
library(brglm2)

fit_b2 <- glm(cisr_dep ~ splines::ns(infl_value,df = 4) * (sex + bmi24 + audit_c + smk24),
              family = binomial(link = 'logit'),
              method = "brglm_fit",
              data = weight_test,
              weightit = W.out$weights)

p_b2 <- avg_predictions(fit_b2,
                       variables = list(infl_value = values))

ggplot(p_b2, aes(x = infl_value)) +
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3) +
  labs(x = "infl_value", y = "E[Y|A]") +
  theme_bw()

#Sort of in the middle

#So this is interesting because the relationship appears to
#be more non-linear in the weighted version

#How do you go about summarising this kind of non-linear relationship?
#



#So lets roll this out to all our variables
d_weight = 
  d_reg |> 
  
  
  
  select(infl:data)



#Fit a regression model for each inflammatory value 
d_weight = 
  d_weight |>
  
  #Do the weights
  ungroup() |>
  mutate(w_d <- map(data,~weightit(infl_value ~ sex + bmi24 + audit_c + smk24,
                                   data = .x, 
                                   moments = 2,
                                   int = TRUE,
                                   method = "ebal")))

d_weight = 
  d_weight |>
  rename(w_d = `... <- NULL`)  |>

  #Fit a model
  mutate(model = map2(data,w_d, ~glm_weightit(cisr_dep ~ splines::ns(infl_value,df = 4) * (sex + bmi24 + audit_c + smk24),
                                              data = .x,
                                              family = binomial(link = 'logit'),
                                              weightit = .y) )) |>
  
  #calculate the marginal effect of our inflammatory variable on depression
  mutate(dep_contrast = map(model,~avg_comparisons(model = .x,
                                                   variables = "infl_value") |>
                              as_tibble()))



##Tabulate the results
d_weight_result = 
  d_weight |>
  select(-c(data,model,w_d)) |>
  ungroup() |>
  unnest(dep_contrast) |>
  mutate(p.adj = p.adjust(p.value,method = "fdr")) |>
  arrange(p.value) |>
  select(infl,term,estimate,conf.low,conf.high,std.error,statistic,p.value,p.adj) 


#Tabulate the weighted results
d_weight_result |>
  mutate(across(where(is.double),~round(.x,digits = 3))) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)



#P value adjustments as below
pacman::p_load(swfdr,qvalue)

#Just have a look at our p value distributions
d_reg_p <- 
  d_reg |> 
  select(-c(data,model)) |> 
  unnest(cols = me) |> 
  mutate(p.adj = p.adjust(p.value,method = "BH"))

d_weight_result |>
  ggplot(aes(p.value)) +
  geom_histogram(binwidth = 0.05)

#So there is a bit of something going on with our dataset, indicated by the clump of pvalues as the lower end

#Make q values
D0q <- qvalue::qvalue(d_weight_result$p.value)

#So what does q values think? Half a chance that the null hypothesis is true, which is fair
D0q$pi0

summary(D0q)

#Trad bonferroni adjustment gives us no positive results whatsoever
sum(D0q$p < 0.05/dim(d_weight_result)[1])

#Glue on the q values 
d_weight_result <-
  d_weight_result |>
  mutate(qv = D0q$qvalues) 

d_weight_result |>
  filter(qv < 0.05)

#But we can now do a p value adjustment thing based on some covariates. I wondered which covariates to use.
#In the package example they have gwas data and use the N and the allele frequency. We could also use the N
#and we could use the normative variance of the individual variable? But our variances are z-scored, so
#all the SDs are going to be basically 1...So what covariate would actually be useful? Number/proportion of extreme values?

#Lets try with with N anyway
olink_lm_qvalue <- lm_qvalue(d_weight_result$p.value, X=map_dbl(d_reg$data,~nrow(.x)))

d_weight_result <-
  d_weight_result |>
  mutate(sw_qv = olink_lm_qvalue$qvalues) 

d_weight_result |>
  filter(sw_qv < 0.05)

#Nope, that doesn't work


## Bayesian Models ------

#Experiment with doing these as rstanarm models

#Lets start with the models
d_reg = 
  d_reg |>
  mutate(b1 = map(data, ~stan_glm(formula = infl_value ~ cisr_dep,
                                  family = gaussian(link = "identity") , 
                                  data = .x,
                                  prior_intercept = rstanarm::normal(0,1.5),
                                  prior_aux = rstanarm::exponential(rate = 1),
                                  prior = rstanarm::student_t(3,0,0.5),
                                  seed = 1, refresh = 0,
                                  chains = 4, cores = 6,
                                  iter = 4000, warmup = 1000,
                                  refresh = 1,
                                  adapt_delta = 0.995)),
         b2 = map(b1, ~update(.x,formula. = infl_value ~ cisr_dep * (sex + audit_c + smk24))),
         b3 = map(b1, ~update(.x,formula. = infl_value ~ cisr_dep * (sex + audit_c + smk24 + bmi24)))
  )

#Now do the marginal effects
d_reg <- 
  d_reg |>
  pivot_longer(-infl,names_to = "model_name",values_to = "model") |>
  mutate(me = map(model,~avg_comparisons(.x,variables = list(cisr_dep = 0:1),conf_level = 0.95))) |>
  mutate(me_tab = map(me,~as_tibble(.x) |> as.matrix())) |>
  mutate(me_pd  = map_dbl(me,~posterior_draws(.x) |>
                            pull(draw) |>
                            bayestestR::p_direction() |>
                            bayestestR::pd_to_p() |>
                            as.numeric()))

#Extract our model results as a nice table
d_reg_result <- 
  d_reg |>
  mutate(me_tab = map(me_tab, ~as_tibble(.x) |> 
                        select(estimate:conf.high) |>
                        mutate(across(.cols = everything(),.fns = ~as.double(.x))))) |>
  select(-c(model:me)) |>
  unnest(me_tab) |>
  arrange(me_pd)


#### Table of results -----

#Get the top markers for the fully adjusted model
sig_markers = 
  d_reg_result |>
  filter(model_name == "b3") |>
  #arrange(me_pd) |> print(n = 30)
  #filter(infl %nin% c("rbc","hct")) |>
  filter(me_pd < 0.05) |>
  #mutate(p.adj = p.adjust(me_pd,method = 'fdr')) |>
  #filter(p.adj < 0.05) |>
  pull(infl)

#Tabulate  
d_reg_result |>
  ungroup() |>
  filter(infl %in% sig_markers) |>
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(estimate = paste0(estimate," (",conf.low,", ",conf.high,"), pd = ",me_pd)) |>
  arrange(infl) |>
  select(-c(conf.low:me_pd)) |>
  mutate(model_name = case_when(model_name == "b1" ~ "unadjusted",
                                model_name == "b2" ~ "adjusted no bmi",
                                model_name == "b3" ~ "adjusted + bmi")) |>
  pivot_wider(names_from = model_name,values_from = estimate) |>
  relocate(`adjusted + bmi`,.after = `adjusted no bmi`) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


#Note that the model coefficients represent the marginal change in mean value
#of the stnadardised inflammatory marker e.g. for IL6 thats - 
d_reg_result |> 
  filter(infl == "il6" & model_name == "b3")|>
  mutate(estimate = round_half_up(estimate, digits = 3)) |> 
  pull(estimate)

#e.g. 0.152 SD units

d_reg |>
  filter(infl == "il6" & model_name == "b3") |>
  pluck('model',1) |>
  avg_predictions(variables = list(cisr_dep = 0:1))


#But of course firstly this is very slow, and secondly doing the full set of MI
#and matching with the Bayesian models is very very slow and therefore it is 
#more tractable to stick with frequentist methods throughout the paper for 
#consistency

#### Plot ------

d_reg_post <- 
  d_reg |>
  filter(infl %in% sig_markers) |>
  select(c(infl,model_name,me)) |>
  mutate(model_name = case_when(model_name == "b1" ~ "unadjusted",
                                model_name == "b2" ~ "adjusted no bmi",
                                model_name == "b3" ~ "adjusted + bmi")) |>
  mutate(model_name = factor(model_name,levels = c("unadjusted","adjusted no bmi","adjusted + bmi"))) |>
  mutate(me_post = map(me,~posterior_draws(.x))) |>
  select(-c(me)) |>
  unnest(me_post) 

p_reg <- 
  d_reg_post|>
  ggplot(aes(y = infl, x = draw, fill = model_name,colour = model_name))+
  stat_halfeye(slab_alpha = 0.5,position = position_dodge(width = 0.75),
               point_interval = "mean_qi")+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
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

## Limma ======

pak::pak('limma')

library(limma)

t <- 
  d_train |>
  select(id,cisr_dep,all_of(covars),contains("_f24")) |>
  drop_na()

npx_df <- 
  t |>
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  select(id,contains("_f24")) |>
  
  #Sort out variable names
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24"))  



pheno_df <- 
  t |>
  select(id,cisr_dep,all_of(covars)) |>
  mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  # mutate(cisr_dep = case_when(cisr_dep == 0 ~ "no_dep",
  #                             cisr_dep == 1 ~ "dep")) |>
  as.data.frame()
  
  rownames(npx_df) <- 
  npx_df$id

npx_df <- 
  npx_df %>% 
  select(-id) |> 
  as.matrix()
  

#Here’s the bit of code I used (npx_df includes id and olink columns, pheno_df is id and everything else):
design <- model.matrix(~ 0 + cisr_dep , pheno_df)

head(design)

# dim(npx_df)
# dim(pheno_df)
# dim(design)

fit <- lmFit(t(npx_df), design)

contr <- makeContrasts(cisr_dep1-cisr_dep0, 
                       levels = design)

fit2 <- 
  fit |>
  contrasts.fit(contr) |>
  eBayes(proportion = 0.3)

sig_tab <- 
  topTable(fit2, lfc=0, n=190,p.value = 0.05,adjust.method = "BH") |>
  as_tibble(rownames = "infl")

print(sig_tab)

volcanoplot(fit2,coef = 1,highlight = dim(sig_tab)[1],names = colnames(npx_df),hl.col = "red4")


#So a limma model set up to be the same as an unadjusted model produces... exactly the same
#results as my univariate pairwise models. Cool cool cool. 

#On this basis it feels reasonable to use the same univariate pairwise modelling approach 
#throughout


## Other forms of p value adjustment -----

# Based on a comment from Ruby and this paper - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1716-1

#Can we do better than the good old BH method? It is worth noting that the paper does not suggest big differences
#and also says that differences only become apparent for some methods with thousands of predictors
# library("BiocManager")
# BiocManager::install("swfdr")

pacman::p_load(swfdr,qvalue)


#Fit models
d_reg <- 
  d |>
  select(id,cisr_dep,all_of(covars),contains("_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  
  #Sort out variable names
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  
  #Put each blood variable into its own dataset
  pivot_longer(-c(id,cisr_dep,all_of(covars)), names_to = "infl", values_to = "infl_value") |> 
  drop_na() |>
  nest_by(infl)

d_reg = 
  d_reg |> 
  mutate(n = dim(data)[1]) |>
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
  pivot_longer(-c(infl,data,n),names_to = "model_name",values_to = "model") |>
  mutate(me = map(model,~avg_comparisons(model = .x,
                                         variables = list(cisr_dep = 0:1),conf_level = 0.95) |>
                    as_tibble() |>
                    as.matrix() |>
                    as_tibble() |>
                    mutate(across(.cols = -c(term:contrast), .fns = ~as.double(.))) |>
                    select(-c(term,contrast,s.value))
  )
  )

#Just have a look at our p value distributions
d_reg_p <- 
  d_reg |> 
  select(-c(data,model)) |> 
  unnest(cols = me) |> 
  filter(model_name == "m1") |>
  mutate(p.adj = p.adjust(p.value,method = "BH"))

d_reg_p |>
  ggplot(aes(p.value)) +
  geom_histogram(binwidth = 0.05)

#So there is a bit of something going on with our dataset, indicated by the clump of pvalues as the lower end

#Make q values
D0q <- qvalue::qvalue(d_reg_p$p.value)


#So what does q values think? Half a chance that the null hypothesis is true, which is fair
D0q$pi0

summary(D0q)

#Trad bonferroni adjustment gives us no positive results whatsoever
sum(D0q$p < 0.05/dim(d_reg_p)[1])

#Glue on the q values 
d_reg_p <-
  d_reg_p |>
  mutate(qv = D0q$qvalues) 

d_reg_p |>
  filter(qv < 0.05)

#But we can now do a p value adjustment thing based on some covariates. I wondered which covariates to use.
#In the package example they have gwas data and use the N and the allele frequency. We could also use the N
#and we could use the normative variance of the individual variable? But our variances are z-scored, so
#all the SDs are going to be basically 1...So what covariate would actually be useful? Number/proportion of extreme values?

#Lets try with with N anyway
olink_lm_qvalue <- lm_qvalue(d_reg_p$p.value, X=d_reg_p[, c("n")])

d_reg_p <-
  d_reg_p |>
  mutate(sw_qv = olink_lm_qvalue$qvalues) 

d_reg_p |>
  filter(sw_qv < 0.05)

#So this gets us a few more hits but does not change the overall picture
#Maybe we try the extreme values

#### EV model -----

d_res <- 
  left_join(
    d |>
      select(id,cisr_dep,contains('_f24'),all_of(covars)) |>
      drop_na(cisr_dep) |>
      filter(cisr_dep == 0) |>
      select(-cisr_dep) |>
      pivot_longer(-c(id,all_of(covars)),names_to = 'infl',values_to = 'infl_value') |>
      drop_na() |>
      nest_by(infl,.key = "data_nd"),
    
    d |>
      select(id,cisr_dep,contains('_f24'),all_of(covars)) |>
      drop_na(cisr_dep) |>
      filter(cisr_dep == 1) |>
      select(-cisr_dep) |>
      pivot_longer(-c(id,all_of(covars)),names_to = 'infl',values_to = 'infl_value') |>
      drop_na() |>
      nest_by(infl,.key = "data_d"),
    by = "infl") |>
  
  mutate(m1 = 
           list(
             lm(infl_value ~ sex + bmi24 + smk24 + audit_c,data = data_nd)
           ),
         #lm(infl_value ~ sex,data = data_nd)
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
  unnest(cols = c(mu_nd,sd_nd)) |>
  select(infl,pred_nd,pred_d,mu_nd,sd_nd) |>
  pivot_longer(starts_with("pred_")) |>
  unnest(value) |>
  
  #Standardise the residuals using the mean and sd of residuals from the non-depressed group
  mutate(resid_z = (resid-mu_nd)/sd_nd) |>
  
  mutate(ex_dev = if_else(abs(resid_z) > 2.6,1,0) |> as.integer())  |>
  select(infl,id,all_of(covars),resid_z,ex_dev) |>
  
  left_join(d |> select(id,cisr_dep), by = "id") |>
  nest_by(infl)


d_res_ev_count_variable <- 
  d_res |>
  unnest(data) |>
  ungroup() |>
  mutate(cisr_dep = case_when(cisr_dep == 0 ~ "No Depression",
                              cisr_dep == 1 ~ "ICD-10 Depression")) |>

  select(id,cisr_dep,infl,ex_dev) |>
  reframe(s = sum(ex_dev,na.rm = T),.by = c('infl')) 


#Glue on the extreme value counts
d_reg_p <- 
  d_reg_p |>
  #select(-c(s.x,s.y)) |>
  left_join(d_res_ev_count_variable |> mutate(infl = str_remove(infl,'_f24')),by = "infl")



olink_lm_qvalue_s <- lm_qvalue(d_reg_p$p.value, X=d_reg_p[, c("n","s")])

d_reg_p <-
  d_reg_p |>
  mutate(sw_qv_s = olink_lm_qvalue_s$qvalues) 

d_reg_p |>
  filter(sw_qv_s < 0.05)

#Does not make a difference. Overall I conclude that although this approach may have
#some utility, especially for GWAS etc with thousands of tests, with only a few hundred
#tests, we are not really seeing much by way of benefit


## Everything all the time -----

library(WeightIt)
library(gamlss)

#What if we did weights, a GAMLSS and the p value adjustment thing
d_gamlss <- 
  d |>
  select(id,cisr_dep,all_of(covars),contains("_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  
  #Sort out variable names
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  
  #Put each blood variable into its own dataset
  pivot_longer(-c(id,cisr_dep,all_of(covars)), names_to = "infl", values_to = "infl_value") |> 
  drop_na() |>
  nest_by(infl)

d_gamlss <- 
  d_gamlss |>
  mutate(model = vector(mode = "list",length = 1),
         me    = vector(mode = "list",length = 1)) |>
  mutate(n     = nrow(data))

for(i in 1:dim(d_gamlss)[1]){
  
  t <-  d_gamlss$data[[i]]
  
  
  #Do a model with covariates but no interactiosn etc
  #fit_s <- 
  
    
  #Do the weightint
  
  W.out <- weightit(infl_value ~ sex + bmi24 + audit_c + smk24,
                    data = t, 
                    moments = 2,
                    int = TRUE,
                    method = "ebal")
  
  
  # love.plot(W.out)
  # summary(W.out)
  
  #This does work which is nice 
  # bal.tab(W.out, stats = c("m", "v"), thresholds = c(m = .05))
  
  #OK so entropy balancing worked nicely...although perhaps too nicely?
  
  #Next we can fit our model with a builtin
  fit_w <- glm_weightit(cisr_dep ~ infl_value * (sex + bmi24 + audit_c + smk24),
                        family = binomial(link = 'logit'),
                        data = t,
                        weightit = W.out)
  
  # avg_comparisons(fit_w,variables = list(infl_value = "sd"))
  # plot_predictions(fit_w,condition = 'infl_value')
  
  # avg_comparisons(fit_w,variables = list(infl_value = "sd"),
  #                 newdata = t,
  #                 type = "response") |>
  #   as_tibble() 
  
  
  
  #Now do the same thing but also with a model on sigma
  fit_g <- gamlss(data   = t,
                formula = cisr_dep ~ infl_value * (sex + bmi24 + audit_c + smk24),
                sigma.formula = ~ infl_value * (sex + bmi24 + audit_c + smk24),
                weights = W.out$weights,
                family  = BI) 

  # plot_predictions(fit_g,what = 'mu',condition = 'infl_value')

  g_comp <- 
    bind_rows(
      #GAMLLS
      marginaleffects::avg_comparisons(fit_g,variables = list(infl_value = "sd"),
                                       newdata = t,
                                       what = 'mu',
                                       type = "response") |>
        as_tibble() |>
        mutate(method = "gamlss"),
      
      #Weighted GLM
      marginaleffects::avg_comparisons(fit_w,variables = list(infl_value = "sd"),
                                       newdata = t,
                                       type = "response") |>
        as_tibble() |>
        mutate(method = "weightit")
    ) |>
    as.matrix() |>
    as_tibble() |>
    mutate(across(.cols = estimate:conf.high,.fns = ~janitor::round_half_up(as.double(.x),digits = 10)))
  
  
  d_gamlss$model[[i]] <- fit_g
  d_gamlss$me[[i]]    <- g_comp

}

#See what we got - for the mu values
d_gamlss_p <- 
  d_gamlss |>
  select(infl,me) |>
  unnest(me) |>
  ungroup() |>
  filter(method == "weightit") |>
  arrange(p.value) |>
  mutate(p.adj = p.adjust(p.value,method = 'fdr')) 

#Nothing withstands correction for multiple comparisons here, but its the same 
#set of values you would expect
 
d_gamlss_p |>
  ggplot(aes(p.value)) +
  geom_histogram(binwidth = 0.05)

#Flatter

#Make q values
D0q <- qvalue::qvalue(d_gamlss_p$p.value)

#So what does q values think? Half a chance that the null hypothesis is true, which is fair
D0q$pi0

summary(D0q)

#Trad bonferroni adjustment gives us no positive results whatsoever
sum(D0q$p < 0.05/dim(d_gamlss_p)[1])

#Glue on the q values 
d_gamlss_p <-
  d_gamlss_p |>
  mutate(qv = D0q$qvalues) 

d_gamlss_p |>
  filter(qv < 0.05)

#Good old IL6

#But we can now do a p value adjustment thing based on some covariates. I wondered which covariates to use.
#In the package example they have gwas data and use the N and the allele frequency. We could also use the N
#and we could use the normative variance of the individual variable? But our variances are z-scored, so
#all the SDs are going to be basically 1...So what covariate would actually be useful? Number/proportion of extreme values?

#Lets try with with N anyway
olink_lm_qvalue <- lm_qvalue(d_gamlss_p$p.value, X=d_gamlss_p$n)

d_gamlss_p <-
  d_gamlss_p |>
  mutate(sw_qv = olink_lm_qvalue$qvalues) 

d_gamlss_p |>
  filter(sw_qv < 0.05)

d_gamlss_p |>
  ggplot(aes(x = p.adj,y = sw_qv)) +
  geom_point() +
  geom_smooth(formula = y ~ x,method = 'lm') +
  geom_abline(intercept = 0,slope = 1,lty = 2) +
  coord_fixed()

#So this squashes our adjusted p values downwards


## The Psychosis -----

#What if we used psychosis as our DV

# Current psychotic experiences (the Pliks variable FKPL2010) - 1 trial binomial /bernoulli model
# Current psychotic disorder (the PLiks variable FKPL2210) - 1 trial binomial /bernoulli model
# Current at risk mental state (the presence of any of the CAARMS variables FKPL2700, FKPL2710 or FKPL2720) - 1 trial binomial /bernoulli model


d_psychosis <- 
  d_0 |>
  select(alnqlet,FKPL2010,FKPL2210,FKPL2700,FKPL2710,FKPL2720) |>
  rename(id = alnqlet) |>
  mutate(across(where(is.double),~if_else(.x < 0,NA_integer_,.x))) |>
  rename(p_exp  = FKPL2010,
         p_dis  = FKPL2210) |>
  transmute(id,p_dis,p_exp,
            p_arms = if_else((FKPL2700 + FKPL2710 + FKPL2720) > 0,1,0)) 


d_psychosis_reg <- 
  d |>
  select(id,all_of(covars),contains("_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  
  #Sort out variable names
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  
  #Put each blood variable into its own dataset
  pivot_longer(-c(id,all_of(covars)), names_to = "infl", values_to = "infl_value") |> 
  
  left_join(d_psychosis,by = "id") |>
  pivot_longer(p_dis:p_arms, names_to = "p_var", values_to = "psychosis") |>
  drop_na() |>
  nest_by(infl,p_var)


d_psychosis_reg = 
  d_psychosis_reg |>
  mutate(m1 = lm(formula = infl_value ~ psychosis, 
                 data = data) |>
           list(),
         m2 = lm(formula = infl_value ~ psychosis * (sex + audit_c + smk24), 
                 data = data) |>
           list(),
         m3 = lm(formula = infl_value ~ psychosis * (sex + audit_c + smk24 + bmi24), 
                 data = data) |>
           list()) 

d_psychosis_me = 
  d_psychosis_reg |>
  select(-c(data)) |>
  pivot_longer(-c(infl,p_var),names_to = "model_name",values_to = "model") |>
  mutate(me = map(model,~avg_comparisons(model = .x,
                                         variables = list(psychosis = 0:1),conf_level = 0.95) |>
                    as_tibble() |>
                    as.matrix()))

#Quickly tabulate our model results
d_psychosis_result <- 
  d_psychosis_me |>
  mutate(me_tab = map(me, ~as_tibble(.x) |> 
                        select(estimate:conf.high) |>
                        mutate(across(.cols = everything(),.fns = ~as.double(.x))))) |>
  select(-c(model:me)) |>
  unnest(me_tab) |>
  arrange(p.value)

d_psychosis_result |>
  nest_by(p_var,model_name) |>
  mutate(sig_vars = list(
    data |> mutate(p.adj = p.adjust(p.value,method = "fdr")) |> filter(p.adj < 0.05) |> pull(infl)
  )) |>
  select(-data) |>
  unnest(sig_vars) |>
  ungroup()

#Actually its the same stuff mostly. IL6 cannot be defeated.


## Clustering on raw Olink Data -----

library(parameters)

#Hierarchical Clustering
clust_cor_f24 <- 
  hclust(as.dist(1 - cor_f24),method = "complete")

#We can ask how many clusters to cut into thus
n <- n_clusters(as_tibble(d_t2), package = "all")
plot(n)

#Elbow
x <- n_clusters_elbow(as_tibble(d_hm))
x
plot(x)


# Gap
x <- n_clusters_gap(as_tibble(d_hm))
x
as.data.frame(x)
plot(x)

# Silhouette
x <- n_clusters_silhouette(as_tibble(d_hm))
x
as.data.frame(x)
plot(x)

#DBSCAN method
x <- n_clusters_dbscan(as_tibble(d_hm), method = "kNN", min_size = 0.05) # 5 percent
x
head(as.data.frame(x))
plot(x)

x <- n_clusters_dbscan(as_tibble(d_hm),min_size = 0.1)
x
head(as.data.frame(x))
plot(x)

#HCLUST method
x <- n_clusters_hclust(as_tibble(d_hm), iterations = 50, ci = 0.90)
x
head(as.data.frame(x), n = 10) # Print 10 first rows
plot(x)

#So 2 clusters?
clu <- cluster_analysis(as_tibble(d_hm),n = 3, method = "pam")

clu

predict(clu)

plot(clu) + theme_bw()


clu2 <- cluster_analysis(as_tibble(d_hm),
                         n = NULL,
                         method = "hclust",
                         iterations = 500,
                         ci = 0.90
)

clu2
plot(clu2) + theme_bw()

clu3 <- cluster_analysis(as_tibble(d_hm), method = "mixture")
plot(clu3) + theme_bw()

#### Metacluster -----

rez_kmeans  <- cluster_analysis(as_tibble(d_hm), n = 3, method = "kmeans")
rez_hclust  <- cluster_analysis(as_tibble(d_hm), n = 3, method = "hclust")
rez_hkmeans <- cluster_analysis(as_tibble(d_hm), n = 3, method = "hkmeans")
rez_pam     <- cluster_analysis(as_tibble(d_hm), n = 3, method = "pam")
rez_hclust2 <- cluster_analysis(as_tibble(d_hm),
                                n = NULL,
                                method = "hclust",
                                iterations = 500,
                                ci = 0.90
)
rez_dbscan  <- cluster_analysis(as_tibble(d_hm), method = "dbscan", dbscan_eps = 3.5)
rez_hdbscan <- cluster_analysis(as_tibble(d_hm), method = "hdbscan")
rez_pamk    <- cluster_analysis(as_tibble(d_hm), method = "pamk")
rez_mixture <- cluster_analysis(as_tibble(d_hm), method = "mixture")

list_of_results <- list(
  rez_kmeans, rez_hclust, rez_hkmeans, rez_pam,
  rez_hclust2, rez_dbscan, rez_hdbscan, rez_pamk,rez_mixture
)

probability_matrix <- cluster_meta(list_of_results)

# Plot the matrix as a reordered heatmap
heatmap(probability_matrix, 
        #Rowv = NA, Colv = NA,
        scale = "none",
        col = grDevices::hcl.colors(256, palette = "inferno")
)

predict(probability_matrix, n = 2) |>
  as_tibble()

d <- as.dist(abs(probability_matrix - 1))
model <- hclust(d)
plot(model, hang = -1)
