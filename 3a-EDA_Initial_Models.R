# Introduction --------

#Sometimes its just not very nice having to make markdown documents. Lets 
#reproduce the alspac-olink-basics.qmd document here but we will save the figures
#in nice ways, and we might work out how to export the tables nicely too

pacman::p_load(broom,
               tidyverse,
               tidymodels,
               rsample,
               brglm2,
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

#Set our covariates
covars = c("sex","bmi24","smk24","audit_c","ph")


# Load our datasets - invoke a script that loads from the secure RDSF folder
source('./Final/0-Common_data.R')

#Load up our pre-prepared data
d <- 
  read_rds(paste(data_dir,"alspac_data_final.rds",sep = '//'))


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
  select(id,cisr_dep,all_of(covars),ph_type) |>
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

#Physical health problems by depression
d_24_cv |>  
  drop_na(ph) |> 
  group_by(cisr_dep) |>  
  summarise(ph_1 = sum(ph == 1),
            ph_0 = sum(ph == 0)) |>
  mutate(ph_prop = ph_1/(ph_0+ph_1))


ph_t <- 
  d_24_cv |>  
  drop_na(ph) |> 
  group_by(cisr_dep) |>  
  summarise(ph_1 = sum(ph == 1),
            ph_0 = sum(ph == 0)) |>
  select(-cisr_dep) |>
  as.matrix() |>
  as.table()

dimnames(ph_t) <- 
  list(cisr_dep = c("no_dep","dep"),
       phy_health = c("ph_prob","no_prob"))


(Xsq <- stats::chisq.test(ph_t))  # Prints test summary


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

#This is combined with the results, below, to make supplementary_table_1.xlsx
d |>
  select(id,cisr_dep,all_of(covars),contains("_f24")) |>
  
  #Un-log some logged variables
  mutate(across(.cols = c(alt_f24,ast_f24,crp_f24,insulin_f24,trig_f24,vldl_f24), .fns = exp)) |>
  
  pivot_longer(contains("_f24"),names_to = "infl",values_to = "infl_value") |>
  mutate(infl = str_remove(infl,'_f24')) |>
  summarise(mu = mean(infl_value,na.rm = T),
            sigma = sd(infl_value,na.rm = T),
            .by = c("infl","sex")) |>
  arrange(infl) |>
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(mu_sig = paste0(mu," (",sigma,")")) |>
  select(-c(mu,sigma)) |>
  pivot_wider(names_from = sex,values_from = mu_sig) |>
  left_join(  v_0 |>
                as_tibble() |>
                select(name,lab) |>
                filter(str_detect(name,"_F24")) |>
                mutate(name = tolower(str_remove(name,"_F24"))) |>
                
                rename(infl = name),
              by = "infl") |>

  relocate(lab,infl,Male,Female) |>
  arrange(infl) |>
  rename(`Immunometabolic Variable` = lab,
         `Abbreviated Name` = infl) |> 
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


#Plot physical health problem variable

#Lets also plot our depression variable
p_ph = 
  d_24_cv |>
  drop_na(ph) |>
  mutate(ph = case_when(ph == 0 ~ "No",
                        ph == 1 ~ "Yes")) |>
  count(sex,ph)  |>
  ggplot(aes(x = ph, y = n,fill = sex)) +
  geom_col(colour = "black",position = position_dodge()) +
  theme(panel.grid = element_blank(),
        title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_fill_brewer(palette = 2) +
  labs(x = "Self-Reported Physical Health Problem", y = "n @ Age 24",title = "Physical Health Conditions in Participants")

#Assemble 
p_f24_covars = 
  
  (p_sex|p_cisr|p_cisr_prop)/(p_bmi24|p_drinking|(p_smoking/p_ph))+ 
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect",
              heights = c(1,2))

#And save
ggsave(filename = "./Figures/f24_covars.pdf",
       height = 7,
       width  = 8,
       plot = p_f24_covars)

#This will need a little wrangling in inkscape to get everything nicely aligned


#F24 Models --------

## Prepare Data -----

#Now we fit our models, first as a complete case analysis with the inflammatory
#variable as the IV

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
  nest_by(infl) |>
  mutate(nd = nrow(data))

#So the n for these models varies between 2592 and 2973
d_reg |> 
  ungroup() |>
  reframe(r = range(nd))


## Fit models ------

d_reg = 
  d_reg |>
  mutate(m1 = glm(formula = cisr_dep ~ infl_value,
                  data = data,
                  family = binomial(link = "logit"),
                  method = "brglmFit") |>
           list(),
         m3 = glm(cisr_dep ~ infl_value * (sex + bmi24 + smk24 + audit_c + ph),
                   data = data,
                   family = binomial(link = "logit"),
                   method = "brglmFit") |>
           list())


d_reg = 
  d_reg |>
  select(-c(data)) |>
  pivot_longer(-c(infl,nd),names_to = "model_name",values_to = "model") |>
  rowwise() |>
  mutate(me = avg_comparisons(model = model,
                              variables = list(infl_value = "sd"),
                              by = TRUE,
                              type = "response")  |>
           as_tibble() |>
           as.matrix() |>
           as_tibble() |>
           mutate(across(.cols = -c(term:contrast), .fns = ~as.double(.))) |>
           list()
                  )



## Multiple-comparison adjustment -----

#Have a look at our p value distributions
library(swfdr)

d_reg_p <- 
  d_reg |> 
  select(-c(model)) |> 
  unnest(cols = me) |> 
  filter(model_name == "m1") |>
  mutate(p.adj = p.adjust(p.value,method = "fdr"))

d_reg_p |>
  ggplot(aes(p.value)) +
  geom_histogram(binwidth = 0.05)

#So there is an increase in p values around the 0 bin which is encouraging


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

#Do the swfdr method including N as co-variate given it varies between dataset
olink_lm_qvalue <- swfdr::lm_qvalue(d_reg_p$p.value, X=d_reg_p[, c("nd")])

#Glue on these new values
d_reg_p <-
  d_reg_p |>
  mutate(sw_qv = olink_lm_qvalue$qvalues) 

d_reg_p |>
  filter(sw_qv < 0.05)


## Table of results -----

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
  select(-c(nd,model)) |>
  unnest(cols = me) |>
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(estimate = paste0(estimate," [",conf.low,", ",conf.high,"]")) |>
  arrange(infl) |>
  select(-c(std.error:predicted)) |>
  mutate(model_name = case_when(model_name == "m1" ~ "Model 1: Unadjusted",
                                model_name == "m3" ~ "Model 3: Adjusted")) |>
  pivot_wider(names_from = model_name,values_from = estimate) |>
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
  select(-c(nd,model)) |>
  unnest(cols = me) |>
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(estimate = paste0(estimate," [",conf.low,", ",conf.high,"]")) |>
  arrange(infl) |>
  select(-c(std.error:predicted)) |>
  mutate(model_name = case_when(model_name == "m1" ~ "Model 1: Unadjusted",
                                model_name == "m3" ~ "Model 3: Adjusted (including BMI)")) |>
  pivot_wider(names_from = model_name,values_from = estimate) |>
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


## Plot ----

p_reg <- 
  d_reg |>
  ungroup() |>
  select(-c(nd,model)) |>
  unnest(cols = me) |>
  filter(infl %in% sig_markers) |>
  mutate(model_name = case_when(model_name == "m1" ~ "Unadjusted",
                                model_name == "m3" ~ "Adjusted")) |>
  mutate(model_name = factor(model_name,levels = c("Unadjusted","Adjusted"))) |>

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

# Sensitivity Analyses ------

#We need to add a couple of additional sensitivity analyses which have been suggested by reviewers:

# - Some more on the haemoglobin

## Hb -----

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
      glm(cisr_dep ~ hb_f24  * (audit_c + smk24 + bmi24),
         family = binomial(link = "logit"),
         method = "brglmFit",
         data = data)
      ),
    me = list(
      avg_comparisons(model = m3,
                      variables = list(hb_f24 = "sd"),
                      by = TRUE,
                      type = "response") |>
        as_tibble()
    ))

m_hb_strat |>
  select(sex,me) |>
  unnest(cols = me) |>
  ungroup()

#So within sex groups you see no association between Hb and depression

## ALT ------

#ALT and Hb appear to be correlated e.g.

d |>
  select(id,cisr_dep,all_of(covars),contains("_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
 
   mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  
  #Sort out variable names
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  
  lm(hb ~ il6 + sex, data = _) |>
  marginaleffects::avg_comparisons() |>
  tidy()


d |>
  select(id,cisr_dep,all_of(covars),contains("_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
 # mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  
  #Sort out variable names
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  
  # ggplot(aes(x = hb, y = alt)) +
  # geom_point() +
  # geom_smooth(method = "gam",
  #             formula = y ~ s(x,k=5))+
  # facet_wrap(~sex)
  # 
  # lm(alt ~ hb *  sex, data = _) |>
  mgcv::gam(alt ~ sex + s(hb, bs = "tp",k = 4,by = sex), data = _) |>
  marginaleffects::plot_predictions(condition = list(hb = 90:180,"sex"))



#So more Hb is associated with slightly higher liver enzymes. If we add Hb in a model including other covariates, what happens?


d |>
  select(id,cisr_dep,all_of(covars),contains("_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  
  mutate(across(.cols = c(bmi24,contains("_f24")),.fns = zscore)) |>
  
  #Sort out variable names
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  
  glm(cisr_dep ~ alt * (sex + bmi24 + smk24 + audit_c + ph),
      data = _,
      family = binomial(link = "logit"),
      method = "brglmFit") |>
  avg_comparisons(model = _,
                  variables = list(alt = "sd"),
                  by = TRUE,
                  type = "response")  |>
  as_tibble() 

## Physical Health Status ------

#We are now including the physical health variable as a covariate


## Multiple Imputation Results -------

#See 3b and 3c-MI_Modelling.R for the code that produces these results; it takes
#about 2 days of computation to run so we don't run it here.
d_dep_tab <-
  read_csv("./Models/alspac_imputed_results_table.csv")


d_mi <- 
  d_dep_tab |>
  rename(infl = vars) |>
  mutate(infl = map_chr(infl,~str_remove(.x,"_f24"))) |>
  mutate(model = map_chr(model,~str_remove(.x," +"))) |>
  select(-c(s_value)) |>
  # mutate(model = factor(model,
  #                       levels = c("M1: unadjusted","M2: adjusted + bmi"))) |>
  arrange(infl,model)


### Table -----

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


### Plot ----

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
      select(-c(model)) |>
      unnest(cols = me) |>
      filter(infl %in% sig_markers) |>
      rename(model = model_name) |>
      mutate(mi = "Complete Case")
    ) |>
  select(-c(nd,term,contrast,s.value)) |>
  mutate(model = case_when(model == "m1" ~ "M1: unadjusted",
                                model == "m3" ~ "M3: adjusted + bmi")) |>
  mutate(model = factor(model,levels = c("M1: unadjusted",
                                                   "M3: adjusted + bmi"))) |>
  
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
