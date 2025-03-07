
#Introduction ------

#What I want to do, inspired by https://www.nature.com/articles/s41593-023-01404-6
#is to count up the number of extreme values for inflammatory markers our
#depressing participants have

#Load packages
pacman::p_load(tidyverse,
               tidymodels,
               lme4,
               marginaleffects,
               ggdist,
               scico,
               ComplexHeatmap,
               rsample,
               dbscan,
               patchwork,
               circlize)

tidymodels_prefer()

conflicted::conflict_prefer('select','dplyr')
conflicted::conflict_prefer('filter','dplyr')
conflicted::conflict_prefer('map'   ,'purrr')
conflicted::conflicts_prefer(crayon::`%+%`)

covars = c("sex","bmi24","smk24","audit_c","ph")

# Prepare Data ----

#Load up our pre-prepared data
source('./Final/0-Common_data.R')

d <- 
  read_rds(paste(data_dir,"alspac_data_final.rds",sep = '//'))

d_covars = 
  d_cv |>
  select(id, all_of(covars))

d_ev = 
  d |>
  select(id,cisr_dep,contains('_f24')) |>
  left_join(d_covars,by = "id")


# Fit normative model ----

#Start by making a normative model for the non-depressed participants
d_res <- 
  left_join(
    d_ev |>
      drop_na(cisr_dep) |>
      filter(cisr_dep == 0) |>
      select(-cisr_dep) |>
      pivot_longer(-c(id,all_of(covars)),names_to = 'infl',values_to = 'infl_value') |>
      drop_na() |>
      nest_by(infl,.key = "data_nd"),
    
    d_ev |>
      drop_na(cisr_dep) |>
      filter(cisr_dep == 1) |>
      select(-cisr_dep) |>
      pivot_longer(-c(id,all_of(covars)),names_to = 'infl',values_to = 'infl_value') |>
      drop_na() |>
      nest_by(infl,.key = "data_d"),
    by = "infl") |>
  
  #Fit a model predicting the inflammatory value from the covariates
  mutate(m1 = 
           list(
             lm(infl_value ~ sex + bmi24 + smk24 + audit_c + ph,data = data_nd)
           ),
         
  #Then work out the studentized residual for each individual inflammatory value         
         r1 = 
           list(rstudent(m1))  ) |>
  
  #Predict in the non-depressed group and calculate residual
  mutate(
    pred_nd = 
      list(
        predict.lm(m1,newdata = data_nd) |>
          as_tibble() |>
          rename(pred = value) |>
          bind_cols(data_nd) |>
          mutate(resid = infl_value-pred)
      ),
    
    #Do the same for the depressed individuals
    pred_d = 
      list(
        predict.lm(m1,newdata = data_d) |>
          as_tibble() |>
          rename(pred = value)|>
          bind_cols(data_d) |>
          mutate(resid = infl_value-pred)
      ),
    #Work out the normative mean and SD for each variable
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
  
  #Mark a measurement as an extreme value if it is > 2.6
  mutate(ex_dev = if_else(abs(resid_z) > 2.6,1,0) |> as.integer())  |>
  select(infl,id,all_of(covars),resid_z,ex_dev) |>
  
  left_join(d |> select(id,cisr_dep), by = "id") |>
  nest_by(infl)


## Make plots ----

#Plot the distribution of residuals for an example variable
d_res |>
  filter(infl == "il10ra_f24") |>
  unnest(data) |>
  ggplot(aes(resid_z,fill = cisr_dep)) +
  geom_density(alpha = 0.4)

#So the models look pretty reasonable and the standardized residuals look good


#Plot extreme values per variable
d_res_ev_count_variable <- 
  d_res |>
  unnest(data) |>
  ungroup() |>
  mutate(cisr_dep = case_when(cisr_dep == 0 ~ "No Depression",
                              cisr_dep == 1 ~ "ICD-10 Depression")) |>
  
  #Make an indicator if its a positive or negative
  mutate(ex_dev_sign = case_when(resid_z > 2.6  ~ 'pos',
                                 resid_z < -2.6 ~ 'neg',
                                 TRUE ~ 'none')) |>
  select(id,cisr_dep,infl,ex_dev,ex_dev_sign) |>
  reframe(s = sum(ex_dev,na.rm = T),.by = c('infl','ex_dev_sign')) |>
  filter(ex_dev_sign != "none") |>
  arrange(infl)

d_res_ev_count_variable <- 
  d_res_ev_count_variable|>
  expand(infl,ex_dev_sign) |>
  left_join(d_res_ev_count_variable,by = c("infl","ex_dev_sign")) |>
  mutate(s = if_else(is.na(s),0,s))


#How many had no positive deviations?
d_res_ev_count_variable |> filter(ex_dev_sign == "pos" & s == 0) |> dim()

#How many had no negative deviations?
d_res_ev_count_variable |> filter(ex_dev_sign == "neg" & s == 0) |> dim()

#How many had negative deviations?
d_res_ev_count_variable |> filter(ex_dev_sign == "neg" & s != 0) |> dim()

#Range
d_res_ev_count_variable |>
  reframe(range_ev = range(s))

d_res_ev_count_variable |> slice_max(s)
d_res_ev_count_variable |> slice_min(s) |> print(n = 30)


#Pos-Neg difference
d_res_ev_count_variable |>
  pivot_wider(names_from = ex_dev_sign,values_from = s) |>
  mutate(diff = if_else((pos - neg) > 0,1,0)) |>
  reframe(s_d = sum(diff))



#Plot
p_res_ev_variable <- 
  d_res_ev_count_variable|>
  mutate(ex_dev_sign = case_when(ex_dev_sign == "pos" ~ "Positive",
                                 ex_dev_sign == "neg" ~ "Negative"),
         infl = str_remove( infl,"_f24"),
         infl = fct_reorder( infl,s,max)) |>
  ggplot(aes(y =  infl, x = s,
             color = ex_dev_sign,
             fill = ex_dev_sign)) +
  geom_point() +
  theme(title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_manual(values = c("skyblue3","red3"), name = "Direction")+
  scale_color_manual(values = c("skyblue3","red3"), name = "Direction")+
  labs(y = "Immunometabolic Variable", x = "Count")

p_res_ev_variable


ggsave(filename = "./Figures/f24_ev_variable_count.pdf",
       height = 8,
       width  = 4,
       plot = p_res_ev_variable)


#Plot total number of extreme values person as previously 
p_res_ev <- 
  d_res |>
  unnest(data) |>
  ungroup() |>
  group_by(id,cisr_dep) |>
  reframe(sum_ev = sum(ex_dev)) |>
  mutate(cisr_dep = case_when(cisr_dep == 0 ~ "No Depression",
                              cisr_dep == 1 ~ "ICD-10 Depression")) |>
  group_by(cisr_dep) |>
  count(sum_ev) |>
  mutate(t = sum(n)) |>
  mutate(p = n / t) |>
  
  ggplot(aes(x = sum_ev,y = p,colour = cisr_dep,fill = cisr_dep)) +
  geom_col() +
  facet_wrap(~cisr_dep,ncol = 1)+
  # geom_point() +
  # geom_segment(aes(y = 0, yend = p, x = sum_ev,xend = sum_ev)) +
  theme(panel.grid = element_blank(),
        title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  scale_fill_manual(values = c("#32324B","#E1AF64"))+
  scale_color_manual(values = c("#32324B","#E1AF64"))+
  labs(x = "Number of extreme values", y = "Proportion of Participants")

p_res_ev


## Model ------

#Now we need to sum up the total count of extreme values for each person
d_res |>
  unnest(data) |>
  ungroup() |>
  group_by(id) |>
  reframe(sum_ev = sum(ex_dev))|>
  
  left_join(d_ev |> select(id,cisr_dep, all_of(covars)), by = "id") |>
  summarise(mu = mean(sum_ev),
            sd = sd(sum_ev),
            med = median(sum_ev),
            iq_25 = quantile(sum_ev,probs = 0.25),
            iq_75 = quantile(sum_ev,probs = 0.75),
            .by = "cisr_dep")


m_res <- 
  d_res |>
  unnest(data) |>
  ungroup() |>
  group_by(id) |>
  reframe(sum_ev = sum(ex_dev))|>
  
  left_join(d_ev |> select(id,cisr_dep, all_of(covars)), by = "id") |>
  
  glm(data = _,
      formula = sum_ev ~ cisr_dep * (sex + bmi24 + smk24 + audit_c + ph ),
      family = poisson(link = "log")) 

#Marginal difference
m_res |>
  # marginaleffects::avg_predictions(model = _, 
  #                                  by = c("cisr_dep"),
  #                                  type = "response")
  marginaleffects::avg_comparisons(model = _,
                                   variables = list(cisr_dep = 0:1)) |>
  as_tibble()

#Predicted averages
m_res |>
  marginaleffects::avg_predictions(model = _,
                                   by = c("cisr_dep"),
                                   type = "response")|>
  as_tibble()


#So we get an increase in extreme values per person of 0.214

#The residualisation is a bit like the weighting approach, and we can do a kind of 
#dobule assurance


## Summary Plot -----

p_res_ev_m <- 
  marginaleffects::plot_predictions(m_res,by = "cisr_dep",type = "response") +
  coord_flip() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "none") +
  scale_fill_manual(values = c("#32324B","#E1AF64"))+
  scale_color_manual(values = c("#32324B","#E1AF64"))+
  labs(y = "Average Count of Extreme Values", x = NULL) +
  scale_x_discrete(labels=c("0" = "No Depression","1" =  "ICD-10 Depression"))



# marginaleffects::plot_predictions(m2,condition = list("cisr_dep","sex")) +
#   coord_flip()

p_res_ev_combined <- 
  (p_res_ev|(p_res_ev_m/plot_spacer())) + plot_layout(heights = c(1,1))

p_res_ev_combined

#And save
ggsave(filename = "./Figures/f24_ev.pdf",
       height = 5,
       width  = 3,
       plot = p_res_ev_combined)


# EV by variable cluster -----

#Load up the cluster IDs
clust_infl_names <-  read_csv('./Models/f24_infl_variable_clusters.csv')

#Make the EV measure as before, but this time by cluster

#Attach these to our residualised extreme value data
d_res_clus <- 
  d_res |>
  
  #Glue on cluster ID
  mutate(infl = str_remove(infl,"_f24")) |>
  left_join(clust_infl_names |> select(infl,clus),by = "infl") |>
  mutate(clus = factor(clus)) |>
  unnest(cols = data) |>
  group_by(id,clus) |>
  reframe(s = sum(ex_dev,na.rm = T)) |>
  arrange(-s) |>
  left_join(d_covars,by = "id") |>
  left_join(d_24 |> select(id,cisr_dep),by = "id")

#What do the distributions look like for this new variable?
d_res_clus |>
  ggplot(aes(s,fill = factor(cisr_dep))) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~clus,nrow = 1)

d_res_clus |>
  group_by(factor(cisr_dep)) |>
  summarise(sd = sd(s),
            v  = var(s),
            mu = mean(s))


## Model -----

#So the variance is > the mean, so it is probably overdispersed. LMER can't do 
#quasi-distributions though so we stick with poisson

library(lme4)

d_res_clus <- 
  d_res_clus |>
  mutate(clus = fct_relevel(clus,"2"))


#Use our covariates using the same formula as before but add in the interaction of cluster id with depression
#and a random intercept for ID. This model takes a while to fit.
m_res_clus2 <- glmer(s ~ cisr_dep * (sex + bmi24 + smk24 + audit_c + ph) + cisr_dep * clus + (1|id),data = d_res_clus,
                     family = poisson(link = "log"))

#Overall model tests to get the interaction stat
car::Anova(m_res_clus2,type = 3)

#Marginal differences
marginaleffects::avg_comparisons(m_res_clus2,variables = list(cisr_dep = 0:1), by = "clus",re.form=NA) |>
  as_tibble() 

#So after adjusting for our covariates, interesting we find we only have a 
#significant result for cluster 3, the high inflammation-depression association
#cluster


## Plot -----

p_cluster_ev <- 
  marginaleffects::plot_predictions(m_res_clus2, by = c("clus","cisr_dep"))+
  coord_flip() +
  theme(#panel.grid.major.y = element_blank(),
    #panel.grid.minor.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    title = element_text(size = 8),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_blank()) +
  scale_fill_manual(values = c("#E1AF64","#32324B"))+
  scale_color_manual(values = c("#E1AF64","#32324B"))+
  labs(y = "Average Count of Extreme Values", x = "Variable Cluster")


p_res_ev_combined <- 
  (p_res_ev/p_res_ev_m/p_cluster_ev) + 
  plot_layout(heights = c(2,1,1)) +
  plot_annotation(tag_levels = "A")

p_res_ev_combined

#And save
ggsave(filename = "./Figures/f24_ev.pdf",
       height = 5,
       width  = 3,
       plot = p_res_ev_combined)



# EV by Symptom -----

#Another interesting question might be to ask if the count of extreme values
#is associated with different symptoms (your hypothesis would be is would be
#more associated with those somatic symptoms than psychological ones)
scales_use = c("ach","ftg" ,"con","slp_dec","slp_inc","app_dec","app_inc","irt",
               "sad","morn","dsx","rstl","slow"   ,"anh"    ,"sneg"   ,"mot"    ,"dcog",
               "stb","wor","anx","pho","pan")



d_res_symptom <- 
  d_res |>
  
  #Glue on cluster ID
  mutate(infl = str_remove(infl,"_f24")) |>

  unnest(cols = data) |>
  group_by(id) |>
  reframe(s = sum(ex_dev,na.rm = T)) |>
  arrange(-s) |>
  left_join(d_covars,by = "id") |>
  left_join(d |> select(id,all_of(scales_use)),by = "id") |>
  pivot_longer(all_of(scales_use),names_to = "symptom",values_to = "symptom_score") |>
  mutate(symptom_score = as.integer(symptom_score) -1) |>
  nest_by(symptom)


#Fit models
d_res_symptom <- 
  d_res_symptom |>
  mutate(adjusted = 
           list(
             glmer(s ~ symptom_score * (sex + bmi24 + smk24 + audit_c + ph)  + (1|id),
                   data = data ,
                   family = poisson(link = "log"))
           )
  )

#And its the slope for symptom score that we want
d_res_symptom <- 
  d_res_symptom |>
  ungroup() |>
  mutate(as = map(adjusted,~marginaleffects::avg_slopes(.x, variables = "symptom_score") |>
                    as_tibble()))

#Table
d_res_symptom |>
  select(-c(data, adjusted)) |>
  unnest(as) |>
  select(-c(starts_with("predicted"))) |>
  mutate(p.adj = p.adjust(p.value,method = 'fdr')) |>
  arrange(symptom) |>
  #filter(p.adj < 0.05) |>
  select(-c(term,contrast,s.value)) |>
  mutate(across(where(is.double),~janitor::round_half_up(.x,digits = 3))) |>
  mutate(estimate = paste0(estimate," (",conf.low,", ",conf.high,"), p = ",p.adj)) |>
  select(-c(conf.low,conf.high,std.error,statistic,p.value,p.adj)) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)

#Plot
p_ev_symptom <- 
  d_res_symptom |>
  select(-adjusted) |>
  unnest(as) |>

  ggplot(aes(x = estimate,xmin = conf.low,xmax = conf.high,
             y = symptom)) +
  geom_pointrange(position = position_dodge(width = 0.25)) +
  geom_vline(xintercept = 0,lty = 2) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.position = "none") +
  labs(x = "Average Change in EV Count for 1 Unit Symptom Change", y = "Symptom")

p_ev_symptom

#Save it
ggsave(filename = "./Figures/f24_cisr_symptom_ev_plot.pdf",
       height = 3,
       width  = 3,
       plot = p_ev_symptom)

#Again we could normalise this but then the axis becomes harder to understand

