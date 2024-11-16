
#Introduction ------

#What I want to do, inspired by https://www.nature.com/articles/s41593-023-01404-6
#is to count up the number of extreme values for inflammatory markers our
#depressing participants have

#Load packages
pacman::p_load(tidyverse,
               tidymodels,
               marginaleffects,
               lme4,
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

source('./Final/0-Common_data.R')


covars = c("sex","bmi24","smk24","audit_c")

# Prepare Data ----

#Load up our pre-split data so we can focus on the training data only
d_split <- 
  read_rds(paste(data_dir,"alspac_data_split.rds",sep = '//'))

d_train <- 
  d_split |>
  training()

d_covars = 
  d_train |>
  select(id, all_of(covars))


# Fit normative model ----

#Start by making a normative model for the non-depressed participants
d_res <- 
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
  
  left_join(d_train |> select(id,cisr_dep), by = "id") |>
  nest_by(infl)


## Make plots ----

#Plot the distribution of residuals for an example
d_res |>
  filter(infl == "cxcl10_f24") |>
  unnest(data) |>
  ggplot(aes(resid_z,fill = cisr_dep)) +
  geom_density(alpha = 0.4)

#So the models look pretty reasonable and the stanardised residuals look good
#Although the sign appears to be inverted in a way that is confusing?


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
       height = 7,
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
  
  ggplot(aes(x = sum_ev,y = p,colour = cisr_dep)) +
  geom_point() +
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
m_res <- 
  d_res |>
  unnest(data) |>
  ungroup() |>
  group_by(id) |>
  reframe(sum_ev = sum(ex_dev))|>
  
  #Append covariates
  left_join(d_train |> select(id,cisr_dep, all_of(covars)), by = "id") |>
  
  #Model
  glm(data = _,
      formula = sum_ev ~ cisr_dep + sex + bmi24 + smk24 + audit_c ,
      family = poisson(link = "log")) 


m_res |>
  marginaleffects::avg_comparisons(model = _,
                                   variables = list(cisr_dep = 0:1)) |>
  as_tibble()

#So we get an increase in extreme values per person of 0.417

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


p_res_ev_combined <- 
  (p_res_ev/p_res_ev_m) + plot_layout(heights = c(1,1))

p_res_ev_combined

#And save
ggsave(filename = "./Figures/f24_ev.pdf",
       height = 5,
       width  = 3,
       plot = p_res_ev_combined)


# Add in cluster ID -----

#Load up the cluster IDs

clust_infl_names <-  read_csv('./Models/f24_infl_variable_clusters.csv')

## Model ------

#But this time with clusters

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

#So the variance is > the mean, so it is probably overdispersed. LMER can't do 
#quasi-distributions though so we stick with poisson
d_res_clus <- 
  d_res_clus |>
  mutate(clus = fct_relevel(clus,"2"))


#Model as before
m_res_clus2 <- glmer(s ~ cisr_dep * (sex + bmi24 + smk24 + audit_c) + cisr_dep * clus + (1|id),data = d_res_clus,
                     family = poisson(link = "log"))

#Overall model tests
car::Anova(m_res_clus2,type = 3)


#Marginal differences
marginaleffects::avg_comparisons(m_res_clus2,variables = list(cisr_dep = 0:1), by = "clus",re.form=NA)

#So after adjusting for our covariates, interesting we find we only have a 
#significant result for cluster 3, the high inflammation-depression association
#cluster, with a significant interaction overall


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
  scale_fill_manual(values = c("#32324B","#E1AF64"))+
  scale_color_manual(values = c("#32324B","#E1AF64"))+
  labs(y = "Average Count of Extreme Values", x = "Variable Cluster")


p_res_ev_combined <- 
  (p_res_ev/p_res_ev_m/p_cluster_ev) + 
  plot_layout(heights = c(1,1,1)) +
  plot_annotation(tag_levels = "A")

p_res_ev_combined

#And save
ggsave(filename = "./Figures/f24_ev.pdf",
       height = 5,
       width  = 3,
       plot = p_res_ev_combined)
