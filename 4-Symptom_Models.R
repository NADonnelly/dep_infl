# Introduction -------

#Lets do all symptom x marker models, and then apply some subscale scores based 
#previous studies. We will fit fully adjusted models, and then make heatmaps
#of the resulting marginal effects

pacman::p_load(tidyverse,
               tidymodels,
               marginaleffects,
               ggdist,
               scico,
               ComplexHeatmap,
               rsample,
               dbscan,
               patchwork,
               circlize,
               brglm2)

tidymodels_prefer()

conflicted::conflict_prefer('select','dplyr')
conflicted::conflict_prefer('filter','dplyr')
conflicted::conflict_prefer('map'   ,'purrr')
conflicted::conflicts_prefer(crayon::`%+%`)

source('./Final/0-Common_data.R')

#Load up our pre-split data so we can focus on the training data only
d_split <- 
  read_rds(paste(data_dir,"alspac_data_split.rds",sep = '//'))

#And get the training data for analysis
d_train <- 
  d_split |>
  training()

#Set heatmap colour schemes
col_fun = colorRamp2(seq(from = -0.03,to =  0.03, length.out = 11),
                     rev(scico(11, palette = 'vik',direction = -1)))

## Prepare the immuno-metabolic dataset -----

#Prepare the immunometabolic dataset
d_im = 
  d_train |>
  
  #Lets select from these variables; get rid of the CISR outcome data
  select(id, contains("_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(-c(id),zscore))

im_names = 
  d_im |>
  select(-c(id)) |>
  colnames()

##Prepare our covariates -------

#These are the covariates that our DAG suggests form the minimum adjustment set
covars = c("sex","bmi24","smk24","audit_c")

d_covars = 
  d_train |>
  select(id, all_of(covars)) 


# Diagnosis association -------

## Model =====

#Model the categorical CISR diagnosis here

d_icd = 
  left_join(d_im,
            d_train |> select(id,cisr_dep),
            by = "id") |>
  left_join(d_covars,by = "id") |>
  pivot_longer(all_of(im_names), 
               names_to = "infl", values_to = "infl_value") |> 
  drop_na() |>
  nest_by(infl) 


d_icd = 
  d_icd |>
  mutate(
    model = glm(cisr_dep ~ infl_value * (sex + bmi24 + smk24 + audit_c),
                family = binomial(link = 'logit'),
                data = data,
                method = "brglmFit") |> list(),
    me = avg_comparisons(model,variables = list(infl_value = "sd"),
                         by = TRUE,
                         type = "response") |>
      as_tibble() |>
      list())



## Heatmap -----

#Make a heatmap with clustering of all inflammation x symptom effects
d_hm_icd = 
  d_icd |>
  mutate(infl = str_remove(infl,"_f24")) |>
  select(-c(model,data)) |>
  unnest(me) |>
  ungroup() |>
  select(infl,estimate) |>
  arrange(infl) |>
  select(-infl) |>
  as.matrix()

rownames(d_hm_icd) <- 
  d_icd |>
  mutate(infl = str_remove(infl,"_f24"))  |>
  distinct(infl) |>
  pull(infl)

column_ha0 = 
  HeatmapAnnotation(bar0 = anno_barplot(d_hm_icd |> rowMeans(),ylim = c(-0.015,0.015)))

row_ha0 = 
  rowAnnotation(bar0 = anno_barplot(d_hm_icd |> colMeans(),ylim = c(-0.005,0.005)))

p_hm_icd =
  Heatmap(d_hm_icd |> t(),
          col = col_fun,
          na_col = "grey",
          top_annotation = column_ha0, 
          right_annotation = row_ha0,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          heatmap_legend_param = list(title = "Marginal Effect"))


# Symptom scores -----

## Build symptom scores ------

d_symptom_score <- 
  d_train |>
  select(id,atypical,psychological,somatic,anxiety)
  
#Look at how these synthetic scores are distributed
p_cisr_scores <- 
  d_symptom_score |>
  pivot_longer(-id) |>
  mutate(name = str_to_title(name)) |>
  ggplot(aes(value)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~name,ncol = 4) +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  theme(panel.grid = element_blank(),
        axis.title  = element_text(size = 8),
        axis.text   = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_blank()) +
  labs(x = "Score", 
       y = "Count")

#As each item has been made into a binary, we can fit binomial models 
ss_trial = tibble(ss = c("atypical","psychological","somatic","anxiety"),
                  nt = c(5,4,9,7))


## Model -----

d_ss = 
  left_join(d_im,
            d_symptom_score,
            by = "id") |>
  left_join(d_covars,by = "id") |>
  
  pivot_longer(all_of(im_names), 
               names_to = "infl", values_to = "infl_value") |> 
  pivot_longer(atypical:anxiety, 
               names_to = "ss", values_to = "s_value") |> 
  drop_na() |>
  nest_by(infl,ss) |>
  left_join(ss_trial,by = "ss")


#I think we need to stick to a common modelling approach. 
d_ss = 
  d_ss |>
  mutate(model = glm(cbind(s_value, nt - s_value) ~ infl_value * (sex + bmi24 + smk24 + audit_c),
                     family = binomial(link = 'logit'),
                     data = data,
                     method = "brglmFit") |> list(),
         me = avg_comparisons(model,variables = list(infl_value = "sd"),
                              by = TRUE,
                              type = "response") |>
           as_tibble() |>
           list()
  )


## Heatmap ------

d_hm_ss = 
  d_ss |>
  select(infl,ss,me) |>
  unnest(me) |>
  ungroup() |>
  select(infl,ss,estimate) |>
  pivot_wider(id_cols = "infl",names_from = "ss",values_from = "estimate") %>%
  arrange(infl) |>
  select(-infl) |>
  as.matrix()

rownames(d_hm_ss) <- 
  d_ss |>
  ungroup() |>
  mutate(infl = str_remove(infl,"_f24"))  |>
  distinct(infl) |>
  pull(infl)

column_ha = 
  HeatmapAnnotation(bar1 = anno_barplot(d_hm_ss |> rowMeans(),ylim = c(-0.015,0.015)))

row_ha = 
  rowAnnotation(bar2 = anno_barplot(d_hm_ss |> colMeans(),ylim = c(-0.005,0.005)))

p_hm_ss = 
  Heatmap(d_hm_ss |> t(),
        col = col_fun,
        top_annotation = column_ha, 
        right_annotation = row_ha,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title = "Marginal Effect"))


# Negative Symptoms and Personality Items -----

#Look at how these synthetic scores are distributed for CAPE and SAPAS
d_sapas <- d_train |> select(id,sapas_total)
d_cape  <- d_train |> select(id,cape_total)

p_pers_scales <- 
  bind_rows(d_sapas |>  pivot_longer(-id) ,
            d_cape  |> pivot_longer(-id)) |>
  mutate(name = str_remove(name,"_total")) |>
  filter(name %in% c("sapas","cape"))  |>
  mutate(name = toupper(name)) |>
  ggplot(aes(value)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~name,ncol = 2) +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(0,30,4))  +
  theme(panel.grid = element_blank(),
        axis.title  = element_text(size = 8),
        axis.text   = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_blank()) +
  labs(x = "Score", 
       y = "Count")


## Model -----

#We are going to model from the sum total scores
pers_trial = tibble(pers = c("sapas_total","cape_total"),
                    nt = c(8 ,     30))

d_pers = 
  left_join(d_sapas |> select(id,sapas_total),
            d_cape  |> select(id,cape_total),
            by = "id") |>
  right_join(d_im,
            by = "id") |>
  left_join(d_covars,by = "id") |>
  
  pivot_longer(all_of(im_names), 
               names_to = "infl", values_to = "infl_value") |> 
  pivot_longer(sapas_total:cape_total, 
               names_to = "pers", values_to = "pers_value") |> 
  drop_na() |>
  nest_by(infl,pers) |>
  left_join(pers_trial,by = "pers")


#I think we need to stick to a common modelling approach. The binomial has the
#advantage that the results are a probability so will be on the same scale no
#matter how many items are included. 

d_pers = 
  d_pers |>
  mutate(model = glm(cbind(pers_value, nt - pers_value) ~ infl_value * (sex + bmi24 + smk24 + audit_c),
                     family = binomial(link = 'logit'),
                     data = data,
                     method = "brglmFit") |> list(),
         me = avg_comparisons(model,variables = list(infl_value = "sd"),
                              by = TRUE,
                              type = "response") |>
           as_tibble() |>
           list()
  )

d_pers |>
  select(infl,pers,me) |>
  unnest(me) |>
  ungroup() |>
  mutate(p.adj = p.adjust(p.value,method = 'fdr')) |>
  filter(p.adj < 0.05) |>
  ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high,y = infl)) +
  geom_pointrange() +
  facet_wrap(~pers,ncol = 4) +
  geom_vline(xintercept = 0,lty = 2)


## Heatmap ------

d_hm_pers = 
  d_pers |>
  select(infl,pers,me) |>
  unnest(me) |>
  ungroup() |>
  select(infl,pers,estimate) |>
  pivot_wider(id_cols = "infl",names_from = "pers",values_from = "estimate") %>%
  arrange(infl) |>
  select(-infl) |>
  as.matrix()

rownames(d_hm_pers) <- 
  d_pers |>
  ungroup() |>
  mutate(infl = str_remove(infl,"_f24"))  |>
  distinct(infl) |>
  pull(infl)

column_ha = 
  HeatmapAnnotation(bar1 = anno_barplot(d_hm_pers |> rowMeans(),ylim = c(-0.015,0.015)))

row_ha = 
  rowAnnotation(bar2 = anno_barplot(d_hm_pers |> colMeans(),ylim = c(-0.005,0.005)))


p_hm_pers = 
  Heatmap(d_hm_pers |> t(),
          col = col_fun,
          top_annotation = column_ha, 
          right_annotation = row_ha,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          heatmap_legend_param = list(title = "Marginal Effect"))


# CISR Subscales -----

## Model =====

d_cisr_scales = 
  d_train |>
  select(c(id,som:pan)) |>
  mutate(across(som:pan,~as.integer(.x)-1))

#Quick plot
p_cisr_scales <- 
  d_cisr_scales |>  
  pivot_longer(-id) |>
  ggplot(aes(value)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~name,ncol = 6) +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(0,30,2))  +
  theme(panel.grid = element_blank(),
        axis.title  = element_text(size = 8),
        axis.text   = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_blank()) +
  labs(x = "Score", 
       y = "Count")


#Make the number of trials for binomial regression
cisr_trial = 
  tibble(ss = c("som","ftg","con","slp","irt","dep","did","wor","anx","pho","pan"),
         nt = 4) |>
  mutate(nt = case_when(ss == "did" ~ 5,
                        TRUE ~ nt))


d_cisr_s = 
  left_join(d_im,
            d_cisr_scales |>
              select(id,som,ftg,con,slp,irt,dep,did,wor,anx,pho,pan) ,
            by = "id") |>
  left_join(d_covars,by = "id") |>
  
  pivot_longer(all_of(im_names), 
               names_to = "infl", values_to = "infl_value") |> 
  pivot_longer(som:pan, 
               names_to = "ss", values_to = "s_value") |> 
  drop_na() |>
  nest_by(infl,ss) |>
  left_join(cisr_trial,by = "ss")


d_cisr_s = 
  d_cisr_s |>
  mutate(model = glm(cbind(s_value, nt - s_value) ~ infl_value * (sex + bmi24 + smk24 + audit_c),
                     family = binomial(link = 'logit'),
                     data = data,
                     method = "brglmFit") |> list(),
         me = avg_comparisons(model,variables = list(infl_value = "sd"),
                              by = TRUE,
                              type = "response") |>
           as_tibble() |>
           list()
  )




#Do a point-range plot of the model marginal effects
d_cisr_s |>
  select(infl,ss,me) |>
  unnest(me) |>
  ungroup() |>
  mutate(p.adj = p.adjust(p.value,method = 'fdr')) |>
  filter(p.value < 0.05) |>
  ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high,y = infl)) +
  geom_pointrange() +
  facet_wrap(~ss,ncol = 6) +
  geom_vline(xintercept = 0,lty = 2)


## Heatmap -----

#Make a heatmap with clustering of all inflammation x symptom effects

d_hm = 
  d_cisr_s |>
  select(infl,ss,me) |>
  unnest(me) |>
  ungroup() |>
  # mutate(estimate = case_when(p.value < 0.05 ~ estimate,
  #                             p.value >= 0.05 ~ 0)) |>
  select(infl,ss,estimate) |>
  pivot_wider(id_cols = "infl",names_from = "ss",values_from = "estimate") %>%
  arrange(infl) |>
  select(-infl) |>
  as.matrix()

rownames(d_hm) <- 
  d_cisr_s |>
  ungroup() |>
  mutate(infl = str_remove(infl,"_f24"))  |>
  distinct(infl) |>
  pull(infl)


column_ha2 = 
  HeatmapAnnotation(bar1 = anno_barplot(d_hm |> rowMeans(),ylim = c(-0.015,0.015)))

row_ha2 = 
  rowAnnotation(bar2 = anno_barplot(d_hm |> colMeans(),ylim = c(-0.005,0.005)))


p_hm = 
  d_hm |>
  t() |>
  Heatmap(clustering_method_rows = "complete",
          clustering_distance_rows = "euclidean",
          clustering_method_columns = "complete",
          clustering_distance_columns = "euclidean",
          col = col_fun,
          na_col = "grey",
          top_annotation = column_ha2, 
          right_annotation = row_ha2,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          heatmap_legend_param = list(title = "Marginal Effect"))


# Heatmap Assessbly -----


#What would be ideal is to glue together all the heatmaps. Fortunately the 
#complexheatmap package allows this
ht_list = p_hm %v% p_hm_ss %v% p_hm_pers %v%  p_hm_icd
draw(ht_list)

pdf(("./Figures/f24_heatmap_symptom.pdf"), width = 12, height = 5)
draw(ht_list)
dev.off()

#Plot symptoms ------

#For a supplement we assemble our scale histograms

p_symptoms <- 
  (p_cisr_scales / (p_cisr_scores|p_pers_scales)) + 
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(2,1))


ggsave(filename = "./Figures/f24_symptom_scale_histograms.pdf",
       height = 4,
       width  = 8,
       plot = p_symptoms)
  

# Cluster of inflammatory variables ------

#From our massive heat map, we can see a number of clusters within the 
#inflammation data - lets interrogate that further, and  rustle up the proper 
#variable names from ALSPAC

#I note that the clustering in Heatmap is done with hclust, so we ought to be
#able to produce it here (based on the package code)
d_clust_infl = 
  d_hm |>
  #cbind(d_hm,d_hm_cape) |>
  dist() |>
  hclust(method = "complete") 


plot(d_clust_infl )

#Note that we are cutting into 3 clusters based on manually inspecting the 
#clustering and based on discussion that we seem to get three distinct groups
#of associations

clust_infl_names = 
  d_clust_infl |> 
  stats::cutree(k = 3) |>
  as_tibble(rownames = "variable") |>
  rename(cluster_id = value)

#Have a look inside
clust_infl_names |>
  arrange(-cluster_id) |>
  # print(n = 100) |> 
  count(cluster_id)


#Plot the mean effect size for each cluster/symptom pair (mean over variables
#within a cluster)
p_cisr_by_cluster <- 
  d_hm |> 
  as_tibble(rownames = "infl") |> 
  #left_join(d_hm_icd |> as_tibble(rownames = "infl") |> rename(cisr_dep = estimate),by = "infl") |>
  left_join(clust_infl_names|> rename(infl = variable),by = "infl") |>
  group_by(cluster_id) |>
  summarise_at(vars(anx:wor), mean, na.rm = TRUE) |>
  pivot_longer(-cluster_id) |>
  
  mutate(name = factor(name, 
                       levels = c("slp","did","dep","ftg","som","pho","pan","irt","con","anx","wor"))) |>
  
  rename(estimate = value) |>
  ggplot(aes(x = cluster_id,y = name,fill = estimate)) +
  geom_tile() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x.bottom = element_text(angle = 0,size = 8,family = "sans",colour = "black"),
        axis.text = element_text(angle = 0,size = 8,family = "sans",colour = "black"),
        axis.title = element_text(angle = 0,size = 8,family = "sans",colour = "black"),
        legend.text = element_text(angle = 0,size = 8,family = "sans",colour = "black"),
        legend.title = element_text(angle = 0,size = 8,family = "sans",colour = "black")) +
  scico::scale_fill_scico(palette = "vik",direction = 1,limits = c(-0.03,0.03)) +
  labs(x = "Cluster", y = "CISR Symptom Scale") + 
  coord_fixed()

#So cluster 3 is the one with the positive associations with symptoms, 
#particularly fatigue and sleep; cluster 1 seems particularly associated with anxiety/
#worry

clust_3 = 
  clust_infl_names|>
  filter(cluster_id == 3)

#cluster 3 is our inflammation cells + friends cluster

clust_1 = 
  clust_infl_names|>
  filter(cluster_id == 1)

#Cluster 1 is our liver enzymes and crp


## Cluster 3 -----

#Make a table of the names of the variables in our high association cluster
tab_cluster_3 <- 
  v_0 |> 
  filter(str_detect(name,"_F24")) |> 
  select(name, lab) |> 
  as_tibble() |>
  mutate(name = str_remove(name,"_F24") |> tolower()) |>
  filter(name %in% clust_3$variable) |>
  arrange(name) |>
  rename(Variable = name,
         Description = lab) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)

tab_cluster_3

#Make a plot of these variables
p_cisr_symptom_cluster_3 =
  d_cisr_s |>
  mutate(infl = str_remove(infl,"_f24")) |>
  filter(infl %in% clust_3$variable) |>
  select(infl,ss,me) |>
  unnest(me) |>
  ungroup() |>
  mutate(p.adj = p.adjust(p.value,method = 'fdr')) |>
  filter(p.adj < 0.05) |>
  select(infl,ss,estimate) |>
  bind_rows(tibble(infl = "wbc",ss = "irt",estimate = NA_integer_)) |>
  mutate(ss = factor(ss, 
                     levels = c("slp","did","dep","ftg","som","pho","pan","irt","con","anx","wor"))) |>
  arrange(ss) |>
  ggplot(aes(x = infl, y = ss, fill = estimate)) +
  geom_tile() +
  theme(panel.background = element_rect(fill = "grey95"),
        panel.grid = element_blank(),
        axis.text.x.bottom = element_text(angle = 90,size = 8,family = "sans",colour = "black"),
        axis.text = element_text(angle = 0,size = 8,family = "sans",colour = "black"),
        axis.title = element_text(angle = 0,size = 8,family = "sans",colour = "black"),
        legend.text = element_text(angle = 0,size = 8,family = "sans",colour = "black"),
        legend.title = element_text(angle = 0,size = 8,family = "sans",colour = "black")) +
  scale_fill_scico(palette = "vik",direction = 1,limits = c(-0.03,0.03),oob = scales::squish) +
  labs(x = "Immunometabolic Variable", y = "CISR Symptom Scale") + 
  coord_fixed()

## Cluster 1 ------

#Make a table of the names of the variables in our negative association cluster
tab_cluster_1 <- 
  v_0 |> 
  filter(str_detect(name,"_F24")) |> 
  select(name, lab) |> 
  as_tibble() |>
  mutate(name = str_remove(name,"_F24") |> tolower()) |>
  filter(name %in% clust_1$variable) |>
  arrange(name) |>
  rename(Variable = name,
         Description = lab) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)

tab_cluster_1

#Make a plot of these variables
p_cisr_symptom_cluster_1 =
  d_cisr_s |>
  mutate(infl = str_remove(infl,'_f24')) |>
  filter(infl %in% clust_1$variable) |>
  select(infl,ss,me) |>
  unnest(me) |>
  ungroup() |>
  mutate(p.adj = p.adjust(p.value,method = 'fdr')) |>
  filter(p.adj < 0.05) |>
  mutate(ss = factor(ss, 
                     levels = c("slp","did","dep","ftg","som","pho","pan","irt","con","anx","wor"))) |>
  arrange(ss) |>
  ggplot(aes(x = infl, y = ss, fill = estimate)) +
  geom_tile() +
  theme(panel.background = element_rect(fill = "grey95"),
        panel.grid = element_blank(),
        axis.text.x.bottom = element_text(angle = 90,size = 8,family = "sans",colour = "black"),
        axis.text = element_text(angle = 0,size = 8,family = "sans",colour = "black"),
        axis.title = element_text(angle = 0,size = 8,family = "sans",colour = "black"),
        legend.text = element_text(angle = 0,size = 8,family = "sans",colour = "black"),
        legend.title = element_text(angle = 0,size = 8,family = "sans",colour = "black")) +
  scale_fill_scico(palette = "vik",direction = 1,limits = c(-0.03,0.03)) +
  labs(x = "Immunometabolic Variable", y = "CISR Symptom Scale") +
  coord_fixed()

## Plot Clustering ------

p_cisr_symptom_cluster <- 
  p_cisr_by_cluster+ p_cisr_symptom_cluster_1 + p_cisr_symptom_cluster_3  + 
  plot_layout(guides = "collect",widths = c(3,5,22)) +
  plot_annotation(tag_levels = "A")


ggsave(filename = "./Figures/f24_cisr_cluster.pdf",
       height = 4,
       width  = 8,
       plot = p_cisr_symptom_cluster)


## Save clusters -----

# Lets also save our cluster identities for future use
v_0 |> 
  filter(str_detect(name,"_F24")) |> 
  select(name, lab) |> 
  as_tibble() |>
  mutate(name = str_remove(name,"_F24") |> tolower()) |>
  right_join(clust_infl_names |> rename(name = variable), by = 'name') |>
  arrange(name) |>
  rename(infl = name,
         desc = lab,
         clus = cluster_id) |>
  
  write_csv('./Models/f24_infl_variable_clusters.csv')