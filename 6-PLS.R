
#Introduction ------

#We are going to do PLS, only with the variables that we identified as being 
#associated with depressive symptoms

pacman::p_load(tidyverse,
               tidymodels,
               marginaleffects,
               patchwork,
               ggdist,
               rsample,
               mixOmics)

conflicted::conflict_prefer('select', 'dplyr')
conflicted::conflict_prefer("rcc"   , "mixOmics")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("map"   , "purrr")

#Standard z score
zscore = function(x){
  z = (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  
  return(z)
}


source('./Final/0-Common_data.R')

#Load up our pre-split data so we can focus on the training data only
d_split <- 
  read_rds(paste(data_dir,"alspac_data_split.rds",sep = '//'))

d_train <- 
  d_split |>
  training()


#We are going to do PLS between a subset of the depression items from CISR and
#the full set of Olink proteins

## Prepare the immuno-metabolic dataset -----

#Prepare the F24 immunometabolic dataset
d_im = 
  d_train |>
  select(id,cisr_dep,contains("_f24")) |>
  rename_with(.cols = -c(id, cisr_dep),.fn = ~str_remove(.x,"_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(-c(id,cisr_dep),zscore))


##Prepare our covariates -------

#We might need some covariates?
#These are the covariates that our DAG suggests form the minimum adjustment set
covars = c("sex","bmi24","smk24","audit_c")

d_covars = 
  d_cv |>
  select(id, all_of(covars)) 


## Prepare variable names =========

#We need to make a set of names for our sets of variables

#Get the names of the immunometabolic variables
vars_metab = d_im |> select(-c(id,cisr_dep)) |> colnames()

infl_clus <- 
  read_csv('./Models/f24_infl_variable_clusters.csv')


select_bloods <- 
  infl_clus |> filter(clus != 2) |> pull(infl)


# Quick correlation plot -------

#Make a correlation plot of the F24 bloods in clusters 1 and 3 (get rid of the
#non-selected variables)
cor_f24 = 
  d_train |>
  select(contains("_f24")) |>
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  dplyr::select(all_of(infl_clus |> filter(clus != 2) |> pull(infl))) |>
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

#Save
ggsave(filename = "./Figures/f24_selected_Vars_correlation.pdf",
       height = 8,
       width  = 10,
       plot = p_f24_cor)

# PLS2 -----

#With the implementation of PLS in mixOmics we can do everything all at once, 
#which they call PLS2. There is also a paper that points out how CCA and PLS are 
#very similar apart from the implementation of different regularization types
#(L1 vs L2)

## Build the symptom scales ======

#All subscales are scored 0-4, except depressive ideas, which ranges 0-5
d_cisr_scales = 
  d_train |> 
  select(c(id,som:pan)) |> 
  mutate(across(som:pan,~as.integer(.x)-1))
 
#Make plots of the scores
d_cisr_scales |>
  select(id,som,ftg,con,slp,irt,dep,did,wor,anx,pho,pan) |>
  pivot_longer(-id) |>
  
  ggplot(aes(value)) +
  geom_histogram(binwidth = 1) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~name,ncol = 4)


## Assemble Datasets -----

#Do some setup to avoid problems with clashing function names
vars_cisr = c("som","ftg","con","slp","irt","dep","did","wor","anx","pho","pan")

d_complete =
  left_join(d_im,
            d_cisr_scales,
            by = "id")|>
  
  #Complete cases only
  drop_na() |>
  select(id,c(all_of(select_bloods),all_of(vars_cisr)))

x = 
  d_complete |>
  select(all_of(select_bloods)) |>
  as.matrix()

y = 
  d_complete |>
  select(all_of(vars_cisr)) |>
  as.matrix()

#So we have 2551 individuals contributing complete case data after imputing as above
#dim(x)
#dim(y)

## Explore dataset =====

# produce a heat map of the cross correlation matrix
imgCor(x, y, sideColors = c("purple", "green"), color = color.jet(21, alpha = 1)) 

#Following in the footsteps of the mixomics case study, we will
#start with PCA.
pca.cisr  <- pca(y,ncomp = 10,center = FALSE,scale = FALSE)
pca.metab <- pca(x,ncomp = 10,center = FALSE,scale = FALSE)

plot(pca.cisr)  # arguably 1 component 
plot(pca.metab) # arguably 1 or 3 components


#Look at the clustering of individuals in PC subspace
plotIndiv(pca.cisr,comp = c(1,2),
          title = 'CISR, PCA comp 1 - 2')

#There is that slightly odd stripey clustering we see elsewhere due to the
#cisr variables not being continuous

plotIndiv(pca.metab,comp = c(1,2),
          title = 'Metabolic, PCA comp 1 - 2')

#PCs look continuously distributed and non-correlated, which we would expect


## Initial, untuned PLS model =====

#We need to think about what mode to run pls in. Are we trying to use
#the metabolic data to predict the cisr data, or to use the cisr data
#to predict the metabolic data? In biology both are possibilities e.g.
#increases to inflammation increase MH problems -  the converse is probably
#less likely in cross sectional data - you could imagine something different
#over time - maybe if you are more depressed you smoke more and that causes more
#inflammation in the longer term or something. 

#We are going to use "regression" mode because we are most interested 
#in predicting depressive symptoms (i.e. the Ys) from immunometabolic
#variables (i.e. the Xs). If we were interested in both being equivalent
#we could use "canonical" mode 

spls.cisr <- spls(X = x, Y = y, ncomp = 5,mode = 'regression',scale = FALSE)

#So this is our basic model but it needs tuning and sparsifying

## Tune model ======

set.seed(06042023)

#First lets select the number of components to include

perf.spls.cisr <- perf(spls.cisr,validation = 'Mfold',
                        folds = 10, nrepeat = 10)

#Various methods exist to compare performance at different numbers
#of components

#The mixomics case study suggest we use the Q total measure
plot(perf.spls.cisr,criterion = 'Q2.total')
plot(perf.spls.cisr,criterion = 'Q2')

#We are below the line for all of Q2.total. For the individual Y measures ("Q2")
#we see that above one component we get a drop below the line, so we take two
#forward as in the package case study

#Now we ask how many variables to select

#There are two measures to use - the correlation between predicted and actual
#omponents, and the residual ss  (RSS). RSS gives more weights to large errors 
#and is sensitive to outliers but it does select less features on Y than the cor
#measure.

#We are more interested in robustness to outliers and less interested in
#deflation so lets go with cor as they do in the case study

#Set a range of test values for number of variables to keep for the x dataframe

#The choice of grid will to some extent determine the outcome - is there a 
#minimum number of variables we want to keep from each  set? Lets say we might 
#want >1 from the cisr data to grant us  the ability to look further into the
#nature of there variables selected. We want the same - >1  for the metabolic 
#variables to represent different groupings within them

#We are going to set our seed at each iteration of this process in order to
#have reproducibility
set.seed(06042023)

list.keepX <- c(seq(2,38,1))
list.keepY <- c(seq(2,11,1))

#Set up for parallel processing (still very slow I'm afraid)
BPPARAM <- BiocParallel::SnowParam(workers = 15)

#We are going to start with lots of components, and then see what is optimal
tune.spls.cisr <- tune.spls(X = x, Y = y, ncomp = 2,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 2, folds = 10,
                             mode = 'regression',
                             measure = 'cor',
                             BPPARAM = BPPARAM)


plot(tune.spls.cisr)

#The two measures are not that dissimilar, which is encouraging
tune.spls.cisr$choice.keepY
tune.spls.cisr$choice.keepX


#Now we store the optimal number of x and y components
optimal.keepX <- tune.spls.cisr$choice.keepX
optimal.keepY <- tune.spls.cisr$choice.keepY
optimal.ncomp <- length(optimal.keepX) # extract optimal number of components


## Fit Final Model ======

#Now we can fit the final model

final.spls.cisr <- spls(X = x,
                        Y = y, 
                         ncomp = optimal.ncomp,
                         keepX = optimal.keepX,
                         keepY = optimal.keepY,
                         mode = 'regression')

#Lets save this model given how long it takes to get all the optimised pls
#component counts. Given, however, that the model structure contains the ALSPAC
#data itself, rather than just coefficients etc, we save is securely - if you
#reuse the code you would need to do the same


## Save final model =======

write_rds(final.spls.cisr,paste(data_dir,"cisr_pls2_model.rds",sep = '/'))
#final.spls.cisr = read_rds(paste(data_dir,"cisr_pls2_model.rds",sep = '/'))


## Inspect final model =======

#Now we can interrogate this optimised model

#Explained variance - not so bad
final.spls.cisr$prop_expl_var |>
  as_tibble() |>
  mutate(component = 1:n()) |>
  transmute(component,
            cumulative_metab = cumsum(X),
            cumulative_cisr  = cumsum(Y)) |>
  mutate(across(where(is.double),~round(.x*100,digits = 1))) 

### Loading Plots -----

#Plot loadings our way

#Note that we are going to flip the loadings and the exported values for 
#components 1 and 3 because we want to have the main direction of loadings be
#positive rather than negative. The valence of the PLS components is arbitrary
#so there is no specific problem with doing this, its purely for interpretability
p_loadings_2 = 
  bind_rows(
    final.spls.cisr$loadings$X |>
    as_tibble(rownames = "variable") |>
    mutate(set = "infl"),
    final.spls.cisr$loadings$Y |>
      as_tibble(rownames = "variable") |>
      mutate(set = "cisr"))    |>

  #Do the flip
  mutate(
    across(.cols = starts_with("comp1"),.fns = ~-.x)) |>

  pivot_longer(-c(set,variable)) |>
  filter(value != 0) |>
  nest_by(set,name) |>
  mutate(plot_col = case_when(name == "comp1" ~ "#1b9e77",
                              name == "comp2" ~ "#d95f02")) |>
  mutate(plot_title = case_when(name == "comp1" ~ "Component 1",
                                name == "comp2" ~ "Component 2")) |>
  mutate(plots = list(
    data |>
      mutate(variable = fct_reorder(variable,value,max)) |> 
      ggplot(aes(x = variable, y = value)) + 
      geom_col(fill = plot_col) + 
      geom_hline(yintercept = 0,lty = 1) +
      theme(#panel.grid = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x.bottom = element_text(angle = -90,vjust = 0.5,size = 6),
        #axis.title.x.bottom = element_text(size = 6),
        axis.title.x.bottom = element_blank(),
        axis.text.y = element_text(size = 6),
        title = element_text(size = 6),
        axis.title = element_blank()) +
      coord_cartesian(ylim = c(-0.7,0.6))+
      labs(title = plot_title))) %>%
  pluck('plots') |>
  patchwork::wrap_plots(ncol = 2) 

p_loadings_2

#Loading Relationships

# - Comp 1 is fatigue/sleep and our acute inflammatory variables
# - Comp 2 is the anxiety and worry vs alt and ada 


## Use Scores -----

### Covariates -----

#Explore PLS score - covariate relationships

#First get the scores
d_pls2 <- 
  bind_cols(final.spls.cisr$variates$X |>
              as_tibble() |>
              rename_with(~paste("infl",.,sep = "_")),
            final.spls.cisr$variates$Y |>
              as_tibble() |>
              rename_with(~paste("cisr",.,sep = "_"))) |>
  bind_cols(d_complete |> select(id) |> left_join(d_cv,by = "id")) |>
  
  #Do the flip
  mutate(infl_comp1 = -infl_comp1,
         cisr_comp1 = -cisr_comp1) |>

  left_join(d_24 |> select(id,cisr_dep), by = "id") |>
  left_join(d_train |> select(id, trauma_total,ace_total_classic), by = "id")|>
  select(-c(eth,bmi9,bp24,cisr_dep,starts_with("smfq")))


#Lets try and regress the PLS components with covariates
p_pls2_covars <- 
  d_pls2  |>
  pivot_longer(infl_comp1:cisr_comp2,names_to = "comp",values_to = "pls_score") |>
  nest_by(comp) |>
  mutate(m1 = 
           list(
             lm(pls_score ~ sex + bmi24 + smk24 + audit_c + ace_total_classic + trauma_total + md + me + sc,
                data = data)
           ),
         mt = list(
           marginaleffects::avg_comparisons(m1) |> as_tibble()
           )) |>
  select(comp,mt) |>
  unnest(cols = mt) |>
  ungroup() |>
  filter(str_detect(comp,"infl")) |>
  mutate(comp = case_when(comp == "infl_comp1" ~ "Component 1",
                          comp == "infl_comp2" ~ "Component 2")) |>
  mutate(term = case_when(term == "ace_total_classic" ~ "Total ACEs",
                          term == "audit_c" ~ "AUDIT-C@24",
                          term == "bmi24" ~ "BMI@24",
                          term == "md" ~ "Maternal EPDS",
                          term == "me" ~ "Maternal Education",
                          term == "sc" ~ "Maternal Social Class",
                          term == "sex" ~ "Sex",
                          term == "smk24" ~ "Smoking@24",
                          term == "trauma_total" ~ "Total Traumatic Events")) |>
  mutate(p.adj = p.adjust(p.value,method = "fdr")) |>
  mutate(p_star = case_when(p.adj  < 0.05 ~ '*',
                            p.adj >= 0.05 ~ ' ')) |>
  ggplot(aes(x = term,y = comp,fill = statistic,label = p_star)) +
  geom_tile() +
  geom_text() +
  scico::scale_fill_scico(palette = "vik",limits = c(-25,25))+
  theme(axis.text.x.bottom = element_text(angle = -90,size = 6),
        axis.text.y.left = element_text(size = 6),
        axis.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank()) +
  coord_equal() +
  labs(x = NULL, y = NULL)


### SMFQ ------

#Future SMFQs

d_smfq <- 
  d_pls2 |>
  select(id,contains("_comp"),all_of(covars)) |>
  left_join(d_train |>
              select(id,smfq_25,smfq_c1,smfq_c2,smfq_c4,smfq_c5),
            by = "id") |>
  pivot_longer(infl_comp1:cisr_comp2,names_to = "comp",values_to = "pls_score") |>
  pivot_longer(starts_with("smfq_"),values_to = "smfq_score",names_to = "smfq_wave") |>
  nest_by(smfq_wave,comp) |>
  mutate(m1 = 
           list(
             glm(smfq_score ~ pls_score * (sex + bmi24 + smk24 + audit_c),
                 family = poisson(link = "log"),
                data = data)
           ),
         mt = list(
           avg_comparisons(model = m1,
                           variables = list(pls_score = "sd"),
                           type = "response",
                           conf_level = 0.95) |>
             as_tibble()
         )) 

#Or we could do lines for this - 
p_pls2_smfq <- 
  d_smfq|>
  select(smfq_wave,comp,mt) |>
  unnest(cols = mt) |>
  ungroup() |>
  filter(str_detect(comp,"infl_")) |>
  mutate(comp = case_when(comp == "infl_comp1" ~ "Component 1",
                          comp == "infl_comp2" ~ "Component 2",
                          comp == "infl_comp3" ~ "Component 3")) |>
  mutate(p.adj = p.adjust(p.value,method = "fdr")) |>
  # filter(p.adj < 0.05) |>
  mutate(p_star = case_when(p.adj  < 0.05 ~ '*',
                            p.adj >= 0.05 ~ ' ')) |>
  ggplot(aes(x = smfq_wave,y = estimate,ymin = conf.low,ymax = conf.high,
             colour = comp,
             label = p_star,
             group = comp)) +
  geom_line() +
  geom_pointrange()+
  geom_text(colour = "black",size = 6,nudge_y = 0.25) +
  geom_hline(yintercept = 0,lty = 2) +
  facet_wrap(~comp,ncol = 3) +
  scale_colour_brewer(palette = "Dark2")+
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        strip.text = element_text(size = 6),
        legend.position = "none",
        panel.grid = element_blank())+
  labs(x = "SMFQ Wave",y = "Marginal Change in SMFQ \nfor 1 unit PLS Score")

#So our scores are associated with our future SMFQ scores, particularly comp1.


## Make Figure ======

#We need to get some colours for the PLS components and use them consistently
#across the figure. Also need to make them different to the cluster colours
#that we used in an earlier figure

#So we also need to get the loading plots...

p_pls2_mega <- 
  (p_loadings_2/p_pls2_covars/p_pls2_smfq) +
  plot_annotation(tag_levels = 'A') + 
  plot_layout(heights = c(2,1,1))

ggsave(filename = "./Figures/f24_pls2_mega.pdf",
       height = 7,
       width  = 6,
       plot = p_pls2_mega)


# Test PLS on the test data -----

#So we made our PLS model on the training data, what if we extract scores from
#the test data and see if the same relationships to covariates etc are present?

#Extract same components from the test data

d_test <- 
  d_split |>
  testing()

y_test <- 
  d_test |>
  select(id,cisr_dep,contains("_f24")) |>
  rename_with(.cols = -c(id, cisr_dep),.fn = ~str_remove(.x,"_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(-c(id,cisr_dep),zscore)) |>
  select(all_of(select_bloods)) |>
  as.matrix()


pls_test <- 
  predict(final.spls.cisr,
          newdata = y_test)

d_pls_test <- 
  pls_test$variates |>
  as_tibble() |>
  rename_with(.cols = everything(),.fn = ~str_replace(.x,"dim","infl_comp")) |>
  bind_cols(d_test |> select(id) |> left_join(d_cv,by = "id")) |>
  
  #Do the flip
  mutate(infl_comp1 = -infl_comp1) |>
  
  left_join(d_24 |> select(id,cisr_dep), by = "id") |>
  left_join(d_test |> select(id, trauma_total,ace_total_classic), by = "id")|>
  select(-c(eth,bmi9,bp24,cisr_dep,starts_with("smfq")))


p_pls_test_12 <- 
  d_pls_test |>
  ggplot(aes(x = infl_comp1,y = infl_comp2)) +
  geom_point(alpha = 0.2) +
  geom_rug() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey95"),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6))  +
  labs(x = "PLS2 Component 1", y = "PLS2 Component 2")+
  coord_fixed(xlim = c(-7.5,7.5),ylim = c(-7.5,7.5))


#Well those look OK. What about some models

p_pls_test_covars_pair <- 
  d_pls_test  |>
  pivot_longer(infl_comp1:infl_comp2,names_to = "comp",values_to = "pls_score") |>
  nest_by(comp) |>
  mutate(m1 = 
           list(
             lm(pls_score ~ sex + bmi24 + smk24 + audit_c + ace_total_classic + trauma_total + md + me + sc,
                data = data)
           ),
         mt = list(
           marginaleffects::avg_comparisons(m1) |> as_tibble()
         )) |>
  select(comp,mt) |>
  unnest(cols = mt) |>
  ungroup() |>
  filter(str_detect(comp,"infl")) |>
  mutate(comp = case_when(comp == "infl_comp1" ~ "Component 1",
                          comp == "infl_comp2" ~ "Component 2")) |>
  mutate(term = case_when(term == "ace_total_classic" ~ "Total ACEs",
                          term == "audit_c" ~ "AUDIT-C@24",
                          term == "bmi24" ~ "BMI@24",
                          term == "md" ~ "Maternal EPDS",
                          term == "me" ~ "Maternal Education",
                          term == "sc" ~ "Maternal Social Class",
                          term == "sex" ~ "Sex",
                          term == "smk24" ~ "Smoking@24",
                          term == "trauma_total" ~ "Total Traumatic Events")) |>
  mutate(p.adj = p.adjust(p.value,method = "fdr")) |>
  mutate(p_star = case_when(p.adj < 0.05 ~ '*',
                            p.adj >= 0.05 ~ ' ')) |>
  ggplot(aes(x = term,y = comp,fill = statistic,label = p_star)) +
  geom_tile() +
  geom_text() +
  scico::scale_fill_scico(palette = "vik",limits = c(-25,25))+
  theme(axis.text.x.bottom = element_text(angle = -90,size = 6),
        axis.text.y.left = element_text(size = 6),
        axis.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank()) +
  coord_equal() +
  labs(x = NULL, y = NULL)

#So the relationships are attenuated but mostly remain; the direction of the 
#associations remains too

#The predictive SMFQ is probably the most important
d_smfq_test <- 
  d_pls_test |>
  select(id,contains("_comp"),all_of(covars)) |>
  left_join(d_test |>
              select(id,smfq_25,smfq_c1,smfq_c2,smfq_c4,smfq_c5),
            by = "id") |>
  pivot_longer(infl_comp1:infl_comp2,names_to = "comp",values_to = "pls_score") |>
  pivot_longer(starts_with("smfq_"),values_to = "smfq_score",names_to = "smfq_wave") |>
  nest_by(smfq_wave,comp) |>
  mutate(m1 = 
           list(
             glm(smfq_score ~ pls_score *(sex + bmi24 + smk24 + audit_c),
                 family = poisson(link = "log"),
                 data = data)
           ),
         mt = list(
           avg_comparisons(model = m1,
                           variables = list(pls_score = "sd"),
                           type = "response",
                           conf_level = 0.95) |>
             as_tibble()
         )) 

#Or we could do lines for this - 
p_pls_test_smfq <- 
  d_smfq_test|>
  select(smfq_wave,comp,mt) |>
  unnest(cols = mt) |>
  ungroup() |>
  mutate(comp = case_when(comp == "infl_comp1" ~ "Component 1",
                          comp == "infl_comp2" ~ "Component 2")) |>
  mutate(p.adj = p.adjust(p.value,method = "fdr")) |>
  mutate(p_star = case_when(p.adj < 0.05 ~ '*',
                            p.adj >= 0.05 ~ ' ')) |>
  ggplot(aes(x = smfq_wave,y = estimate,ymin = conf.low,ymax = conf.high,
             colour = comp,
             label = p_star,
             group = comp)) +
  geom_line() +
  geom_pointrange()+
  geom_text(colour = "black",size = 6,nudge_y = 0.25) +
  geom_hline(yintercept = 0,lty = 2) +
    facet_wrap(~comp,ncol = 1) +
  scale_colour_brewer(palette = "Dark2")+
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        strip.text = element_text(size = 6),
        legend.position = "none",
        panel.grid = element_blank())+
  labs(x = "SMFQ Wave",y = "Marginal Change in SMFQ for 1 unit PLS Score")

#It still works!!!
p_pls_test <- 
  (p_pls_test_covars_pair / p_pls_test_smfq) + 
  plot_layout(heights = c(1,2)) +
  plot_annotation(tag_levels = "A")

ggsave(filename = "./Figures/f24_pls_test_smfq.pdf",
       height = 6,
       width  = 5,
       plot = p_pls_test)


#Or tablise
d_smfq_test|>
  select(smfq_wave,comp,mt) |>
  unnest(cols = mt) |>
  ungroup() |>
  filter(str_detect(comp,"infl_")) |>
  mutate(comp = str_remove(comp,"infl_")) |>
  mutate(p.adj = p.adjust(p.value,method = "fdr")) |>
  filter(p.adj < 0.05) 

#So component 1 and 2 do consistently appear to be associated with the future
#SMFQ score