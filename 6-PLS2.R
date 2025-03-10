
#Introduction ------

#We are going to do PLS, only with the variables that we identified as being associated with depressive symptoms

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

#Load up our pre-prepared data
d <- 
  read_rds(paste(data_dir,"alspac_data_final.rds",sep = '//'))


# Assemble Data ------

## Covariates -------

#We might need some covariates?
#These are the covariates that our DAG suggests form the minimum adjustment set
covars = c("sex","bmi24","smk24","audit_c","ph")

d_covars = 
  d_cv |>
  select(id, all_of(covars)) 

## Immuno-metabolic dataset -----

#Prepare the F24 immunometabolic dataset
d_im = 
  d |>
  select(id,cisr_dep,contains("_f24")) |>
  rename_with(.cols = -c(id, cisr_dep),.fn = ~str_remove(.x,"_f24")) |>
  
  #Lets z score all the inflammatory markers so we get standardised coefficients
  mutate(across(-c(id,cisr_dep),zscore))


#We need to make a set of names for our sets of variables

#Get the names of the immunometabolic variables
vars_metab = d_im |> select(-c(id,cisr_dep)) |> colnames()

infl_clus <-
  read_csv('./Models/f24_infl_variable_clusters.csv')


select_bloods <-
  infl_clus |> filter(clus != 2) |> pull(infl)

## Symptom Scales ======

#Which scales do we use?
scales_use = c("ach","ftg" ,"con","slp_dec","slp_inc","app_dec","app_inc","irt",
               "sad","morn","dsx","rstl","slow"   ,"anh"    ,"sneg"   ,"mot"    ,"dcog",
               "stb","wor","anx","pho","pan")


d_cisr_scales = 
  d |> 
  select(c(id,all_of(scales_use))) |> 
  mutate(across(all_of(scales_use),~as.integer(.x)-1))


#Make plots of the scores
d_cisr_scales |>
  select(id,all_of(scales_use)) |>
  pivot_longer(-id) |>
  
  ggplot(aes(value)) +
  geom_histogram(binwidth = 1) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~name,ncol = 4)



## Assemble Datasets -----

#Do some setup to avoid problems with clashing function names
d_complete =
  left_join(d_im,
            d_cisr_scales,
            by = "id")|>
  left_join(d_covars, by = "id") |>
  
  #Complete cases only
  drop_na() |>
  select(id,c(all_of(select_bloods),all_of(scales_use),all_of(covars)))

#Note based on a biological psychiatry reviewer comment, we are going to include the covariates in the x matrix

x = 
  d_complete |>
  select(all_of(select_bloods),all_of(covars)) |>
  mutate(sex = if_else(sex == "Male",0,1),
         bmi24 = zscore(bmi24),
         audit_c = zscore(audit_c)) |>
  as.matrix()

y = 
  d_complete |>
  select(all_of(scales_use)) |>
  as.matrix()

#So we have 2311 individuals contributing complete case data 
#dim(x)
#dim(y)


# Correlations ----

## Infl correlation plot -------

#Make a correlation plot of the F24 bloods in clusters 1 and 3 (get rid of the
#low association cluster 2 variables)
cor_f24 = 
  d_im |>
  select(-c(id,cisr_dep)) |>
  rename_with(.cols = everything(),.fn = ~str_remove(.x,"_f24")) |>
  dplyr::select(all_of(select_bloods)) |>
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

#I think we save it at this point, although its not going into a figure
#in the final manuscript
ggsave(filename = "./Figures/f24_selected_Vars_correlation.pdf",
       height = 8,
       width  = 10,
       plot = p_f24_cor)

#There are a lot of correlations, mostly positive

#Look at the correlations 
cor_f24 |>
  as_tibble(rownames  = "var1" ) |>
  pivot_longer(-var1,names_to = "var2",values_to = "cor") |>
  mutate(sc = if_else(var1 == var2,T,F)) |>
  filter(!sc)
select(-sc) |>
  
  summarise(mu = mean(cor),
            sig = sd(cor),
            min_c = min(cor),
            max_c = max(cor))

## Symptom correlation plot ------

cor_sym = 
  d_cisr_scales |>
  select(-c(id)) |>
  dplyr::select(all_of(scales_use)) |>
  cor(use = "pairwise.complete.obs",
      method = "spearman") 

#Do our own clustering so we can save the variable names for later - I have 
#taken this code directly from the ggcorrplot github so it should produce 
#identical results to using the hc.ord method
clust_cor_sym <- 
  hclust(as.dist(1 - cor_sym),method = "complete")

#Extract the ordering based on the clustering
ord_sym <- clust_cor_sym$order

#Save variable names for later
vars_sym_ordered <- clust_cor_sym$labels[clust_cor_sym$order]

#Sort the variables
cor_sym <- cor_sym[ord_sym, ord_sym]

#Now make the plot
p_cor_sym <- 
  cor_sym %>% 
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

p_cor_sym

#I think we save it at this point, although its not going into a figure
#in the final manuscript
ggsave(filename = "./Figures/symptom_correlation.pdf",
       height = 8,
       width  = 10,
       plot = p_cor_sym)




#Look at the correlations 
cor_sym |>
  as_tibble(rownames  = "var1" ) |>
  pivot_longer(-var1,names_to = "var2",values_to = "cor") |>
  mutate(sc = if_else(var1 == var2,T,F)) |>
  filter(!sc) |>
  select(-sc) |>
  
  summarise(mu = mean(cor),
            sig = sd(cor),
            min_c = min(cor),
            max_c = max(cor))

# PCAs =====

# produce a heat map of the cross correlation matrix
imgCor(x, y, sideColors = c("purple", "green"), color = color.jet(21, alpha = 1)) 

#Correlations across datasets are lower than within

#Following in the footsteps of the mixomics case study, we will
#start with PCA.

## Immuno-metabolic variables -----

nc_x = 
  d_complete |>
  select(all_of(select_bloods)) |>
  parameters::n_components(type = "PCA")

pca.x <- 
  d_complete |>
  select(all_of(select_bloods)) |>
  as.matrix() |>
  pca(X = _,ncomp = 10,center = FALSE,scale = FALSE)

#Explained variance
pca.x$prop_expl_var$X |> as_tibble() |> mutate(cs = cumsum(value))

#Plot
p_x <- 
  pca.x$prop_expl_var$X |> 
  as_tibble() |>
  mutate(comp = 1:n(),
         value = value * 100) |>
  ggplot(aes(x = comp, y = value)) +
  geom_col() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_x_continuous(breaks = seq(2,10,2)) +
  labs(x = "Principal Component", y = "Explained Variance (%)",title = "Immuno-metabolic PCA")

#Plot scores
plotIndiv(pca.x,comp = c(1,2),
          title = 'Immuno-Metabolic, PCA comp 1 - 2')

#One big clump really

## Symptom variables ------

nc_y = 
  d_complete |>
  select(all_of(scales_use))|>
  parameters::n_components(type = "PCA")


#Run the PCA
pca.y  <- pca(y,ncomp = 10,center = FALSE,scale = FALSE)

p_y <- 
  pca.y$prop_expl_var$X |> 
  as_tibble() |>
  mutate(comp = 1:n(),
         value = value * 100) |>
  ggplot(aes(x = comp, y = value)) +
  geom_col() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_x_continuous(breaks = seq(2,10,2)) +
  labs(x = "Principal Component", y = "Explained Variance (%)",title = "Symptom PCA")


pca.y$prop_expl_var$X |> as_tibble() |> mutate(cs = cumsum(value))

#Look at the clustering of individuals in PC subspace

plotIndiv(pca.y,comp = c(1,2),
          title = 'Symptom, PCA comp 1 - 2')

#There is that slightly odd stripey clustering we see elsewhere due to the
#variables not being continuous

## Make a plot -----

#Assemble 
p_pca = 
  ((p_f24_cor|p_cor_sym) + plot_layout(guides = "collect"))/
  (p_x|plot_spacer()|p_y|plot_spacer()) +
  plot_layout(heights = c(4,1)) +
  plot_annotation(tag_levels = "A")

#And save
ggsave(filename = "./Figures/pca_plot.pdf",
       height = 7,
       width  = 8,
       plot = p_pca)



# Basic, untuned PLS model =====

#We need to think about what mode to run pls in. Are we trying to use
#the metabolic data to predict the cisr data, or to use the cisr data
#to predict the metabolic data? In biology both are possibilities e.g.
#increases to inflammation increase MH problems -  the converse is probably
#less likely in cross sectional data - you could imagine something different
#over time - maybe if you are more depressed you smoke more and that causes more
#inflammation in the longer term or something. 

#We are going to use "regression" mode because we are mostly interested 
#in predicting depressive symptoms (e.g. the Ys) from immunometabolic
#variables (e.g. the Xs). If we were interested in both being equivalent
#we might use "canonical" mode 

spls.cisr <- spls(X = x, Y = y, ncomp = 5,mode = 'regression',scale = FALSE)


#So this is our basic model but it needs tuning and sparsifying

# Tune model ======

set.seed(06042023)

#First lets select the number of components to include

perf.spls.cisr <- perf(spls.cisr,validation = 'Mfold',
                       folds = 10, nrepeat = 10)

#Various methods exist to compare performance at different numbers
#of components

#The mixomics case study suggest we use the Q total measure
plot(perf.spls.cisr,criterion = 'Q2.total')
plot(perf.spls.cisr,criterion = 'Q2')

#We are below the line after 1 for all of Q2.total. For the individual Y measures ("Q2")
#we see that above one component we get a drop below the line, this suggests
#that one or two components are sufficient

#Now we ask how many variables to select

#There are two measures to use - the correlation between
#predicted and actual components, and the residual ss
#RSS gives more weights to large errors and is sensitive to outliers
#but it does select less features on Y than the cor measure.

#We are probably more interested in robustness to outliers and less interests in
#deflation so lets go with cor as they do in the case study

#Set a range of test values for number of variables to keep
#for the x dataframe

#The choice of grid will to some extent determine the outcome
#is there a minimum number of variables we want to keep from each 
#set? Lets say we might want >1 from the cisr data to grant us 
#the ability to look further into the nature of there variables selected
#And we want at least 4 of the metabolic variables to represent different
#groupings within them

#We are going to set our seed at each iteration of this process in order to
#have reproducibility
set.seed(06042023)

list.keepX <- c(seq(2,50,1))
list.keepY <- c(seq(2,22,1))

#Set up for parallel processing (still very slow and I'm not sure if it works)
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

#How many variables do we keep on the Y (symptom) side?
tune.spls.cisr$choice.keepY

#What about on the metabolic side?
tune.spls.cisr$choice.keepX


#Now we store the optimal number of x and y components
optimal.keepX <- tune.spls.cisr$choice.keepX
optimal.keepY <- tune.spls.cisr$choice.keepY
optimal.ncomp <- length(optimal.keepX) # extract optimal number of components


# Fit Final Model ======

#Now we can fit the final model

final.spls.cisr <- spls(X = x,
                        Y = y, 
                        ncomp = optimal.ncomp,
                        keepX = optimal.keepX,
                        keepY = optimal.keepY,
                        mode = 'regression')

#Lets save this model given how long it takes to get all the optimised pls
#component counts. Given, however, that the model structure contains the ALSPAC
#data itself, rather than just coefficients etc, I think we would be best off
#saving it into the RDSF


# Save final model =======

write_rds(final.spls.cisr,"./Models/cisr_pls2_model.rds")
#final.spls.cisr = read_rds("./Models/cisr_pls2_model.rds")

# Inspect final model =======

## Loading Plots -----

#Now we can interrogate this optimised model

#Explained variance - not so bad for the cisr, not great for the immuno-metabolic
#datasets

final.spls.cisr$prop_expl_var |>
  as_tibble() |>
  mutate(component = 1:n()) |>
  transmute(component,
            cumulative_metab = cumsum(X),
            cumulative_cisr  = cumsum(Y)) |>
  mutate(across(where(is.double),~round(.x*100,digits = 1))) 

#So we can get up to 22.1% of the variance in the IM data and
#48.5% of the variance in the CISR data


#Plot loadings our way

#Note that we are going to flip the loadings and the exported values for 
#components 1 and 3 because we want to have the main direction of loadings be
#positive rahter than negative. The valence of the PLS components is arbitrary
#so there is no specific problem with doing this, its purely for interpetability
p_loadings_2 = 
  bind_rows(
    final.spls.cisr$loadings$X |>
      as_tibble(rownames = "variable") |>
      mutate(set = "infl"),
    final.spls.cisr$loadings$Y |>
      as_tibble(rownames = "variable") |>
      mutate(set = "cisr"))    |>
  
  #Do the flip on both components
  mutate(
    across(.cols = starts_with("comp"),.fns = ~-.x)) |>
  
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
        axis.text.x.bottom = element_text(angle = -90,vjust = 0.5,size = 8),
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

# - Comp 1 is fatigue/som/sleep and our acute inflammatory variables, plus sex, smoking and bmi
# - Comp 2 is the anxiety, worry and depressive ideations vs crp


# Use Scores as predictors of SMFQ -----

#I guess we model with them

#First get the scores
d_pls2 <- 
  bind_cols(final.spls.cisr$variates$X |>
              as_tibble() |>
              rename_with(~paste("infl",.,sep = "_")),
            final.spls.cisr$variates$Y |>
              as_tibble() |>
              rename_with(~paste("cisr",.,sep = "_"))) |>
  bind_cols(d_complete |> select(id) |> left_join(d_cv,by = "id")) |>
  
  #Do the flip (sign is meaningless in PLS but for interpretation its better to have > score == > symptom)
  mutate(    across(.cols = contains("comp"),.fns = ~-.x)) |>
  
  left_join(d_24 |> select(id,cisr_dep), by = "id") |>
  select(-c(eth,bmi9,bp24,cisr_dep,starts_with("smfq")))



#Get future SMFQ scores, model
d_smfq <- 
  d_pls2 |>
  select(id,contains("_comp"),all_of(covars)) |>
  left_join(d |>
              select(id,smfq_25,smfq_c1,smfq_c2,smfq_c4,smfq_c5),
            by = "id") |>
  pivot_longer(infl_comp1:cisr_comp2,names_to = "comp",values_to = "pls_score") |>
  pivot_longer(starts_with("smfq_"),values_to = "smfq_score",names_to = "smfq_wave") |>
  nest_by(smfq_wave,comp) 

#Do a glmm:
library(lme4)

d_glmm <- 
  d_pls2 |>
  select(id,contains("_comp"),all_of(covars)) |>
  left_join(d |>
              select(id,smfq_25,smfq_c1,smfq_c2,smfq_c4,smfq_c5),
            by = "id") |>
  pivot_longer(infl_comp1:cisr_comp2,names_to = "comp",values_to = "pls_score") |>
  pivot_longer(starts_with("smfq_"),values_to = "smfq_score",names_to = "smfq_wave") |>
  mutate(smfq_wave = case_when(smfq_wave == "smfq_25" ~ 1,
                               smfq_wave == "smfq_c1" ~ 2,
                               smfq_wave == "smfq_c2" ~ 3,
                               smfq_wave == "smfq_c4" ~ 4,
                               smfq_wave == "smfq_c5" ~ 5)) |>
  nest_by(comp) |>
  filter(str_detect(comp,"infl_")) |>
  mutate(m1 = 
           list(
             glmer(smfq_score ~ 1 + pls_score + smfq_wave + pls_score:smfq_wave + (smfq_wave|id),
                   family = poisson(link = "log"),
                   data = data)
           ),
         mt = list(
           avg_comparisons(model = m1,
                           variables = list(pls_score = "sd"),
                           type = "response",
                           re.form=NA,
                           conf_level = 0.95) |>
             as_tibble()
         )) 

#Look at the model output
d_glmm |> select(comp,mt) |> unnest(cols = mt)

#So the average (over all SMFQ waves) effect of a 1 SD in PLS component score is 0.855 SMFQ
#points for comp 1, and 0.393 for comp 2. This varies e.g. 


d_glmm_me <- 
  d_glmm |> 
  select(comp,m1) |> 
  ungroup() |>
  mutate(wave = c(1,5)) |> 
  tidyr::expand(comp,wave) |> 
  left_join(d_glmm |> 
              ungroup() |> 
              select(comp,m1),by = "comp") |> 
  rowwise() |>
  mutate(me = avg_comparisons(model = m1,
                              variables = list(pls_score = "sd"),
                              newdata = datagrid(smfq_wave = wave, grid_type = "counterfactual"),
                              type = "response",
                              re.form=NA,
                              conf_level = 0.95)|>
           as_tibble() |>
           list())


d_glmm_me |>
  dplyr::select(comp,wave,me) |>
  ungroup() |>
  unnest(me) |>
  dplyr::select(comp,wave,contrast,estimate:conf.high)

#Component 1- 

#At wave 1 the difference is:
d_glmm$m1[[1]] |>
  avg_comparisons(model = _,
                  variables = list(pls_score = "sd"),
                  newdata = datagrid(smfq_wave = 1, grid_type = "counterfactual"),
                  type = "response",
                  re.form=NA,
                  conf_level = 0.95)|>
  as_tibble()

#Whereas at wave 5, it is
d_glmm$m1[[1]] |>
  avg_comparisons(model = _,
                  variables = list(pls_score = "sd"),
                  newdata = datagrid(smfq_wave = 5, grid_type = "counterfactual"),
                  type = "response",
                  re.form=NA,
                  conf_level = 0.95)|>
  as_tibble()

#Component 2- 

#At wave 1 the difference is:
d_glmm$m1[[2]] |>
  avg_comparisons(model = _,
                  variables = list(pls_score = "sd"),
                  newdata = datagrid(smfq_wave = 1, grid_type = "counterfactual"),
                  type = "response",
                  re.form=NA,
                  conf_level = 0.95)|>
  as_tibble()

#Whereas at wave 5, it is
d_glmm$m1[[2]] |>
  avg_comparisons(model = _,
                  variables = list(pls_score = "sd"),
                  newdata = datagrid(smfq_wave = 5, grid_type = "counterfactual"),
                  type = "response",
                  re.form=NA,
                  conf_level = 0.95)|>
  as_tibble()



#Plot it
p_glmm1 <- 
  d_glmm$m1[[1]] |>  
  plot_comparisons(variables = list("pls_score" = "sd"),
                   re.form=NA,
                   condition = c("smfq_wave")) +
  geom_hline(yintercept = 0,lty = 2) +
  theme(#panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x.bottom = element_text(angle = -90,vjust = 0.5,size = 6),
    #axis.title.x.bottom = element_text(size = 6),
    axis.title.x.bottom = element_blank(),
    axis.text.y = element_text(size = 6),
    title = element_text(size = 6),
    axis.title = element_blank()) +
  labs(title = "Component 1") +
  coord_cartesian(ylim = c(0,1.6))

p_glmm2 <- 
  d_glmm$m1[[2]] |>  
  plot_comparisons(variables = list("pls_score" = "sd"),
                   re.form=NA,
                   condition = c("smfq_wave"))+
  geom_hline(yintercept = 0,lty = 2)+
  theme(#panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x.bottom = element_text(angle = -90,vjust = 0.5,size = 6),
    #axis.title.x.bottom = element_text(size = 6),
    axis.title.x.bottom = element_blank(),
    axis.text.y = element_text(size = 6),
    title = element_text(size = 6),
    axis.title = element_blank()) +
  labs(title = "Component 2") +
  coord_cartesian(ylim = c(0,1.6))


p_pls2_smfq/(p_glmm1|p_glmm2)

#These results are similar, which is nice

#Could you do a model with both component scores together?

d_glmm_2 <- 
  d_pls2 |>
  select(id,contains("infl_comp"),all_of(covars)) |>
  left_join(d |>
              select(id,smfq_25,smfq_c1,smfq_c2,smfq_c4,smfq_c5),
            by = "id") |>
  pivot_longer(starts_with("smfq_"),values_to = "smfq_score",names_to = "smfq_wave") |>
  mutate(smfq_wave = case_when(smfq_wave == "smfq_25" ~ 1,
                               smfq_wave == "smfq_c1" ~ 2,
                               smfq_wave == "smfq_c2" ~ 3,
                               smfq_wave == "smfq_c4" ~ 4,
                               smfq_wave == "smfq_c5" ~ 5)) 

#This is a strange model. We don't need covariates as they are in the components
m_2 <- 
  glmer(smfq_score ~ 1 + infl_comp1 + infl_comp2 + smfq_wave + infl_comp1:smfq_wave + infl_comp2:smfq_wave + infl_comp1:infl_comp2 +  infl_comp1:infl_comp2:smfq_wave + + (smfq_wave|id),
        family = poisson(link = "log"),
      data = d_glmm_2)


avg_comparisons(model = m_2,
                variables = list(infl_comp1 = "sd",
                                 infl_comp2 = "sd"),
                type = "response",
                re.form=NA,
                conf_level = 0.95) |>
  as_tibble()

m_2 |>  
  plot_comparisons(variables = list("infl_comp2" = "sd"),
                   re.form=NA,
                   condition = c("smfq_wave")) +
  geom_hline(yintercept = 0,lty = 2) +
  theme(#panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x.bottom = element_text(angle = -90,vjust = 0.5,size = 6),
    #axis.title.x.bottom = element_text(size = 6),
    axis.title.x.bottom = element_blank(),
    axis.text.y = element_text(size = 6),
    title = element_text(size = 6),
    axis.title = element_blank()) +
  labs(title = "Component 1") +
  coord_cartesian(ylim = c(0,1.6))


## Make Figure ======

#We need to get some colours for the PLS components and use them consistently
#across the figure. Also nee dot make them different to the cluster colours
#that we used in an earlier figure


#So we also need to get the loading plots...


p_pls2_mega <- 
  (p_loadings_2/(p_glmm1|p_glmm2)) +
  plot_annotation(tag_levels = 'A') + 
  plot_layout(heights = c(2,1))

ggsave(filename = "./Figures/f24_pls2_mega.pdf",
       height = 7,
       width  = 6,
       plot = p_pls2_mega)
