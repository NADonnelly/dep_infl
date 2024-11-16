
# DAG Life -------

#In this script, based on a drawing in Daggity, we will draw the proposed project 
#DAG and explore the minimum adjustment set required for our variables

pacman::p_load(tidyverse,
               dagitty)

source('./Final/0-Common_data.R')


## Cross-sectional DAG  -------

d <- 
  d_24 |> 
  select(id,cisr_dep,il6) |> 
  left_join(d_cv,by = "id") |>
  select(-c(id,smfq_9,smfq_10,smfq_18,md,eth)) |>
  rename(Infl = il6,
         Dep = cisr_dep,
         Sex = sex,
         SES = sc,
         BMI = bmi24,
         smk = smk24,
         Alc = audit_c) |>
  drop_na(Dep,Infl) |> 
  mutate(SES = if_else(SES == "nonmanual",0,1)) |>
  mutate(Sex = if_else(Sex == "Female",0,1)) |>
  mutate(Infl = zscore(Infl))
  

g <-
  dagitty('dag {
                Alc [pos="-0.891,-0.650"]
                BMI [pos="-0.442,-0.650"]
                Dep [outcome,pos="-0.173,-0.103"]
                Infl [exposure,pos="-1.274,-0.103"]
                SES [pos="-1.274,-1.101"]
                Sex [pos="-0.173,-1.101"]
                Smk [pos="-1.274,-0.650"]
                Alc -> BMI
                Alc -> Dep
                Alc -> Infl
                Alc <-> Smk
                BMI -> Dep
                BMI -> Infl
                Infl -> Dep
                SES -> BMI
                SES -> Dep
                SES -> Smk
                Sex -> Alc
                Sex -> Dep
                Sex -> Infl
                Smk -> Dep
                Smk -> Infl
                }')

plot(g)



png(file="./Figures/dep_dag.png",
    units = "cm",
    res = 500,
    height = 10,
    width = 10)
plot(g)
dev.off()

#I would prefer this as a PDF for importing into inkscape
pdf(file="./Figures/dep_dag.pdf",
    height = 3,
    width = 5)
plot(g)
dev.off()


#Ask the DAG which covariates we need to adjust for to shut all
#those pesky backdoor paths - note that we are actually not going
#to include SES as it is not in a path that needs to be blocked
adjustmentSets(g,
               exposure = "Infl",
               outcome  = "Dep")

r <- localTests(g,d |> select(-me) |> drop_na() |> as.matrix())


plotLocalTestResults(r)

#This suggests that alcohol and smoking are not independent as we have implied in 
#our DAG, we caninclude a bidirection link between alcohol and smoking.
