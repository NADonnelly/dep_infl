
# DAG Life -------

#In this script, based on a drawing in Daggity, we will draw the proposed project 
#DAG and explore the minimum adjustment set required for our variables

pacman::p_load(tidyverse,
               dagitty)

source('./Final/0-Common_data.R')


## Cross-sectional DAG  -------

d_dag <- 
  d |> 
  select(id,cisr_dep,il6_f24) |> 
  left_join(d_cv,by = "id") |>
  
  #Make the physical health variable too
  left_join(d_0 |> select(alnqlet,FKDQ2510) |> rename(id = alnqlet),
            by = 'id') |>
  
  #So assign this physical health condition marker to be 1 for all of the listed non-MH
  #conditions
  mutate(ph = case_when(FKDQ2510 %in% c(1:7) ~ 1,
                        FKDQ2510 == 8|FKDQ2510 == 9 ~ 0,
                        FKDQ2510 < 1 ~ NA_integer_)) |>
  select(-FKDQ2510) |>


  select(-c(id,smfq_9,smfq_10,smfq_18,md,eth)) |>
  rename(Infl = il6_f24,
         Dep = cisr_dep,
         Sex = sex,
         SES = sc,
         BMI = bmi24,
         Smk = smk24,
         Alc = audit_c,
         CD  = ph) |>
  drop_na(Dep,Infl) |> 
  mutate(SES = if_else(SES == "nonmanual",0,1)) |>
  mutate(Sex = if_else(Sex == "Female",0,1)) |>
  mutate(Dep = as.numeric(Dep)-1) |>
  mutate(Infl = zscore(Infl)) |>
  
  select(-c(me,bmi9,bp24)) |>
  drop_na()



g <-
  dagitty('dag {
                Alc [pos="-0.891,-0.650"]
                BMI [pos="-0.442,-0.650"]
                Dep [outcome,pos="-0.173,-0.103"]
                Infl [exposure,pos="-1.274,-0.103"]
                SES [pos="-1.274,-1.101"]
                Sex [pos="-0.173,-1.101"]
                Smk [pos="-1.274,-0.650"]
                CD [pos="-0.65,-0.3"]
                Alc -> BMI
                Alc -> Dep
                Alc -> Infl
                Alc -> CD
                Alc <-> Smk
                BMI -> Dep
                BMI -> Infl
                BMI -> CD
                Infl -> Dep
                SES -> BMI
                SES -> Dep
                SES -> Smk
                SES -> CD
                Sex -> Alc
                Sex -> Dep
                Sex -> Infl
                Sex -> CD
                Smk -> Dep
                Smk -> Infl
                Smk -> CD
                CD -> Dep
                CD -> Infl
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


r <- 
  d_dag |> 
  # filter(if_any(everything(), ~is.infinite(.x)))
  # filter(if_any(everything(), ~is.na(.x)))
  # as.matrix() |>
  dagitty::localTests(x = g,
                      data = _)

r <- localTests(g,)


plotLocalTestResults(r)

#This suggests that alcohol and smoking are not independent as we have implied in 
#our DAG, therefore we include a bidirection link between alcohol and smoking

