
# DAG Life -------

#In this script, based on a drawing in Daggity, we will draw the proposed project 
#DAG and explore the minimum adjustment set required for our variables

pacman::p_load(tidyverse,
               dagitty)

source('./Final/0-Common_data.R')


## Cross-sectional DAG  -------

d_dag <- 
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

r <- localTests(g,d_dag |> select(-me) |> drop_na() |> as.matrix())


plotLocalTestResults(r)

#This suggests that alcohol and smoking are not independent as we have implied in 
#our DAG

# I suspect there is an unobserved latent variable that is "substance misuse"
#here, or we could include a bidirection link between alcohol and smoking?


#Equivalent DAGS

ec <- equivalentDAGs(g)

dagitty(ec[[1]]) |> plot()
dagitty(ec[[2]]) |> plot()



## Longitudinal DAG ======

#We ought to also make a DAG for doing longitudinal analysis 
d_l = 
  d_cv |>
  select(id, audit_c,smk24,sex,sc,md,smfq_9,bmi24) |>
  rename(alc24 = audit_c,
         ses = sc,
         matdep = md,
         smfq9 = smfq_9) |>
  full_join(
    d_0 |>
      select(alnqlet,f9ms026a) |>
      rename(id = alnqlet,
             bmi9 = f9ms026a),
    by = 'id') |>
  right_join(d_9 |>
               select(id,il6) |>
               rename(infl9 = il6),
             by = 'id') |>
  right_join(d_24 |>
               select(id,cisr_dep,il6) |>
               rename(infl24 = il6,
                      dep24 = cisr_dep),
             by = 'id') |>
  
  #I recall the can only have numerical variables 
  mutate(sex = if_else(sex == "Female",0,1),
         ses = case_when(ses == "I" ~ 6,
                         ses == "II" ~ 5,
                         ses == "III non-manual" ~ 4,
                         ses == "III manual" ~ 3,
                         ses == "IV" ~ 2,
                         ses == "V" ~ 1,
                         TRUE ~ NA_integer_)) |>
  drop_na()


g <-
  dagitty('dag {
                alc24 [pos="0.078,-1.15"]
                bmi24 [pos="1.25,-0.546"]
                bmi9 [pos="-2.25,-1.520"]
                dep24 [outcome,pos="1.25,0.6"]
                infl24 [pos="0.503,0.01"]
                infl9 [exposure,pos="-1.114,0.6"]
                matdep [pos="-1.917,-0.25"]
                ses [pos="-2.25,0.1"]
                sex [pos="-1.268,1.354"]
                smfq9 [pos="-1.116,-0.25"]
                smk24 [pos="0.722,-1.15"]
                alc24 -> bmi24
                alc24 -> dep24
                alc24 -> infl24
                alc24 <-> smk24
                bmi24 -> dep24
                bmi24 -> infl24
                bmi9 -> bmi24
                bmi9 -> infl9
                bmi9 -> smk24
                infl24 -> dep24
                infl9 -> bmi24
                infl9 -> dep24
                infl9 -> infl24
                infl9 -> smfq9
                matdep -> alc24
                matdep -> dep24
                matdep -> smfq9
                ses -> bmi24
                ses -> bmi9
                ses -> dep24
                ses -> infl9
                ses -> matdep
                ses -> smk24
                sex -> alc24
                sex -> bmi9
                sex -> dep24
                sex -> infl9
                smfq9 -> alc24
                smfq9 -> dep24
                smfq9 -> smk24
                smk24 -> dep24
                smk24 -> infl24
                }
')

plot(g)



png(file="./Figures/dep_dag_longitudinal.png",
    units = "cm",
    res = 500,
    height = 20,
    width = 20)
plot(g)
dev.off()



#Ask the DAG which covariates we need to adjust for to shut all
#those pesky backdoor paths - note that we are actually not going
#to include SES as it is not in a path that needs to be blocked
adjustmentSets(g,
               exposure = "infl9",
               outcome  = "dep24")

#So a simple set is bmi9, ses and sex

r <- localTests(g,d_l |> select(-id) |> drop_na() |> as.matrix())

plotLocalTestResults(r)

#This suggests many complicated things

# - BMI at age 24 is apparently not conditionally independent of infl9 after 
#adjustment for bmi9, ses and sex - so lets put in a link between infl9 and
#bmi24

# - Alcohol at 24 is not independent of smfq9 at 9 after infl9 and matdep
# apparently higher smfq is linked to lower alcohol

lm(alc24 ~ smfq9 + infl9 + matdep,data = d_l) |>
  tidy()

lm(alc24 ~ smfq9 + infl9 + matdep + ses,data = d_l) |>
  tidy(conf.int = T)

glm(smk24 ~ bmi9 + infl9 + ses + smfq9, data = d_l, family = binomial(link = 'probit')) |>
  broom::tidy(conf.int = T, exponentiate = T)


#OK we have now added the missing links and we get the adjustment sets we need

#Equivalent DAGS

ec <- equivalentDAGs(g)

dagitty(ec[[1]]) |> plot()
dagitty(ec[[2]]) |> plot()



