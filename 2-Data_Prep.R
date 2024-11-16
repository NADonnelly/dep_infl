# Introduction ------

#We are going to load the data, assemble a single dataset of all the variables 
#we are going to use for the rest of the project. Based on comments from a previous round of reviews, 
#we are not doing a test/train split but we will use nested CV for the machine learning later

pacman::p_load(tidyverse,
               tidymodels,
               patchwork,
               crayon)

tidymodels_prefer()

conflicted::conflicts_prefer(crayon::`%+%`)

source('./Final/0-Common_data.R')

#Do some clearing up
rm(list = c("wdn","valid_vars_f24","valid_vars_f9","v_0"))


# Data Preparation ------

## CISR Data ======

### CISR Diagnosis -----

#CISR diagnostic data
d_cisr = 
  d_0 |> 
  select(alnqlet,starts_with("FKDQ")) |> 
  rename(id  = alnqlet) |>
  select(id, 
         FKDQ1000:FKDQ1020)  |>
  mutate(across(FKDQ1000:FKDQ1020,~if_else(. < 0,NA_integer_,.)))|>
  transmute(id,
            
            #You need the depression diagnostic variables here too
            cisr_dep     = FKDQ1000,
            cisr_dep_mod = FKDQ1010,
            cisr_dep_sev = FKDQ1020) 

### CISR Subscales ------

#Build the CISR subscale scoes
d_cisr_scales = 
  d_0 |> 
  select(alnqlet,starts_with("FKDQ")) |> 
  rename(id  = alnqlet) |>
  
  #Get rid of missing data on the basis that if you have a depression diagnosis
  #you must have completed the CISR
  filter(FKDQ1000 >=0) |>
  
  select(id, FKDQ2000:FKDQ7000) |>
  
  transmute(id,
            
            #Somatic symptom score
            
            #Somatic symptom 1
            som1 = if_else(FKDQ2620 == 3,1,0) + if_else(FKDQ2720 == 3,1,0),
            som1 = if_else(som1 > 0,1,0),
            
            #Somatic symptom 2
            som2 = if_else(FKDQ2630 == 2,1,0) + if_else(FKDQ2730 == 2,1,0),
            som2 = if_else(som2 > 0,1,0),
            
            #Somatic symptom 3
            som3 = if_else(FKDQ2640 >= 3,1,0) + if_else(FKDQ2740 >= 3,1,0),
            som3 = if_else(som3 > 0,1,0),
            
            #Somatic symptom 4
            som4 = if_else(FKDQ2650 == 2,1,0) + if_else(FKDQ2750 == 2,1,0),
            som4 = if_else(som4 > 0,1,0),  
            
            #Total Somatic symptom score
            som = som1+som2+som3+som4,
            
            
            #Fatigue symptom score
            
            #Fatigue symptom 1
            ftg1 = if_else(FKDQ3020 == 3,1,0) + if_else(FKDQ3120 == 3,1,0),
            ftg1 = if_else(ftg1 > 0,1,0),
            
            #Fatigue symptom 2
            ftg2 = if_else(FKDQ3030 == 2,1,0) + if_else(FKDQ3130 == 2,1,0),
            ftg2 = if_else(ftg2 > 0,1,0),
            
            #Fatigue symptom 3
            ftg3 = if_else(FKDQ3040 == 2,1,0) + if_else(FKDQ3140 == 2,1,0),
            ftg3 = if_else(ftg3 > 0,1,0),
            
            #Fatigue symptom 4
            ftg4 = if_else(FKDQ3050 == 2,1,0) + if_else(FKDQ3150 == 2,1,0),
            ftg4 = if_else(ftg4 > 0,1,0),
            
            #Total Fatigue symptom score
            ftg = ftg1+ftg2+ftg3+ftg4,
            
            
            #Concentration and forgetfulness score
            
            #Concentration symptom 1
            con1 = if_else(FKDQ3520 == 3,1,0),
            
            #Concentration symptom 2
            con2 = if_else(FKDQ3530 == 2,1,0),
            
            #Concentration symptom 3
            con3 = if_else(FKDQ3540 == 2,1,0),
            
            #Concentration symptom 4
            con4 = if_else(FKDQ3560 == 2,1,0),
            
            #Total Concentration symptom score
            con = con1+con2+con3+con4,
            
            
            #Sleep problems score
            
            #To get this onto a 0-4 scale, we will unite the hypersomnia
            #and insomnia items, rather than treating them separately
            
            
            #Sleep symptom 1: nights out of past 7 with problems with sleep
            slp1 = if_else(FKDQ4010 == 3,1,0), 
            
            
            #Sleep symptom 2: nights/7 with either > 3 hours trying to get to 
            #sleep OR > 3 hours excess sleep
            slp2 = if_else(FKDQ4030 == 3,1,0) + if_else(FKDQ4130 == 3,1,0),
            slp2 = if_else(slp2 > 0,1,0),
            
            #Sleep symptom 3: night w/least sleep - how long trying to sleep
            slp3 = case_when(FKDQ4020 == 2 ~ 1,
                             FKDQ4020 == 3 ~ 2,
                             FKDQ4020 == 4 ~ 2,
                             TRUE ~ 0),
            
            #Sleep symptom 4: night w/most sleep - how much longer than usual
            slp4 = case_when(FKDQ4120 == 2 ~ 1,
                             FKDQ4120 == 3 ~ 2,
                             FKDQ4120 == 4 ~ 2,
                             TRUE ~ 0),
            
            #Now we do some logic to merge
            slp5 = case_when(slp3 == 0 & slp4 == 0 ~ 0,
                             slp3 == 1 & slp4 == 0 ~ 1,
                             slp3 == 0 & slp4 == 1 ~ 1,
                             slp3 == 1 & slp4 == 1 ~ 2,
                             slp3 == 2 | slp4 == 2 ~ 2),
            
            #Total Sleep Problem score
            slp = slp1 + slp2 + slp5,
            
            
            #Irritability
            
            #irritability symptom 1
            irt1 = if_else(FKDQ4520 == 3,1,0),
            
            #irritability symptom 2
            irt2 = if_else(FKDQ4530 == 2,1,0),
            
            #irritability symptom 3
            irt3 = if_else(FKDQ4540 >= 2,1,0),
            
            #irritability symptom 4
            irt4 = if_else(FKDQ4550 == 3,1,0),
            
            #Total irritability score
            irt = irt1 + irt2 + irt3 + irt4,
            
            
            #Depression
            
            #depression symptom 1 - anhedonia
            dep1 = if_else(FKDQ5030 >= 2,1,0),
            
            #depression symptom 2 - feeling sad days/7
            dep2 = if_else(FKDQ5040 == 3,1,0),
            
            #depression symptom 3 sad/anhedonia >3 hours/day
            dep3 = if_else(FKDQ5050 == 2,1,0),
            
            #depression symptom 4 nothing cheered up
            dep4 = if_else(FKDQ5070 == 3,1,0),
            
            #Total depression score
            dep = dep1 + dep2 + dep3 + dep4,
            
            
            #Depressive ideas (scores 0-5)
            
            #depressive ideas symptom 1 - Guilt
            did1 = if_else(FKDQ5140 >= 3,1,0),
            
            #depressive ideas symptom 2 - Inadequacy 
            did2 = if_else(FKDQ5150 == 2,1,0),
            
            #depressive ideas symptom 3 - Hopelessness
            did3 = if_else(FKDQ5160 == 2,1,0),
            
            #depressive ideas symptom 4 - Life not worth living
            did4 = if_else(FKDQ5170 >= 2,1,0),
            
            #depressive ideas symptom 5 - Suicidal ideation
            did5 = if_else(FKDQ5180 == 2,1,0),            
            
            #Total depressive ideas score
            did = did1 + did2 + did3 + did4 + did5,
            
            
            #Worry
            
            #worry symptom 1
            wor1 = if_else(FKDQ6020 == 3,1,0),
            
            #worry symptom 2
            wor2 = if_else(FKDQ6030 == 2,1,0),
            
            #worry symptom 3
            wor3 = if_else(FKDQ6040 >= 3,1,0),
            
            #worry symptom 4
            wor4 = if_else(FKDQ6050 == 2,1,0),
            
            #Total worry score
            wor = wor1 + wor2 + wor3 + wor4,
            
            
            #anxiety
            
            #anxiety symptom 1
            anx1 = if_else(FKDQ6300 == 3,1,0),
            
            #anxiety symptom 2
            anx2 = if_else(FKDQ6310 >= 3,1,0),
            
            #anxiety symptom 3
            anx3 = 
              if_else(FKDQ6320 == 2,1,0) + 
              if_else(FKDQ6330 == 2,1,0) +
              if_else(FKDQ6340 == 2,1,0) +
              if_else(FKDQ6350 == 2,1,0) + 
              if_else(FKDQ6360 == 2,1,0) +
              if_else(FKDQ6370 == 2,1,0) +
              if_else(FKDQ6380 == 2,1,0) +
              if_else(FKDQ6390 == 2,1,0),
            anx3 = if_else(anx3 > 0,1,0),
            
            #anxiety symptom 4
            anx4 = if_else(FKDQ6400 == 2,1,0),
            
            #Total anxiety score
            anx = anx1 + anx2 + anx3 + anx4,
            
            
            #phobias
            
            #phobias symptom 1
            pho1 = if_else(FKDQ6520 == 3,1,0),
            
            #phobias symptom 2
            pho2 =               
              if_else(FKDQ6530 == 2,1,0) + 
              if_else(FKDQ6540 == 2,1,0) +
              if_else(FKDQ6550 == 2,1,0) +
              if_else(FKDQ6560 == 2,1,0) + 
              if_else(FKDQ6570 == 2,1,0) +
              if_else(FKDQ6580 == 2,1,0) +
              if_else(FKDQ6590 == 2,1,0) +
              if_else(FKDQ6600 == 2,1,0),
            
            pho2 = if_else(pho2 > 0,1,0),
            
            #phobias symptom 3
            pho3 = case_when(FKDQ6620 == 2 ~ 1,
                             FKDQ6620 == 3 ~ 2,
                             TRUE ~ 0),
            
            #Total phobias symptom score
            pho = pho1+pho2+pho3,
            
            
            #panic
            
            #panic symptom 1
            pan1 = case_when(FKDQ6710 == 2 ~ 1,
                             FKDQ6710 == 3 ~ 2,
                             TRUE ~ 0),
            
            #panic symptom 2
            pan2 = if_else(FKDQ6720 >= 2,1,0),
            
            #panic symptom 3
            pan3 = if_else(FKDQ6730 == 2,1,0),
            
            #Total panic score
            pan = pan1 + pan2 + pan3
            
  ) |>
  
  select(id,som,ftg,con,slp,irt,dep,did,wor,anx,pho,pan)


### CISR Domain Scales -----

#using CIS-R symptom data create specific symptom scores (somatic, 
#psychological, atypical). These will be created by summing individual symptoms. 
#In Ye et al EClinMed, we created somatic and psychological scores using PHQ-9 
#items. For atypical, you will find papers by Yuri Milaneshi, Brenda Penninx and 
#Femke Lamers from Amsterdam describing how this is computed. 

#OK so in Milaneschi https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8873022/ 
#they make a "somatic" score from PHQ items like this - 

#DEPRESSION Somatic = D.3 + D.4 + D.5 + D.7 + D.8

#3. Trouble sleeping, or sleeping too much
#4. Recent feelings of tiredness or low energy
#5. Recent poor appetite or overeating
#7. Recent trouble concentration
#8. Recent changes in moving/speaking

#What are the nearest CISR variableS?
#For 3, FKDQ4000 (troubling getting to sleep, binary), FKDQ4100 (sleeping more 
#       then normal, 1-3)
#For 4, FKDQ3000 (getting tired, binary), FKDQ3100 (lacking energy, binary)
#For 5, FKDQ2000 (loss of appetite, binary), FKDQ2100 (increased appetite, binary)
#For 7, FKDQ3500 (problems concentrating, binary)
#For 8, This is a problem because we only have contingent items which not all 
#       participants answered and zeroing people who didn't get that question 
#       introduces bias. The items are FKDQ5120 and FKDQ5130. 

# They make a "psychological" score from -

#DEPRESSION Psychological = D.1 + D.2 + D.6 + D.9

#1. Recent lack of interest or pleasure in doing things
#2. Recent feelings of depression
#6. Recent feelings of inadequacy
#9. Recent thoughts of suicide or self-harm

#What are the nearest CISR variables?

#For 1, FKDQ5020 (able to enjoy things, 1-3)
#For 2, FKDQ5000 (spell of feeling sad, binary)
#For 6, FKDQ5150 (contingent variable, not as good as others, binary)
#For 9, FKDQ5180 (contingent variable, thoughts of harming themselves, binary)

#And they also use a summary of the GAD7 items, we can create a version
#of this from the CISR too

#A1 - FKDQ6100 (anxious/nervous,binary)
#A2 - FKDQ6000 (worried more than needed to about things, 1-3)
#A3 - FKDQ6030 (contingent variable, worrying too much, binary)
#A4 - FKDQ6110 (muscles tense or could not relax, 1-3)
#A5 - FKDQ5120 (restless, but a contingent variable from the sad section, binary)
#A6 - FKDQ4500 (irritable, binary)
#A7 - FKDQ6200 (anxious about situations where no real danger, binary)


# For atypical, we have previously done a combo of:
#increased appetite, increased weight, hypersomnia and low energy

#Lets do it
d_symptom_score = 
  d_0 |> 
  select(alnqlet,starts_with("FKDQ")) |> 
  rename(id  = alnqlet) |>
  
  #Get rid of missing data on the basis that if you have a depression diagnosis
  #you must have completed the CISR
  filter(FKDQ1000 >=0) |>
  select(id, FKDQ2000:FKDQ7000) |>
  
  #Make all the <0 codes into NA
  mutate(across(FKDQ2000:FKDQ7000,~if_else(. < 0,NA_integer_,.))) |>
  
  #Now lets turn any 1-3 variables into binary and make sure we are coding the
  #no symptom response as zero...
  mutate(FKDQ2000 = case_when(FKDQ2000 == 1 ~ 0,
                              FKDQ2000 == 2 ~ 1),
         FKDQ2100 = case_when(FKDQ2100 == 1 ~ 0,
                              FKDQ2100 == 2 ~ 1),
         FKDQ2110 = case_when(FKDQ2110 == 1 ~ 0,
                              FKDQ2110 == 2 ~ 1),
         FKDQ2120 = case_when(FKDQ2120 == 1 ~ 0,
                              FKDQ2120 == 2 ~ 1,
                              FKDQ2120 == 3 ~ 1),
         FKDQ3000 = case_when(FKDQ3000 == 1 ~ 0,
                              FKDQ3000 == 2 ~ 1),
         FKDQ3100 = case_when(FKDQ3100 == 1 ~ 0,
                              FKDQ3100 == 2 ~ 1),
         FKDQ3500 = case_when(FKDQ3500 == 1 ~ 0,
                              FKDQ3500 == 2 ~ 1),
         FKDQ4000 = case_when(FKDQ4000 == 1 ~ 0,
                              FKDQ4000 == 2 ~ 1),
         FKDQ4100 = case_when(FKDQ4100 == 1 ~ 0,
                              FKDQ4100 == 2 ~ 0,
                              FKDQ4100 == 3 ~ 1),
         FKDQ4500 = case_when(FKDQ4500 == 1 ~ 0,
                              FKDQ4500 == 2 ~ 1),
         FKDQ5000 = case_when(FKDQ5000 == 1 ~ 0,
                              FKDQ5000 == 2 ~ 1),
         FKDQ5020 = case_when(FKDQ5020 == 1 ~ 0,
                              FKDQ5020 == 2 ~ 1,
                              FKDQ5020 == 3 ~ 1),
         FKDQ5120 = case_when(FKDQ5120 == 1 ~ 0,
                              FKDQ5120 == 2 ~ 1),
         FKDQ5150 = case_when(FKDQ5150 == 1 ~ 0,
                              FKDQ5150 == 2 ~ 1),
         FKDQ5130 = case_when(FKDQ5130 == 1 ~ 0,
                              FKDQ5130 == 2 ~ 1),
         FKDQ5180 = case_when(FKDQ5180 == 1 ~ 0,
                              FKDQ5180 == 2 ~ 1),
         FKDQ6000 = case_when(FKDQ6000 == 1 ~ 0,
                              FKDQ6000 == 2 ~ 1,
                              FKDQ6000 == 3 ~ 1),
         FKDQ6030 = case_when(FKDQ6030 == 1 ~ 0,
                              FKDQ6030 == 2 ~ 1), 
         FKDQ6100 = case_when(FKDQ6100 == 1 ~ 0,
                              FKDQ6100 == 2 ~ 1), 
         FKDQ6110 = case_when(FKDQ6110 == 1 ~ 0,
                              FKDQ6110 == 2 ~ 1,
                              FKDQ6110 == 3 ~ 1),
         FKDQ6200 = case_when(FKDQ6200 == 1 ~ 0,
                              FKDQ6200 == 2 ~ 1)
         
  ) |>
  
  #Make a common weight gain variable from the ALSPAC sex-split variables
  mutate(FKDQ2115 = rowSums(across(FKDQ2110:FKDQ2120),na.rm = T))|> 
  
  #Build the scores using the variables above
  transmute(id,
            atypical      = rowSums(across(c(FKDQ2100,FKDQ2115,FKDQ4100,FKDQ3100)),na.rm = T),
            somatic       = rowSums(across(c(FKDQ4000,FKDQ4100,
                                             FKDQ3000,FKDQ3100,
                                             FKDQ2000,FKDQ2100,
                                             FKDQ3500,
                                             FKDQ5120,FKDQ5130)),na.rm = T),
            psychological = rowSums(across(c(FKDQ5020,FKDQ5000,FKDQ5150,FKDQ5180)),na.rm = T),
            anxiety       = rowSums(across(c(FKDQ6100,FKDQ6000,FKDQ6030,FKDQ6110,
                                             FKDQ5120,FKDQ4500,FKDQ6200)),na.rm = T))



## SMFQ Data ========== 

#First lets select the SMFQ data

d_smfq = 
  d_0 |> 
  select(alnqlet,
         
         #KU
         ku673b,  # Age 9 SMFQ total score (adjusted for missing data) - parent rated
         
         
         #F10
         fddp130, # Age 10 SMFQ- clinic based with the child posting things to indicate repsonses
         
         
         #KW9
         kw6100b, # Age 11.5 SMFQ - parent completed
         
         
         #TF1 - age 12.5  - no total score - we can make it from
         ff6500, # - Feel miserable or unhappy 
         ff6502, # - Didn't enjoy anything
         ff6503, # - So tired I sat around and did nothing
         ff6504, # - Felt very restless
         ff6505, # - Felt I was no good anymore
         ff6506, # - Cried a lot
         ff6508, # - Found it hard to think properly or concentrate
         ff6509, # - Hated myself
         ff6511, # - felt I was a bad person
         ff6512, # - Felt lonely
         ff6513, # - Thought nobody really loved me
         ff6514, # - Thought I would never be as good as other kids
         ff6515, # - Did everything wrong
         
         
         #TA - age 13 - parent completed (also has what looks like SAPAS?)
         ta5020, # - Feel miserable or unhappy 
         ta5021, # - Didn't enjoy anything
         ta5022, # - So tired I sat around and did nothing
         ta5023, # - Felt very restless
         ta5024, # - Felt I was no good anymore
         ta5025, # - Cried a lot
         ta5026, # - Found it hard to think properly or concentrate
         ta5027, # - Hated myself
         ta5028, # - felt I was a bad person
         ta5029, # - Felt lonely
         ta5030, # - Thought nobody really loved me
         ta5031, # - Thought I would never be as good as other kids
         ta5032, # - Did everything wrong
         
         
         #TF2
         fg7226, #TF2 i.e. clinic based Age 13.5 SMFQ total score - 
         
         
         #CCS - age 16.5 - no total score - we can make it from - 
         ccs4500, # - Feel miserable or unhappy 
         ccs4502, # - Didn't enjoy anything
         ccs4503, # - So tired I sat around and did nothing
         ccs4504, # - Felt very restless
         ccs4505, # - Felt I was no good anymore
         ccs4506, # - Cried a lot
         ccs4508, # - Found it hard to think properly or concentrate
         ccs4509, # - Hated myself
         ccs4511, # - felt I was a bad person
         ccs4512, # - Felt lonely
         ccs4513, # - Thought nobody really loved me
         ccs4514, # - Thought I would never be as good as other kids
         ccs4515, # - Did everything wrong
         
         
         #CCXD - age 17.5
         CCXD917, # - sum total score
         CCXD918, # - flag for respondents just putting "somewwhat true" for everyhitng - 
         # - relevant codes are 1,3,5,7
         
         
         #CCT - age 18
         cct2715, # -  sum total score
         
         
         #YPA - age 21 - need to make a total score - from
         YPA2000, # - Feel miserable or unhappy 
         YPA2010, # - Didn't enjoy anything
         YPA2020, # - So tired I sat around and did nothing
         YPA2030, # - Felt very restless
         YPA2040, # - Felt I was no good anymore
         YPA2050, # - Cried a lot
         YPA2060, # - Found it hard to think properly or concentrate
         YPA2070, # - Hated myself
         YPA2080, # - felt I was a bad person
         YPA2090, # - Felt lonely
         YPA2100, # - Thought nobody really loved me
         YPA2110, # - Thought I would never be as good as other kids
         YPA2120, # - Did everything wrong
         
         
         #YPB - age 22 - total score
         YPB5180,
         
         
         #YPC - age 23 
         YPC1650, # - Feel miserable or unhappy 
         YPC1651, # - Didn't enjoy anything
         YPC1653, # - So tired I sat around and did nothing
         YPC1654, # - Felt very restless
         YPC1655, # - Felt I was no good anymore
         YPC1656, # - Cried a lot
         YPC1658, # - Found it hard to think properly or concentrate
         YPC1659, # - Hated myself
         YPC1660, # - felt I was a bad person
         YPC1662, # - Felt lonely
         YPC1663, # - Thought nobody really loved me
         YPC1665, # - Thought I would never be as good as other kids
         YPC1667, # - Did everything wrong
         
         
         #YPE - age 25
         YPE4080, # - Feel miserable or unhappy 
         YPE4082, # - Didn't enjoy anything
         YPE4083, # - So tired I sat around and did nothing
         YPE4084, # - Felt very restless
         YPE4085, # - Felt I was no good anymore
         YPE4086, # - Cried a lot
         YPE4088, # - Found it hard to think properly or concentrate
         YPE4089, # - Hated myself
         YPE4091, # - felt I was a bad person
         YPE4092, # - Felt lonely
         YPE4093, # - Thought nobody really loved me
         YPE4094, # - Thought I would never be as good as other kids
         YPE4095, # - Did everything wrong
         
         #COVID era questionnaires - 1, 2, 4, 5
         #These were done around age 29
         
         #COVID1
         covid1yp_4050, # - Feel miserable or unhappy 
         covid1yp_4051, # - Didn't enjoy anything
         covid1yp_4052, # - So tired I sat around and did nothing
         covid1yp_4053, # - Felt very restless
         covid1yp_4054, # - Felt I was no good anymore
         covid1yp_4055, # - Cried a lot
         covid1yp_4056, # - Found it hard to think properly or concentrate
         covid1yp_4057, # - Hated myself
         covid1yp_4058, # - felt I was a bad person
         covid1yp_4059, # - Felt lonely
         covid1yp_4060, # - Thought nobody really loved me
         covid1yp_4061, # - Thought I would never be as good as other kids
         covid1yp_4062, # - Did everything wrong
         
         #COVID2
         covid2yp_4050, # - Feel miserable or unhappy 
         covid2yp_4051, # - Didn't enjoy anything
         covid2yp_4052, # - So tired I sat around and did nothing
         covid2yp_4053, # - Felt very restless
         covid2yp_4054, # - Felt I was no good anymore
         covid2yp_4055, # - Cried a lot
         covid2yp_4056, # - Found it hard to think properly or concentrate
         covid2yp_4057, # - Hated myself
         covid2yp_4058, # - felt I was a bad person
         covid2yp_4059, # - Felt lonely
         covid2yp_4060, # - Thought nobody really loved me
         covid2yp_4061, # - Thought I would never be as good as other kids
         covid2yp_4062, # - Did everything wrong
         
         #COVID4
         covid4yp_4050, # - Feel miserable or unhappy 
         covid4yp_4051, # - Didn't enjoy anything
         covid4yp_4052, # - So tired I sat around and did nothing
         covid4yp_4053, # - Felt very restless
         covid4yp_4054, # - Felt I was no good anymore
         covid4yp_4055, # - Cried a lot
         covid4yp_4056, # - Found it hard to think properly or concentrate
         covid4yp_4057, # - Hated myself
         covid4yp_4058, # - felt I was a bad person
         covid4yp_4059, # - Felt lonely
         covid4yp_4060, # - Thought nobody really loved me
         covid4yp_4061, # - Thought I would never be as good as other kids
         covid4yp_4062, # - Did everything wrong
         
         #COVID5
         covid5yp_4050, # - Feel miserable or unhappy 
         covid5yp_4051, # - Didn't enjoy anything
         covid5yp_4052, # - So tired I sat around and did nothing
         covid5yp_4053, # - Felt very restless
         covid5yp_4054, # - Felt I was no good anymore
         covid5yp_4055, # - Cried a lot
         covid5yp_4056, # - Found it hard to think properly or concentrate
         covid5yp_4057, # - Hated myself
         covid5yp_4058, # - felt I was a bad person
         covid5yp_4059, # - Felt lonely
         covid5yp_4060, # - Thought nobody really loved me
         covid5yp_4061, # - Thought I would never be as good as other kids
         covid5yp_4062  # - Did everything wrong
         
  ) |>
  rename(id  = alnqlet) |>
  
  #Do some recoding specifically for the age 17.5 data where some participants 
  #put the middle response for everything - if flagged, NA them 
  mutate(CCXD917 = case_when(CCXD918 == 1 ~ NA_integer_,
                             CCXD918 == 3 ~ NA_integer_,
                             CCXD918 == 5 ~ NA_integer_,
                             CCXD918 == 7 ~ NA_integer_,
                             TRUE ~ CCXD917)) |>
  
  mutate(across(-id,
                ~case_when( . < 0 ~ NA_integer_ ,
                            TRUE ~ .))) |>
  
  #Recode variables where the coding scheme seems to be reversed
  mutate(across(starts_with("ff") | starts_with("ta")|starts_with("ccs")| starts_with("YPA"),
                ~case_when( . < 0 ~ NA_integer_ ,
                            . == 3 ~ 0,
                            . == 2 ~ 1,
                            . == 1 ~ 2,
                            TRUE ~ .))) |>
  
  mutate(across(starts_with("YPC"),
                ~case_when(. == 1 ~ 0,
                           . == 2 ~ 1,
                           . == 3 ~ 2,
                           TRUE ~ .))) |>
  
  #Make sum scores
  mutate(smfq_12.5 = rowSums(across(starts_with("ff"))),
         smfq_13   = rowSums(across(starts_with("ta"))),
         smfq_16.5 = rowSums(across(starts_with("ccs"))),
         smfq_21   = rowSums(across(starts_with("YPA"))),
         smfq_23   = rowSums(across(starts_with("YPC"))),
         smfq_25   = rowSums(across(starts_with("YPE"))),
         smfq_c1   = rowSums(across(starts_with("covid1yp_"))),
         smfq_c2   = rowSums(across(starts_with("covid2yp_"))),
         smfq_c4   = rowSums(across(starts_with("covid4yp_"))),
         smfq_c5   = rowSums(across(starts_with("covid5yp_")))) |>
  
  #Rename our variables that came as sum scores already into a common format
  rename(smfq_9    = ku673b,
         smfq_10   = fddp130,
         smfq_11.5 = kw6100b,
         smfq_13.5 = fg7226,
         smfq_17.5 = CCXD917,
         smfq_18   = cct2715,
         smfq_22   = YPB5180) |>

  select(id, 
         smfq_9,smfq_10,smfq_11.5,smfq_12.5,smfq_13, smfq_13.5,
         smfq_16.5,smfq_17.5,smfq_18,smfq_21,smfq_22,smfq_23,
         smfq_25,smfq_c1,smfq_c2,smfq_c4,smfq_c5)


# d_smfq |>
#   select(-c(smfq_9,smfq_11.5,smfq_13)) |>
#   pivot_longer(-id,names_to = "age",values_to = "smfq_total") |>
#   drop_na() |>
#   mutate(age = str_remove(age,"smfq_")) |>
#   
#   ggplot(aes(x = age,y = smfq_total)) +
#   geom_boxplot()

## SAPAS ------

#We also include negative symptoms in keeping with our leptin paper

#SAPAS is 8x binary items, so is binomial from 8 which is simple - 
d_sapas <- 
  d_0 |> 
  select(alnqlet,FKDQ1000,starts_with("FKPE1")) |> 
  rename(id  = alnqlet) |>
  
  #Get rid of missing data on the basis that if you have a depression diagnosis
  #you must have completed the CISR
  filter(FKDQ1000 >=0) |>
  select(id, starts_with("FKPE1")) |>

  
  #Make all the <0 codes into NA
  mutate(across(-id,~if_else(. < 0,NA_integer_,.))) |>
  
  #Do some reverse coding so that higher scores == more personality difficulty
  mutate(FKPE1020 = case_when(FKPE1020 == 0 ~ 1,
                              FKPE1020 == 1 ~ 0,
                              TRUE ~ FKPE1020))

d_sapas <- 
  d_sapas |> 
  transmute(id,sapas_total = rowSums(across(where(is.numeric))))


## CAPE ------

d_cape <- 
  d_0 |> 
  select(alnqlet,FKDQ1000,starts_with("FKPE2")) |> 
  rename(id  = alnqlet) |>
  filter(FKDQ1000 >=0) |>
  select(id, starts_with("FKPE2")) |>
  
  #Make all the <0 codes into NA
  mutate(across(-id,~if_else(. < 0,NA_integer_,.))) |>
  
  #Reverse code
  mutate(across(-id,~4 - .x)) |> 
  transmute(id,cape_total = rowSums(across(where(is.numeric))))



## Auxillary Variables ------

#Variables that are available from earlier rounds which might contain useful 
#data for imputation
d_aux = 
  d_0 |> 
  select(
    alnqlet, # ID variable
    
    #BMI variables
    f7ms026a,    #F7 BMI
    f9ms026a,    #F9 BMI
    fdms026a,    #F10 BMI
    fems026a,    #F11 BMI
    ff2039,      #TF1 BMI (12.5)
    fg3139,      #TF2 BMI (13.5)
    fh3019,      #TF3 BMI (15.5)
    FJMR022a,    #TF4 BMI (17)
    
    #SDQ total difficulties
    kq348f,  #KQ (6.75) SDQ total difficulties (prorated)
    ku710b,  #KU (9.5)  SDQ total difficulties (prorated)
    kw6605b, #KW (11.5) SDQ total difficulties (prorated)
    ta7025f, #TA (13)   SDQ total difficulties (prorated)
    tc4025f, #TC (16.5) SDQ total difficulties (prorated)
    
    #CISR 17 data
    FJCI050,  # CISR Total score (out of 57)
    FJCI1000, # CISR Depression symptom score (out of 21)
    FJCI603,  # Mild Depression
    FJCI608,  # Moderate Depression
    FJCI609,  # Severe Depression
    
    FJCI104, # Loss of appetite
    FJCI105, # Increased appetite
    FJCI106, # Somatic Symptoms total score
    FJCI150, # Fatigue
    FJCI200, # Problems concentrating
    FJCI250, # Sleep symptom score
    FJCI253, # Sleeping too much
    FJCI300, # Irritable/short tempered
    FJCI350, # Feeling sad/miserable/depressed,
    FJCI352, # Anhedonia
    #  FJCI400, # worry score (removed due to low n)
    #  FJCI500, # physical anxiety symptoms (removed due to low n)
    FJCI501, # anxiety score
    #  FJCI550, # panic score (removed due to low n)
    
    #Earlier Blood Test Data - there is age 7 data  - n varies but for Hb its 5581
    #Lets take those variables with n > 5000 as they may provide useful information
    CHOL_F7,  # Cholesterol F7
    hb_F7  ,  # Haemoglobin F7
    IGE_F7 ,  # IgE F7
    HDL_F7 ,  # HDL F7
    LDL_F7 ,  # LDL F7
    TRIG_F7,  # Triglycerides F7
    TSH_F7 ,  # Thyroid stimulating Hormone F7
    
    #Other F7 variables have much lower n so are not going to be much use
    
    #Similarly, we have the F11 data, n ~ 4700 for hb
    hb_F11,  # Haemoglobin F11
    
    #All the other F11 data are regarding infections, which probably is not good
    #material for imputation
    
    #TF3 - age ~ 15 - N ~ 3000
    crp_TF3,     # CRP age 15
    hdl_TF3,     # HDL age 15
    ldl_TF3,     # LDL
    vldl_TF3,    # VLDL 
    trig_TF3,    # Triglycerides
    insulin_TF3, # Insulin
    GGT_TF3,     # Gamma GT
    glucose_TF3, # Glucose
    
    #TF4 - age ~ 17 - N ~ 3000
    CRP_TF4,     # CRP age 17
    HDL_TF4,     # HDL age 17
    LDL_TF4,     # LDL
    VLDL_TF4,    # VLDL
    TRIG_TF4,    # Triglycerides
    insulin_TF4, # Insulin
    glucose_TF4, # Glucose
    
    
    #Cognitive variables from earlier clinics

    #early WISC variables
    f8ws053, # F8 WISC vocabulary score
    f8ws055, # F8 WISC - Digit span scaled score: 
    f8ws110, # F8 WISC - Verbal IQ: F8
    f8ws111, # F8 WISC - Performance IQ: F8
    f8ws112, # F8 WISC - Total IQ: F8
    
    #WAIS
    fh6280   # TF4 WAIS total IQ
    
  )|>
  rename(id  = alnqlet,
         
         bmi7  = f7ms026a,    #F7 BMI
         bmi9  = f9ms026a,    #F9 BMI
         bmi10 = fdms026a,    #F10 BMI
         bmi11 = fems026a,    #F11 BMI
         bmi12 = ff2039,      #TF1 BMI (12.5)
         bmi13 = fg3139,      #TF2 BMI (13.5)
         bmi15 = fh3019,      #TF3 BMI (15.5)
         bmi17 = FJMR022a,    #TF4 BMI (17)
         
         sdq7  = kq348f,  #KQ (6.75) SDQ total difficulties (prorated)
         sdq9  = ku710b,  #KU (9.5)  SDQ total difficulties (prorated)
         sdq11 = kw6605b, #KW (11.5) SDQ total difficulties (prorated)
         sdq13 = ta7025f, #TA (13)   SDQ total difficulties (prorated)
         sdq16 = tc4025f, #TC (16.5) SDQ total difficulties (prorated)
         
         cisr_total_17 = FJCI050,  # CISR Total score (out of 57)
         cisr_deps_17  = FJCI1000, # CISR Depression symptom score (out of 21)
         cisr_mild_17  = FJCI603,  # Mild Depression
         cisr_mod_17   = FJCI608,  # Moderate Depression
         cisr_sev_17   = FJCI609,  # Severe Depression
         loss_app_17   = FJCI104, # Loss of appetite
         inc_app_17    = FJCI105, # Increased appetite
         somatic_17    = FJCI106, # Somatic Symptoms total score
         fatigue_17    = FJCI150, # Fatigue
         conc_17       = FJCI200, # Problems concentrating
         sleeps_17     = FJCI250, # Sleep symptom score
         sleep_inc_17  = FJCI253, # Sleeping too much
         irritable_17  = FJCI300, # Irritable/short tempered
         sad_17        = FJCI350, # Feeling sad/miserable/depressed,
         anhedonia_17  = FJCI352, # Anhedonia
         # worry_17      = FJCI400, # worry score = Removed due to excessive missingness
         # phys_anx_17   = FJCI500, # physical anxiety symptoms = Removed due to excessive missingness
         anx_17        = FJCI501, # anxiety score
         # panic_17      = FJCI550, # panic score = Removed due to excessive missingness
   
         
         #early WISC variables
         wisc_vocab_8     = f8ws053, # F8 WISC vocabulary score
         wisc_digit_8     = f8ws055, # F8 WISC - Digit span scaled score: 
         wisc_viq_8       = f8ws110, # F8 WISC - Verbal IQ: F8
         wisc_piq_8       = f8ws111, # F8 WISC - Performance IQ: F8
         wisc_tiq_9       = f8ws112, # F8 WISC - Total IQ: F8
         
         #WAIS
         wais_iq_17       = fh6280   # TF4 WAIS total IQ
         
         
  ) |> 
  janitor::clean_names() |>
  
  #Do some recoding
  mutate(across(-id,~case_when(. < 0 ~ NA_integer_,
                               TRUE ~ .))) |>
  
  #Transform variables we know are non-normally distributed
  mutate(across(.cols = c(tsh_f7,
                          crp_tf3, insulin_tf3, ggt_tf3,
                          crp_tf4, insulin_tf4), .fns = log)) 


# d_aux|>
# 
#   #Plot
#   select(-id) |>
#   summarytools::dfSummary() |>
#   summarytools::view()

# #Plot
# d_aux|>
#   select(-id) |>
#   visdat::vis_miss(warn_large_data = FALSE)

#Lots of missingness over the whole ALSPAC cohort

# Stick it all together -------


#Unite datasets
d_mi = 
  
  #Start with the f24 bloods
  d_24 |>
  select(-c(cisr_dep)) |>
  rename_with(~paste0(.x,'_f24'),.cols = -id) |>
  
  #Add the f9 bloods
  full_join(d_9 |> rename_with(~paste0(.x,'_f9'),.cols = -id),by = "id") |>
  
  #Add CISR diagnostic Data
  full_join(d_cisr, by = "id") |>
  
  #Add CISR Subscales
  full_join(d_cisr_scales, by = "id") |>
  
  #Add CISR symptom scores
  full_join(d_symptom_score, by = "id") |>
  
  #Add the SAPAS and CAPE data
  full_join(d_sapas, by = "id") |>
  full_join(d_cape, by = "id") |>
  
  #Add historic SMFQ Data
  full_join(d_smfq, by = "id") |>
  
  #Glue in covariates
  left_join(d_cv |> select(-c(bmi9,starts_with("smfq"))), by = "id") |>
  
  #Glue in auxillary variables
  left_join(d_aux, by = "id") 


#Visualise
# d_mi|>
#   #Remove the ID variable for plotting
#   select(-c(id)) |>
#   visdat::vis_miss(warn_large_data = FALSE) +
#   theme(axis.text.x.top = element_text(angle = 90,size = 6))

#The individuals who have bloods are definitely more present in the dataset than
#those who dropped out before then


#Next remove individuals who are NA for both the F24 and F9 bloods. Through MI we can
#potentially get up to ~41% more individuals included, which is a big deal,


#Lets start with a criteria that you have to have at least one blood test 
#(i.e. at F9 or F24)
d_mi = 
  d_mi |>
  filter(!is.na(il6_f24) | !is.na(il6_f9)) |>   #Gives 4164 participants
  #filter(!is.na(cisr_dep) | !is.na(il6_f24)) |> #Gives 3998 participants
  #drop_na(cisr_dep,il6_f24) |>                  #Gives 2954 participants
  
  #Lets also remove individuals with missing sex data as we probably ought
  #not to be imputing that
  filter(!is.na(sex)) #This doesn't actually change our n

#Inspect distributions
# d_mi|>
#   select(-id) |>
#   summarytools::dfSummary() |>
#   summarytools::view()

#Looks like smfq 18 and 21 have lots of missing data (>50%)
#sleep increased item @17 is >50%

#Interestingly, and helpfully, the group of participants we have with at least 
#one olink blood have only 45% missingness on the ace_total_classic variable
#which is at least a slight improvement (although still not ideal)

#Visualise again
# d_mi |>
#   visdat::vis_miss(warn_large_data = FALSE) +
#   theme(axis.text.x.top = element_text(angle = 90,size = 6))

#Overall from this sample, 28.8% of data is missing


#Lets do some sanity checking of variables for weird ALSPAC codes
# skimr::skim(d_mi)

#Looks like at this point we are free of any such problems, phew!


## - Wrangle ------

#Make sure our outcomes are coded as factors 
d_mi <- 
  d_mi |>
  relocate(id,cisr_dep,som:pan,atypical,psychological,somatic,anxiety,sapas_total,cape_total) |>
  mutate(across(.cols = cisr_dep:pan,.fns = factor)) 

#This is the point at which we set our dataset in terms of participants - we 
#remove some variables but from now onward we don't remove any rows apart from
#when dropping NAs for specific analysis


#Remove variables with > 50% missingness, with some exceptions
var_keep <- 
  d_mi |>
  summarise(across(everything(), ~.x |> is.na() |> sum()/nrow(d_mi) )) |>
  pivot_longer(everything(),names_to = "Variable",values_to = "prop_miss") |>
  arrange(-prop_miss) |>
  filter(prop_miss < 0.5 |
           str_detect(Variable,"smfq_c")) |>
  pull(Variable)
  
#We also need to make an exception for the age 28 AD variables as they only 
#apply to a subset of participants by design

#This removes the smfq18, smfq21, the increased sleep variable from
#the age 17 cisr,  two affect recognition task at 17 measures and ggt_tf3
  
#I also make an exception for the covid smfq data as its in the future compared
#to the f24 data and we aren't going to use it for anything except testing 
#predictive models

d_mi <- 
  d_mi |>
  select(all_of(var_keep)) |>
  relocate(id,cisr_dep,som:pan,atypical,psychological,somatic,anxiety,sapas_total,cape_total) 



#Much like the CNV project, we need to make sure all ordinal variables
#are coded as integers
is_int_var = 
  d_mi |> 
  select(-c(id:pan,sc,me,sex,eth,smk24)) |>
  pivot_longer(everything()) |> 
  group_by(name) |> 
  nest() |>
  mutate(all_int = map_lgl(data, ~ .x |>
                             drop_na() |> 
                             mutate(frac = map_lgl(value,~.x%%1 == 0)) |> 
                             summarise(is_int = all(frac)) |>
                             pull())) |>
  select(name,all_int) |> 
  filter(all_int == T) |> 
  pull(name)


#Now change to integers wherever appropriate
d_mi <- 
  d_mi |>
  mutate(across(all_of(is_int_var), as.integer)) 



# Save -----

#Lets save the split data object for use throughout the rest of the project
write_rds(d_mi,paste(data_dir,"alspac_data_final.rds",sep = '//'))

d <- d_mi

#Do some clearing up
rm(list = c("d_aux","d_cape","d_cisr","d_cisr_scales","d_corr_pairs","d_sapas","d_smfq","d_split","d_symptom_score","d_train","d_mi"))

