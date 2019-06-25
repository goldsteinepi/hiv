#################
# HIV simulation
# Citation: LeVasseur MT, Goldstein ND, Tabb LP, Olivieri-Mui BL, Welles SL. The Effect of PrEP on HIV Incidence Among Men Who Have Sex With Men in the Context of Condom Use, Treatment as Prevention, and Seroadaptive Practices. Manuscript in preparation.
# Note: Simulation dataset may be downloaded from: https://drive.google.com/open?id=1pCV3OMhfpIQcwDIuC1oIIIYrC1OuMeDz
# 8/22/16 -- Michael LeVasseur and Neal Goldstein
#################


### FUNCTIONS ###

prevention_paradigm = function(rand_vector, prevention_var, dataset, prep_level=NA)
{
  #TAPVL (viral load under TAP, 1=suppressed, 2=not suppressed), PROBCOND (condom probability), SEROTYPE (1=seropositioning/insertive, 2=seropositioning/receptive, 3=serosorting, 4=no seroadaptive behavior), ONPREP (0=not on prep, 1=on prep)
  if (prevention_var=="TAPVL") {
    starting_vals = ifelse(dataset$STATUS==3 & rand_vector<=0.42, 1, ifelse(dataset$STATUS==3 & rand_vector>0.42, 2, ifelse(dataset$STATUS==2 & dataset$HIV==1, 2, -1)))
  } else if (prevention_var=="PROBCOND") {
    starting_vals = ifelse(rand_vector<=0.25, 1, ifelse(rand_vector<=0.49, 0, ifelse(rand_vector<=0.70, 0.75, ifelse(rand_vector<=0.86, 0.50, 0.25))))
  } else if (prevention_var=="SEROTYPE") {
    starting_vals = ifelse(dataset$STATUS==3 & rand_vector<=0.133, 2, ifelse(dataset$STATUS==3 & rand_vector<=0.347, 3, ifelse(dataset$STATUS==3 & rand_vector>0.347, 4, ifelse(dataset$STATUS==2 & rand_vector<=0.053, 1, ifelse(dataset$STATUS==2 & rand_vector<=0.374, 3, ifelse(dataset$STATUS==2 & rand_vector>0.374, 4, ifelse(dataset$STATUS==1 & rand_vector<=0.095, 1, ifelse(dataset$STATUS==1 & rand_vector<=0.402, 3, ifelse(dataset$STATUS==1 & rand_vector>0.402, 4, NA)))))))))
  } else if (prevention_var=="ONPREP") {
    starting_vals = ifelse(dataset$STATUS==1 & rand_vector<=(prep_level/((1-0.19)-((1-0.19)*0.369))), 1, 0)
  }
  
  return(starting_vals)
}

### STEP0: GLOBAL INITIATION ###

#ensuring that all datasets are same
set.seed(777)

#ensure that each of 20 simulations is unique
seed_val_population = sample(1:1000,20,replace=F)

#initialize container for results of all simulations
prep_sims_20 = list(NA)

#20 simulations for precision estimates
for (sims in 1:20)
{
  cat("\n\n************** ","Simulation: ",sims," **************\n",sep="")
  
  ### STEP1: INITIALIZING BASE POPULATION ###
  
  #counterfactual within each simulation: ensures that characteristics of the population are the same within each simulation
  set.seed(seed_val_population[sims])
  
  #generate an initial population of 10,000 men, with same characteristics for all prevention scenarios
  baseline_pop = data.frame("ID"=1:10000, "HIV"=NA, "VL"=NA, "STATUS"=NA, "CIRC"=NA, "SEX_POSITION"=NA, "TEST_DAY"=NA, "MAX_PARTNERS"=NA, "ROUND"=0, "SEX_PARTNERS"=0, "ORIGINAL_STATUS"=NA, "PREP_PREVENT"=0, "TAP_PREVENT"=0, "SERO_PREVENT"=0, "COND_PREVENT"=0, "OVERALL_PREVENT"=0, "DISCORDANT"=0, "CAUSE_INFECT"=0, "FINISHED_SEX"=0, "DAYS_KNOWN_POSITIVE"=-1, "INCIDENT_DAYS"=-1, stringsAsFactors=F)
  
  #infect initial population with 19% prev(HIV); data from CDC
  baseline_pop$HIV = ifelse(runif(10000,0,1)<=0.19, 1, 0)
  
  #viral load category: 1=chronic, 2=acute (everyone at baseline is chronic)
  baseline_pop$VL = ifelse(baseline_pop$HIV==1, 1, -1)
  #baseline_pop$VL = ifelse(baseline_pop$HIV==1, 2, -1); baseline_pop$INCIDENT_DAYS = ifelse(baseline_pop$HIV==1, 0, -1) #for R0 calculation
  
  #knowledge of HIV status: 1=HIV negative and aware (63%), 2=unknown, 3=HIV positive and aware (66%); data from CDC
  baseline_pop$STATUS = ifelse(baseline_pop$HIV==0 & runif(10000,0,1)<=0.369, 2, ifelse(baseline_pop$HIV==0, 1, NA))
  baseline_pop$STATUS = ifelse(is.na(baseline_pop$STATUS) & baseline_pop$HIV==1 & runif(10000,0,1)<=0.440, 2, ifelse(baseline_pop$HIV==1, 3, baseline_pop$STATUS))
  
  #circumcision status; hospital record data from births in 1980
  baseline_pop$CIRC = ifelse(runif(10000,0,1)<=0.647, 1, 0)
  
  #sex position (scale of dominance), probability of being insertive
  baseline_pop$SEX_POSITION = runif(10000,0,1)
  
  #day of year will be tested for HIV, based on a future test probability of 1=tested 1/3 of year, 2=tested 2/3 of year, 3=tested 3/3 of year, 4=not tested; data from Kaiser Family Foundation
  testprob = runif(10000,0,1)
  baseline_pop$futuretest = ifelse(baseline_pop$STATUS==3, 4, ifelse(testprob<=0.19, 1, ifelse(testprob<=0.30, 2, ifelse(testprob<=0.70, 3, 4))))
  probdaytest = runif(10000,0,1)
  baseline_pop$TEST_DAY = ifelse(baseline_pop$futuretest==1, round(1+(110-1)*probdaytest), ifelse(baseline_pop$futuretest==2, round(111+(220-111)*probdaytest), ifelse(baseline_pop$futuretest==3, round(221+(330-221)*probdaytest), -1)))
  rm(testprob,probdaytest)
  baseline_pop$futuretest = NULL
  
  #partner distribution, based on a gamma distribution
  baseline_pop$MAX_PARTNERS = round(rgamma(10000, shape=0.5, scale=10)+1) #using Project Male-Call
  
  #starting HIV status for beginning of 1 year period
  baseline_pop$ORIGINAL_STATUS = baseline_pop$HIV
  
  
  ### STEP2: ASSIGNMENT OF PREVENTION PARADIGMS ###
  
  randtap = runif(10000,0,1)
  randcond = runif(10000,0,1)
  randsero = runif(10000,0,1)
  randprep = runif(10000,0,1)
  
  #create prevention paradigm datasets: abcd where a=prep, b=treatment as prevention, c=condom use, d=seroadaption
  paradigm1000 = baseline_pop
  paradigm1100 = baseline_pop
  paradigm1010 = baseline_pop
  paradigm1001 = baseline_pop
  paradigm1110 = baseline_pop
  paradigm1101 = baseline_pop
  paradigm1011 = baseline_pop
  paradigm1111 = baseline_pop
  
  #PrEP
  paradigm1000$PREP = 1
  paradigm1000$TAP = 0
  paradigm1000$COND = 0
  paradigm1000$SERO = 0
  paradigm1000$TAPVL = NA
  paradigm1000$PROBCOND = NA
  paradigm1000$SEROTYPE = NA
  
  #PrEP & TasP
  paradigm1100$PREP = 1
  paradigm1100$TAP = 1
  paradigm1100$COND = 0
  paradigm1100$SERO = 0
  paradigm1100$TAPVL = prevention_paradigm(randtap,"TAPVL",paradigm1100)
  paradigm1100$PROBCOND = NA
  paradigm1100$SEROTYPE = NA
  
  #PrEP & condoms
  paradigm1010$PREP = 1
  paradigm1010$TAP = 0
  paradigm1010$COND = 1
  paradigm1010$SERO = 0
  paradigm1010$TAPVL = NA
  paradigm1010$PROBCOND = prevention_paradigm(randcond,"PROBCOND",paradigm1010)
  paradigm1010$SEROTYPE = NA
  
  #PrEP & seroadaption
  paradigm1001$PREP = 1
  paradigm1001$TAP = 0
  paradigm1001$COND = 0
  paradigm1001$SERO = 1
  paradigm1001$TAPVL = NA
  paradigm1001$PROBCOND = NA
  paradigm1001$SEROTYPE = prevention_paradigm(randsero,"SEROTYPE",paradigm1001)
  
  #PrEP, TasP & condoms
  paradigm1110$PREP = 1
  paradigm1110$TAP = 1
  paradigm1110$COND = 1
  paradigm1110$SERO = 0
  paradigm1110$TAPVL = prevention_paradigm(randtap,"TAPVL",paradigm1110)
  paradigm1110$PROBCOND = prevention_paradigm(randcond,"PROBCOND",paradigm1110)
  paradigm1110$SEROTYPE = NA
  
  #PrEP, TasP & seroadaption
  paradigm1101$PREP = 1
  paradigm1101$TAP = 1
  paradigm1101$COND = 0
  paradigm1101$SERO = 1
  paradigm1101$TAPVL = prevention_paradigm(randtap,"TAPVL",paradigm1101)
  paradigm1101$PROBCOND = NA
  paradigm1101$SEROTYPE = prevention_paradigm(randsero,"SEROTYPE",paradigm1101)
  
  #PrEP, condoms & seroadaption
  paradigm1011$PREP = 1
  paradigm1011$TAP = 0
  paradigm1011$COND = 1
  paradigm1011$SERO = 1
  paradigm1011$TAPVL = NA
  paradigm1011$PROBCOND = prevention_paradigm(randcond,"PROBCOND",paradigm1011)
  paradigm1011$SEROTYPE = prevention_paradigm(randsero,"SEROTYPE",paradigm1011)
  
  #PrEP, TasP, condoms & seroadaption
  paradigm1111$PREP = 1
  paradigm1111$TAP = 1
  paradigm1111$COND = 1
  paradigm1111$SERO = 1
  paradigm1111$TAPVL = prevention_paradigm(randtap,"TAPVL",paradigm1111)
  paradigm1111$PROBCOND = prevention_paradigm(randcond,"PROBCOND",paradigm1111)
  paradigm1111$SEROTYPE = prevention_paradigm(randsero,"SEROTYPE",paradigm1111)
  
  #assign PrEP into lists
  prep_sims = list("Paradigm1000PrEP0"=cbind(paradigm1000,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1000,0)),
                   "Paradigm1000PrEP1"=cbind(paradigm1000,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1000,0.01)),
                   "Paradigm1000PrEP5"=cbind(paradigm1000,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1000,0.05)),
                   "Paradigm1000PrEP10"=cbind(paradigm1000,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1000,0.10)),
                   "Paradigm1000PrEP15"=cbind(paradigm1000,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1000,0.15)),
                   "Paradigm1000PrEP20"=cbind(paradigm1000,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1000,0.20)),
                   "Paradigm1000PrEP25"=cbind(paradigm1000,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1000,0.25)),
                   "Paradigm1100PrEP0"=cbind(paradigm1100,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1100,0)),
                   "Paradigm1100PrEP1"=cbind(paradigm1100,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1100,0.01)),
                   "Paradigm1100PrEP5"=cbind(paradigm1100,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1100,0.05)),
                   "Paradigm1100PrEP10"=cbind(paradigm1100,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1100,0.10)),
                   "Paradigm1100PrEP15"=cbind(paradigm1100,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1100,0.15)),
                   "Paradigm1100PrEP20"=cbind(paradigm1100,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1100,0.20)),
                   "Paradigm1100PrEP25"=cbind(paradigm1100,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1100,0.25)),
                   "Paradigm1010PrEP0"=cbind(paradigm1010,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1010,0)),
                   "Paradigm1010PrEP1"=cbind(paradigm1010,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1010,0.01)),
                   "Paradigm1010PrEP5"=cbind(paradigm1010,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1010,0.05)),
                   "Paradigm1010PrEP10"=cbind(paradigm1010,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1010,0.10)),
                   "Paradigm1010PrEP15"=cbind(paradigm1010,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1010,0.15)),
                   "Paradigm1010PrEP20"=cbind(paradigm1010,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1010,0.20)),
                   "Paradigm1010PrEP25"=cbind(paradigm1010,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1010,0.25)),
                   "Paradigm1001PrEP0"=cbind(paradigm1001,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1001,0)),
                   "Paradigm1001PrEP1"=cbind(paradigm1001,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1001,0.01)),
                   "Paradigm1001PrEP5"=cbind(paradigm1001,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1001,0.05)),
                   "Paradigm1001PrEP10"=cbind(paradigm1001,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1001,0.10)),
                   "Paradigm1001PrEP15"=cbind(paradigm1001,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1001,0.15)),
                   "Paradigm1001PrEP20"=cbind(paradigm1001,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1001,0.20)),
                   "Paradigm1001PrEP25"=cbind(paradigm1001,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1001,0.25)),
                   "Paradigm1110PrEP0"=cbind(paradigm1110,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1110,0)),
                   "Paradigm1110PrEP1"=cbind(paradigm1110,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1110,0.01)),
                   "Paradigm1110PrEP5"=cbind(paradigm1110,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1110,0.05)),
                   "Paradigm1110PrEP10"=cbind(paradigm1110,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1110,0.10)),
                   "Paradigm1110PrEP15"=cbind(paradigm1110,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1110,0.15)),
                   "Paradigm1110PrEP20"=cbind(paradigm1110,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1110,0.20)),
                   "Paradigm1110PrEP25"=cbind(paradigm1110,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1110,0.25)),
                   "Paradigm1101PrEP0"=cbind(paradigm1101,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1101,0)),
                   "Paradigm1101PrEP1"=cbind(paradigm1101,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1101,0.01)),
                   "Paradigm1101PrEP5"=cbind(paradigm1101,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1101,0.05)),
                   "Paradigm1101PrEP10"=cbind(paradigm1101,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1101,0.10)),
                   "Paradigm1101PrEP15"=cbind(paradigm1101,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1101,0.15)),
                   "Paradigm1101PrEP20"=cbind(paradigm1101,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1101,0.20)),
                   "Paradigm1101PrEP25"=cbind(paradigm1101,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1101,0.25)),
                   "Paradigm1011PrEP0"=cbind(paradigm1011,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1011,0)),
                   "Paradigm1011PrEP1"=cbind(paradigm1011,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1011,0.01)),
                   "Paradigm1011PrEP5"=cbind(paradigm1011,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1011,0.05)),
                   "Paradigm1011PrEP10"=cbind(paradigm1011,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1011,0.10)),
                   "Paradigm1011PrEP15"=cbind(paradigm1011,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1011,0.15)),
                   "Paradigm1011PrEP20"=cbind(paradigm1011,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1011,0.20)),
                   "Paradigm1011PrEP25"=cbind(paradigm1011,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1011,0.25)),
                   "Paradigm1111PrEP0"=cbind(paradigm1111,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1111,0)),
                   "Paradigm1111PrEP1"=cbind(paradigm1111,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1111,0.01)),
                   "Paradigm1111PrEP5"=cbind(paradigm1111,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1111,0.05)),
                   "Paradigm1111PrEP10"=cbind(paradigm1111,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1111,0.10)),
                   "Paradigm1111PrEP15"=cbind(paradigm1111,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1111,0.15)),
                   "Paradigm1111PrEP20"=cbind(paradigm1111,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1111,0.20)),
                   "Paradigm1111PrEP25"=cbind(paradigm1111,"ONPREP"=prevention_paradigm(randprep,"ONPREP",paradigm1111,0.25)))
  
  #cleanup
  rm(randtap,randcond,randsero,randprep,baseline_pop,paradigm1000,paradigm1100,paradigm1010,paradigm1001,paradigm1110,paradigm1101,paradigm1011,paradigm1111)
  gc()
  
  ### STEP3: SIMULATE ###
  
  #ensure that each day is unique
  seed_val_day = sample(1:1000,365,replace=F)
  
  for (day in 1:365)
  {
    cat("\n\n************** ","Day: ",day," **************\n",sep="")
    
    #go through each scenario
    for (current_data in 1:length(prep_sims))
    {
      #cat("\n\n************** ","Scenario: ",current_data," **************\n",sep="")
      
      ## SET DAILY VARIABLES ##
      
      #ensure each scenario is identical for a given day
      set.seed(seed_val_day[day])
      
      #check if reached max number of partners
      prep_sims[[current_data]]$FINISHED_SEX = ifelse(prep_sims[[current_data]]$SEX_PARTNERS>=prep_sims[[current_data]]$MAX_PARTNERS, 1, 0)
      
      #HIV testing today?
      prep_sims[[current_data]]$STATUS = ifelse(day==prep_sims[[current_data]]$TEST_DAY & prep_sims[[current_data]]$HIV==1, 3, prep_sims[[current_data]]$STATUS)
      
      #check if testing identify as positive, used for TAP, after 30 days you may be virally suppressed
      prep_sims[[current_data]]$DAYS_KNOWN_POSITIVE = ifelse(day==prep_sims[[current_data]]$TEST_DAY & prep_sims[[current_data]]$HIV==1, 1, prep_sims[[current_data]]$DAYS_KNOWN_POSITIVE)
      
      #increment days positive, if positive
      prep_sims[[current_data]]$DAYS_KNOWN_POSITIVE = ifelse(prep_sims[[current_data]]$DAYS_KNOWN_POSITIVE>=0, prep_sims[[current_data]]$DAYS_KNOWN_POSITIVE+1, prep_sims[[current_data]]$DAYS_KNOWN_POSITIVE)
      
      #in a TAP scenario, if put on treatment, after 30 days viral load goes to undetectable
      prep_sims[[current_data]]$TAPVL = ifelse(prep_sims[[current_data]]$TAP==1 & prep_sims[[current_data]]$DAYS_KNOWN_POSITIVE==30 & prep_sims[[current_data]]$HIV==1 & runif(nrow(prep_sims[[current_data]]),0,1)<=0.42, 1, prep_sims[[current_data]]$TAPVL)
      
      #increment incident days
      prep_sims[[current_data]]$INCIDENT_DAYS = ifelse(prep_sims[[current_data]]$INCIDENT_DAYS>=0, prep_sims[[current_data]]$INCIDENT_DAYS+1, prep_sims[[current_data]]$INCIDENT_DAYS)
      
      #resolve acute phase of HIV infection
      prep_sims[[current_data]]$VL = ifelse(prep_sims[[current_data]]$VL==2 & prep_sims[[current_data]]$INCIDENT_DAYS==30, 1, prep_sims[[current_data]]$VL)
      
      #if HIV positive and aware and set as receptive
      prep_sims[[current_data]]$SEROTYPE = ifelse(prep_sims[[current_data]]$STATUS==3 & prep_sims[[current_data]]$SEROTYPE==1, 2, prep_sims[[current_data]]$SEROTYPE) 
      
      ## ENGAGE IN SEX ##
      
      #determine who has sex today
      sex_list = prep_sims[[current_data]]$ID[prep_sims[[current_data]]$FINISHED_SEX==0 & (((prep_sims[[current_data]]$MAX_PARTNERS-prep_sims[[current_data]]$SEX_PARTNERS)/(366-day))>=runif(nrow(prep_sims[[current_data]]),0,1))]
      
      #now allocate those individuals randomly to two even groups to see who has sex; if sex_list is an odd length, one individual will not have sex
      sex_list1 = sample(sex_list, floor(length(sex_list)/2), replace=F)
      sex_list2 = sample(sex_list[!sex_list %in% sex_list1], length(sex_list1), replace=F)

      #determining concordance and discordance by HIV status
      sex_exposed_hiv = data.frame("Ego"=sex_list1[which((prep_sims[[current_data]]$HIV[sex_list1] + prep_sims[[current_data]]$HIV[sex_list2])==1)], stringsAsFactors=F)
      sex_exposed_hiv$Partner = sex_list2[which((prep_sims[[current_data]]$HIV[sex_list1] + prep_sims[[current_data]]$HIV[sex_list2])==1)]
      sex_unexposed_hiv = c(sex_list1[which((prep_sims[[current_data]]$HIV[sex_list1] + prep_sims[[current_data]]$HIV[sex_list2])!=1)], sex_list2[which((prep_sims[[current_data]]$HIV[sex_list1] + prep_sims[[current_data]]$HIV[sex_list2])!=1)])
      rm(sex_list,sex_list1,sex_list2)
  
      if (nrow(sex_exposed_hiv)>0)
      {
        ## HIV DISCORDANT ##
        
        ## DETERMINE INFECTION STATUS ##
        
        #add HIV info, 0=negative, 1=positive
        sex_exposed_hiv$HIV = prep_sims[[current_data]]$HIV[sex_exposed_hiv$Ego]
        sex_exposed_hiv$partHIV = prep_sims[[current_data]]$HIV[sex_exposed_hiv$Partner]
        
        #check for sex position, 1=ego bottom/partner top, 2=ego top/partner bottom, 3=flip fuck 
        sex_exposed_hiv$simAI = ifelse(abs(prep_sims[[current_data]]$SEX_POSITION[sex_exposed_hiv$Ego]-prep_sims[[current_data]]$SEX_POSITION[sex_exposed_hiv$Partner])<=0.15, 3, ifelse(prep_sims[[current_data]]$SEX_POSITION[sex_exposed_hiv$Ego]<prep_sims[[current_data]]$SEX_POSITION[sex_exposed_hiv$Partner], 1, 2))
    
        #check for who is exposed, 1=ego, 2=partner
        sex_exposed_hiv$exposure = ifelse(prep_sims[[current_data]]$HIV[sex_exposed_hiv$Ego]<prep_sims[[current_data]]$HIV[sex_exposed_hiv$Partner], 1, 2)
    
        #indicators for various metrics
        sex_exposed_hiv$condct = 0       #infection blocked from condom use, ego
        sex_exposed_hiv$prepct = 0       #infection blocked from PrEP, ego
        sex_exposed_hiv$tapct = 0        #infection blocked from TAP, ego
        sex_exposed_hiv$seroct = 0       #infection blocked from seroadapation, ego
        sex_exposed_hiv$prevct = 0       #infection blocked from any prevention strategy, ego
        sex_exposed_hiv$infectct = 0     #ego caused a new HIV infection to partner
        sex_exposed_hiv$newHIV = 0       #new HIV infection, ego
        sex_exposed_hiv$newVL = -1       #set VL to acute if new infection, ego
        sex_exposed_hiv$partcondct = 0   #infection blocked from condom use, partner
        sex_exposed_hiv$partprepct = 0   #infection blocked from PrEP, partner
        sex_exposed_hiv$parttapct = 0    #infection blocked from TAP,partner
        sex_exposed_hiv$partseroct = 0   #infection blocked from seroadapation,partner
        sex_exposed_hiv$partprevct = 0   #infection blocked from any prevention strategy, partner
        sex_exposed_hiv$partinfectct = 0 #partner caused a new HIV infection to ego
        sex_exposed_hiv$partnewHIV = 0   #new HIV infection, partner
        sex_exposed_hiv$partnewVL = -1   #set VL to acute if new infection, partner
        
        #calculate probabilities of infection
        probinfect1 = runif(nrow(sex_exposed_hiv),0,1)
        probinfect2 = runif(nrow(sex_exposed_hiv),0,1)

        #determine infection under no prevention, ego then partner, chronic then acute, 
        #sex_exposed_hiv$infection = 0         #new HIV infection, ego
        #sex_exposed_hiv$partinfection = 0     #new HIV infection, partner
        sex_exposed_hiv$infection = ifelse(sex_exposed_hiv$exposure==1 & sex_exposed_hiv$simAI==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==1 & probinfect1 <= 0.0134, 1, ifelse(sex_exposed_hiv$exposure==1 & sex_exposed_hiv$simAI==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==1 & probinfect1 <= 0.0010, 1, ifelse(sex_exposed_hiv$exposure==1 & sex_exposed_hiv$simAI==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==0 & probinfect1 <= 0.0059, 1, ifelse(sex_exposed_hiv$exposure==1 & sex_exposed_hiv$simAI==3 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==1 & probinfect1 <= 0.0134, 1, ifelse(sex_exposed_hiv$exposure==1 & sex_exposed_hiv$simAI==3 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==1 & probinfect2 <= 0.0010, 1, ifelse(sex_exposed_hiv$exposure==1 & sex_exposed_hiv$simAI==3 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==0 & probinfect2 <= 0.0059, 1, 0))))))
        sex_exposed_hiv$infection = ifelse(sex_exposed_hiv$exposure==1 & sex_exposed_hiv$simAI==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==2 & probinfect1 <= 0.1284, 1, ifelse(sex_exposed_hiv$exposure==1 & sex_exposed_hiv$simAI==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==1 & probinfect1 <= 0.0101, 1, ifelse(sex_exposed_hiv$exposure==1 & sex_exposed_hiv$simAI==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==0 & probinfect1 <= 0.0569, 1, ifelse(sex_exposed_hiv$exposure==1 & sex_exposed_hiv$simAI==3 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==2 & probinfect1 <= 0.1284, 1, ifelse(sex_exposed_hiv$exposure==1 & sex_exposed_hiv$simAI==3 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==1 & probinfect2 <= 0.0101, 1, ifelse(sex_exposed_hiv$exposure==1 & sex_exposed_hiv$simAI==3 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==0 & probinfect2 <= 0.0569, 1, sex_exposed_hiv$infection))))))
        sex_exposed_hiv$partinfection = ifelse(sex_exposed_hiv$exposure==2 & sex_exposed_hiv$simAI==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==1 & probinfect1 <= 0.0134, 1, ifelse(sex_exposed_hiv$exposure==2 & sex_exposed_hiv$simAI==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==1 & probinfect1 <= 0.0010, 1, ifelse(sex_exposed_hiv$exposure==2 & sex_exposed_hiv$simAI==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==0 & probinfect1 <= 0.0059, 1, ifelse(sex_exposed_hiv$exposure==2 & sex_exposed_hiv$simAI==3 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==1 & probinfect1 <= 0.0134, 1, ifelse(sex_exposed_hiv$exposure==2 & sex_exposed_hiv$simAI==3 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==1 & probinfect2 <= 0.0010, 1, ifelse(sex_exposed_hiv$exposure==2 & sex_exposed_hiv$simAI==3 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==0 & probinfect2 <= 0.0059, 1, 0))))))
        sex_exposed_hiv$partinfection = ifelse(sex_exposed_hiv$exposure==2 & sex_exposed_hiv$simAI==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==2 & probinfect1 <= 0.1284, 1, ifelse(sex_exposed_hiv$exposure==2 & sex_exposed_hiv$simAI==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==1 & probinfect1 <= 0.0101, 1, ifelse(sex_exposed_hiv$exposure==2 & sex_exposed_hiv$simAI==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==0 & probinfect1 <= 0.0569, 1, ifelse(sex_exposed_hiv$exposure==2 & sex_exposed_hiv$simAI==3 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==2 & probinfect1 <= 0.1284, 1, ifelse(sex_exposed_hiv$exposure==2 & sex_exposed_hiv$simAI==3 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==1 & probinfect2 <= 0.0101, 1, ifelse(sex_exposed_hiv$exposure==2 & sex_exposed_hiv$simAI==3 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==0 & probinfect2 <= 0.0569, 1, sex_exposed_hiv$partinfection))))))
        
        #determine infection under seroadaption scenarios
        if (unique(prep_sims[[current_data]]$SERO)==1)
        {
          #check if avoiding sex
          avoid = ifelse(sex_exposed_hiv$infection==1 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==1, 1, ifelse(sex_exposed_hiv$infection==1 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==3 & prep_sims[[current_data]]$STATUS[sex_exposed_hiv$Partner]==3, 1, ifelse(sex_exposed_hiv$infection==1 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==3 & prep_sims[[current_data]]$STATUS[sex_exposed_hiv$Partner]==3, 1, 0)))
          partavoid = ifelse(sex_exposed_hiv$partinfection==1 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==1, 1, ifelse(sex_exposed_hiv$partinfection==1 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==3 & prep_sims[[current_data]]$STATUS[sex_exposed_hiv$Ego]==3, 1, ifelse(sex_exposed_hiv$partinfection==1 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==3 & prep_sims[[current_data]]$STATUS[sex_exposed_hiv$Ego]==3, 1, 0)))
    
          #seroadaption strategies
          seroblock = ifelse(sex_exposed_hiv$infection==1 & avoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==1 & probinfect1 > 0.0010, 1, ifelse(sex_exposed_hiv$infection==1 & avoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==2 & probinfect1 > 0.0059, 1, ifelse(sex_exposed_hiv$infection==1 & avoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==1 & probinfect1 > 0.0101, 1, ifelse(sex_exposed_hiv$infection==1 & avoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==2 & probinfect1 > 0.0569, 1, ifelse(sex_exposed_hiv$infection==1 & avoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==1 & probinfect1 > 0.0010, 1, ifelse(sex_exposed_hiv$infection==1 & avoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==2 & probinfect1 > 0.0059, 1, ifelse(sex_exposed_hiv$infection==1 & avoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==1 & probinfect1 > 0.0101, 1, ifelse(sex_exposed_hiv$infection==1 & avoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Ego]==2 & probinfect1 > 0.0569, 1, ifelse(sex_exposed_hiv$infection==1 & avoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==1 & probinfect1 > 0.0134, 1, ifelse(sex_exposed_hiv$infection==1 & avoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner]==2 & probinfect1 > 0.1284, 1, 0))))))))))
          partseroblock = ifelse(sex_exposed_hiv$partinfection==1 & partavoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==1 & probinfect1 > 0.0010, 1, ifelse(sex_exposed_hiv$partinfection==1 & partavoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==2 & probinfect1 > 0.0059, 1, ifelse(sex_exposed_hiv$partinfection==1 & partavoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==1 & probinfect1 > 0.0101, 1, ifelse(sex_exposed_hiv$partinfection==1 & partavoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Partner]==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==2 & probinfect1 > 0.0569, 1, ifelse(sex_exposed_hiv$partinfection==1 & partavoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==1 & probinfect1 > 0.0010, 1, ifelse(sex_exposed_hiv$partinfection==1 & partavoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==2 & probinfect1 > 0.0059, 1, ifelse(sex_exposed_hiv$partinfection==1 & partavoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==1 & probinfect1 > 0.0101, 1, ifelse(sex_exposed_hiv$partinfection==1 & partavoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==2 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==2 & prep_sims[[current_data]]$CIRC[sex_exposed_hiv$Partner]==2 & probinfect1 > 0.0569, 1, ifelse(sex_exposed_hiv$partinfection==1 & partavoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==1 & probinfect1 > 0.0134, 1, ifelse(sex_exposed_hiv$partinfection==1 & partavoid==0 & prep_sims[[current_data]]$SEROTYPE[sex_exposed_hiv$Ego]==1 & prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego]==2 & probinfect1 > 0.1284, 1, 0))))))))))

          #check if infection blocked
          sex_exposed_hiv$seroct = ifelse(sex_exposed_hiv$infection==1 & avoid==1, 1, ifelse(sex_exposed_hiv$infection==1 & avoid==0 & seroblock==1, 1, 0))
          sex_exposed_hiv$partseroct = ifelse(sex_exposed_hiv$partinfection==1 & partavoid==1, 1, ifelse(sex_exposed_hiv$partinfection==1 & partavoid==0 & partseroblock==1, 1, 0))
          
          rm(avoid,partavoid,seroblock,partseroblock)
        }
        
        #determine infection under condom scenarios
        if (unique(prep_sims[[current_data]]$COND)==1)
        {
          #check if a condom was worn by either partner
          randcond = runif(nrow(sex_exposed_hiv),0,1)
          condom_worn = (randcond<=prep_sims[[current_data]]$PROBCOND[sex_exposed_hiv$Ego]) | (randcond<=prep_sims[[current_data]]$PROBCOND[sex_exposed_hiv$Partner])

          #check if the condom failed
          randfail = runif(nrow(sex_exposed_hiv),0,1)
          condom_success = condom_worn & (randfail<=0.705)

          #check if infection blocked
          sex_exposed_hiv$condct = ifelse(sex_exposed_hiv$infection==1 & condom_success, 1, 0)
          sex_exposed_hiv$partcondct = ifelse(sex_exposed_hiv$partinfection==1 & condom_success, 1, 0)
          
          rm(randcond,randfail,condom_worn,condom_success)
        }
        
        #determine infection under TAP scenarios
        if (unique(prep_sims[[current_data]]$TAP)==1)
        {
          #check if infection blocked
          sex_exposed_hiv$tapct = ifelse(sex_exposed_hiv$infection==1 & prep_sims[[current_data]]$TAPVL[sex_exposed_hiv$Partner]==1, 1, 0)
          sex_exposed_hiv$parttapct = ifelse(sex_exposed_hiv$partinfection==1 & prep_sims[[current_data]]$TAPVL[sex_exposed_hiv$Ego]==1, 1, 0)
        }
        
        #determine infection under PrEP scenarios
        if (unique(prep_sims[[current_data]]$PREP)==1)
        {
          #check if infection blocked
          sex_exposed_hiv$prepct = ifelse(sex_exposed_hiv$infection==1 & prep_sims[[current_data]]$ONPREP[sex_exposed_hiv$Ego]==1, 1, sex_exposed_hiv$prepct)
          sex_exposed_hiv$partprepct = ifelse(sex_exposed_hiv$partinfection==1 & prep_sims[[current_data]]$ONPREP[sex_exposed_hiv$Partner]==1, 1, sex_exposed_hiv$partprepct)
        }
        
        #resolve infection statistics
        anyprev = sex_exposed_hiv$condct + sex_exposed_hiv$prepct + sex_exposed_hiv$tapct + sex_exposed_hiv$seroct
        anypartprev = sex_exposed_hiv$partcondct + sex_exposed_hiv$partprepct + sex_exposed_hiv$parttapct + sex_exposed_hiv$partseroct
        sex_exposed_hiv$prevct = ifelse(anyprev>0 & sex_exposed_hiv$infection==1, 1, sex_exposed_hiv$prevct)
        sex_exposed_hiv$partprevct = ifelse(anypartprev>0 & sex_exposed_hiv$partinfection==1, 1, sex_exposed_hiv$partprevct)
        sex_exposed_hiv$newHIV = ifelse(anyprev==0 & sex_exposed_hiv$infection==1, 1, sex_exposed_hiv$newHIV)
        sex_exposed_hiv$newVL = ifelse(anyprev==0 & sex_exposed_hiv$infection==1, 2, sex_exposed_hiv$newVL)
        sex_exposed_hiv$partinfectct = ifelse(anyprev==0 & sex_exposed_hiv$infection==1, 1, sex_exposed_hiv$partinfectct)
        sex_exposed_hiv$partnewHIV = ifelse(anypartprev==0 & sex_exposed_hiv$partinfection==1, 1, sex_exposed_hiv$partnewHIV)
        sex_exposed_hiv$partnewVL = ifelse(anypartprev==0 & sex_exposed_hiv$partinfection==1, 2, sex_exposed_hiv$partnewVL)
        sex_exposed_hiv$infectct = ifelse(anypartprev==0 & sex_exposed_hiv$partinfection==1, 1, sex_exposed_hiv$infectct)
        
        ## UPDATE STATS IN MAIN DATASET ##
    
        prep_sims[[current_data]]$PREP_PREVENT[sex_exposed_hiv$Ego] = prep_sims[[current_data]]$PREP_PREVENT[sex_exposed_hiv$Ego] + sex_exposed_hiv$prepct
        prep_sims[[current_data]]$PREP_PREVENT[sex_exposed_hiv$Partner] = prep_sims[[current_data]]$PREP_PREVENT[sex_exposed_hiv$Partner] + sex_exposed_hiv$partprepct
        prep_sims[[current_data]]$SERO_PREVENT[sex_exposed_hiv$Ego] = prep_sims[[current_data]]$SERO_PREVENT[sex_exposed_hiv$Ego] + sex_exposed_hiv$seroct
        prep_sims[[current_data]]$SERO_PREVENT[sex_exposed_hiv$Partner] = prep_sims[[current_data]]$SERO_PREVENT[sex_exposed_hiv$Partner] + sex_exposed_hiv$partseroct
        prep_sims[[current_data]]$COND_PREVENT[sex_exposed_hiv$Ego] = prep_sims[[current_data]]$COND_PREVENT[sex_exposed_hiv$Ego] + sex_exposed_hiv$condct
        prep_sims[[current_data]]$COND_PREVENT[sex_exposed_hiv$Partner] = prep_sims[[current_data]]$COND_PREVENT[sex_exposed_hiv$Partner] + sex_exposed_hiv$partcondct
        prep_sims[[current_data]]$TAP_PREVENT[sex_exposed_hiv$Ego] = prep_sims[[current_data]]$TAP_PREVENT[sex_exposed_hiv$Ego] + sex_exposed_hiv$tapct
        prep_sims[[current_data]]$TAP_PREVENT[sex_exposed_hiv$Partner] = prep_sims[[current_data]]$TAP_PREVENT[sex_exposed_hiv$Partner] + sex_exposed_hiv$parttapct
        prep_sims[[current_data]]$OVERALL_PREVENT[sex_exposed_hiv$Ego] = prep_sims[[current_data]]$OVERALL_PREVENT[sex_exposed_hiv$Ego] + sex_exposed_hiv$prevct
        prep_sims[[current_data]]$OVERALL_PREVENT[sex_exposed_hiv$Partner] = prep_sims[[current_data]]$OVERALL_PREVENT[sex_exposed_hiv$Partner] + sex_exposed_hiv$partprevct
        prep_sims[[current_data]]$SEX_PARTNERS[sex_exposed_hiv$Ego] = prep_sims[[current_data]]$SEX_PARTNERS[sex_exposed_hiv$Ego] + 1
        prep_sims[[current_data]]$SEX_PARTNERS[sex_exposed_hiv$Partner] = prep_sims[[current_data]]$SEX_PARTNERS[sex_exposed_hiv$Partner] + 1
        prep_sims[[current_data]]$SEX_PARTNERS[sex_unexposed_hiv] = prep_sims[[current_data]]$SEX_PARTNERS[sex_unexposed_hiv] + 1
        prep_sims[[current_data]]$CAUSE_INFECT[sex_exposed_hiv$Ego] = ifelse(sex_exposed_hiv$infectct==1, prep_sims[[current_data]]$CAUSE_INFECT[sex_exposed_hiv$Ego] + sex_exposed_hiv$infectct, prep_sims[[current_data]]$CAUSE_INFECT[sex_exposed_hiv$Ego])
        prep_sims[[current_data]]$CAUSE_INFECT[sex_exposed_hiv$Partner] = ifelse(sex_exposed_hiv$partinfectct==1, prep_sims[[current_data]]$CAUSE_INFECT[sex_exposed_hiv$Partner] + sex_exposed_hiv$partinfectct, prep_sims[[current_data]]$CAUSE_INFECT[sex_exposed_hiv$Partner])
        prep_sims[[current_data]]$HIV[sex_exposed_hiv$Ego] = ifelse(sex_exposed_hiv$newHIV==1, 1, prep_sims[[current_data]]$HIV[sex_exposed_hiv$Ego])
        prep_sims[[current_data]]$HIV[sex_exposed_hiv$Partner] = ifelse(sex_exposed_hiv$partnewHIV==1, 1, prep_sims[[current_data]]$HIV[sex_exposed_hiv$Partner])
        prep_sims[[current_data]]$INCIDENT_DAYS[sex_exposed_hiv$Ego] = ifelse(sex_exposed_hiv$newHIV==1, 0, prep_sims[[current_data]]$INCIDENT_DAYS[sex_exposed_hiv$Ego])
        prep_sims[[current_data]]$INCIDENT_DAYS[sex_exposed_hiv$Partner] = ifelse(sex_exposed_hiv$partnewHIV==1, 0, prep_sims[[current_data]]$INCIDENT_DAYS[sex_exposed_hiv$Partner])
        prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego] = ifelse(sex_exposed_hiv$newVL==2, 2, prep_sims[[current_data]]$VL[sex_exposed_hiv$Ego])
        prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner] = ifelse(sex_exposed_hiv$partnewVL==2, 2, prep_sims[[current_data]]$VL[sex_exposed_hiv$Partner])
        prep_sims[[current_data]]$DISCORDANT[sex_exposed_hiv$Ego] = prep_sims[[current_data]]$DISCORDANT[sex_exposed_hiv$Ego] + 1
        prep_sims[[current_data]]$DISCORDANT[sex_exposed_hiv$Partner] = prep_sims[[current_data]]$DISCORDANT[sex_exposed_hiv$Partner] + 1
        
        ## CLEAN UP ##
        
        rm(sex_exposed_hiv,anyprev,anypartprev,probinfect1,probinfect2,sex_unexposed_hiv)
        gc()
        
      } else {
        
        ## HIV CONCORDANT ##
        
        ## UPDATE STATS IN MAIN DATASET ##
        
        prep_sims[[current_data]]$SEX_PARTNERS[sex_unexposed_hiv] = prep_sims[[current_data]]$SEX_PARTNERS[sex_unexposed_hiv] + 1
        
        ## CLEAN UP ##
        
        rm(sex_exposed_hiv,sex_unexposed_hiv)
        gc()
        
      }
      
    }

    ## TRACK DAILY TOTALS FOR R0 CALCULATION ##
    
    #set everyone to primary infection VL=2, eliminate chronic transmission, and create a vector r0 = rep(NA,365)
    #track new primary infection per day
    #r0[day] = sum(prep_sims[["Paradigm000PrEP0"]]$HIV) - sum(prep_sims[["Paradigm000PrEP0"]]$ORIGINAL_STATUS) - sum(prep_sims[["Paradigm000PrEP0"]]$INCIDENT_DAYS[prep_sims[["Paradigm000PrEP0"]]$ORIGINAL_STATUS==0]>30)
    #estimate.R(r0,methods=c("AR"), pop.size=10000, S0=0.81)
    
  }
  rm(day,current_data)

  #add this simulation to the list
  names(prep_sims) = paste("simulation",sims,names(prep_sims),sep="_")
  prep_sims_20 = c(prep_sims_20,prep_sims)
  
}
rm(sims,prep_sims,seed_val_day,seed_val_population,prevention_paradigm)
prep_sims_20[[1]] = NULL
gc()


### STEP4: TALLY RESULTS ###

#scenario labels
paradigm = c("None","TAP","Condom","Seroadaption","TAP and Condom", "TAP and Seroadaption", "Seroadaption and Condom", "All")
prep_level = c(0, 1, 5, 10, 15, 20, 25)
scenarios = length(paradigm)*length(prep_level)

#create mean number of new infections data frame
results_infections = data.frame("Paradigm"=NA,"PrEP"=NA,"Mean"=NA,"SE"=NA,stringsAsFactors=F)
for (i in 0:(length(paradigm)-1))
{
  for (j in 1:length(prep_level))
  {
    current = i*length(prep_level)+j
    statistic = c((sum(prep_sims_20[[current]]$HIV) - sum(prep_sims_20[[current]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*1)]]$HIV) - sum(prep_sims_20[[current+(scenarios*1)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*2)]]$HIV) - sum(prep_sims_20[[current+(scenarios*2)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*3)]]$HIV) - sum(prep_sims_20[[current+(scenarios*3)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*4)]]$HIV) - sum(prep_sims_20[[current+(scenarios*4)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*5)]]$HIV) - sum(prep_sims_20[[current+(scenarios*5)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*6)]]$HIV) - sum(prep_sims_20[[current+(scenarios*6)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*7)]]$HIV) - sum(prep_sims_20[[current+(scenarios*7)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*8)]]$HIV) - sum(prep_sims_20[[current+(scenarios*8)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*9)]]$HIV) - sum(prep_sims_20[[current+(scenarios*9)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*10)]]$HIV) - sum(prep_sims_20[[current+(scenarios*10)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*11)]]$HIV) - sum(prep_sims_20[[current+(scenarios*11)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*12)]]$HIV) - sum(prep_sims_20[[current+(scenarios*12)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*13)]]$HIV) - sum(prep_sims_20[[current+(scenarios*13)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*14)]]$HIV) - sum(prep_sims_20[[current+(scenarios*14)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*15)]]$HIV) - sum(prep_sims_20[[current+(scenarios*15)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*16)]]$HIV) - sum(prep_sims_20[[current+(scenarios*16)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*17)]]$HIV) - sum(prep_sims_20[[current+(scenarios*17)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*18)]]$HIV) - sum(prep_sims_20[[current+(scenarios*18)]]$ORIGINAL_STATUS)),
                  (sum(prep_sims_20[[current+(scenarios*19)]]$HIV) - sum(prep_sims_20[[current+(scenarios*19)]]$ORIGINAL_STATUS)))
    results_infections = rbind(results_infections, data.frame("Paradigm"=paradigm[i+1],"PrEP"=prep_level[j],"Mean"=mean(statistic),"SE"=(sd(statistic)/sqrt(length(statistic))),stringsAsFactors=F))
  }
}
rm(i,j,current,statistic)
results_infections = results_infections[-1, ]

#create mean percent of infections prevented data frame
results_prevented = data.frame("Paradigm"=NA,"PrEP"=NA,"Mean"=NA,"SE"=NA,stringsAsFactors=F)
for (i in 0:(length(paradigm)-1))
{
  for (j in 1:length(prep_level))
  {
    current = i*length(prep_level)+j
    statistic = c((sum(prep_sims_20[[current]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current]]$HIV) - sum(prep_sims_20[[current]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*1)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*1)]]$HIV) - sum(prep_sims_20[[current+(scenarios*1)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*1)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*2)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*2)]]$HIV) - sum(prep_sims_20[[current+(scenarios*2)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*2)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*3)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*3)]]$HIV) - sum(prep_sims_20[[current+(scenarios*3)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*3)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*4)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*4)]]$HIV) - sum(prep_sims_20[[current+(scenarios*4)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*4)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*5)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*5)]]$HIV) - sum(prep_sims_20[[current+(scenarios*5)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*5)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*6)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*6)]]$HIV) - sum(prep_sims_20[[current+(scenarios*6)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*6)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*7)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*7)]]$HIV) - sum(prep_sims_20[[current+(scenarios*7)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*7)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*8)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*8)]]$HIV) - sum(prep_sims_20[[current+(scenarios*8)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*8)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*9)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*9)]]$HIV) - sum(prep_sims_20[[current+(scenarios*9)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*9)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*10)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*10)]]$HIV) - sum(prep_sims_20[[current+(scenarios*10)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*10)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*11)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*11)]]$HIV) - sum(prep_sims_20[[current+(scenarios*11)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*11)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*12)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*12)]]$HIV) - sum(prep_sims_20[[current+(scenarios*12)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*12)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*13)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*13)]]$HIV) - sum(prep_sims_20[[current+(scenarios*13)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*13)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*14)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*14)]]$HIV) - sum(prep_sims_20[[current+(scenarios*14)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*14)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*15)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*15)]]$HIV) - sum(prep_sims_20[[current+(scenarios*15)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*15)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*16)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*16)]]$HIV) - sum(prep_sims_20[[current+(scenarios*16)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*16)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*17)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*17)]]$HIV) - sum(prep_sims_20[[current+(scenarios*17)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*17)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*18)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*18)]]$HIV) - sum(prep_sims_20[[current+(scenarios*18)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*18)]]$OVERALL_PREVENT))),
                  (sum(prep_sims_20[[current+(scenarios*19)]]$OVERALL_PREVENT) / (sum(prep_sims_20[[current+(scenarios*19)]]$HIV) - sum(prep_sims_20[[current+(scenarios*19)]]$ORIGINAL_STATUS) + sum(prep_sims_20[[current+(scenarios*19)]]$OVERALL_PREVENT))))
    results_prevented = rbind(results_prevented, data.frame("Paradigm"=paradigm[i+1],"PrEP"=prep_level[j],"Mean"=mean(statistic)*100,"SE"=(sd(statistic)/sqrt(length(statistic)))*100,stringsAsFactors=F))
  }
}
rm(i,j,current,statistic,scenarios)
results_prevented = results_prevented[-1, ]


### SAVE RESULTS ###

save.image("simulation results.RData")


### PLOT RESULTS ###

#results_infections = read.csv("Infections.csv", stringsAsFactors=F)
#results_prevented = read.csv("Prevented.csv", stringsAsFactors=F)

#output to tif for publication 
tiff("Infections.tif",height=6,width=10,units='in',res=1200) 

#plot infections
barcolors = rev(c(gray.colors(6),"#FFFFFF"))
bars = barplot(matrix(data=results_infections$Mean, nrow=7, ncol=8), col=barcolors, beside=T, ylab="Mean No. Infections", xlab="Additional Prevention Strategy", names.arg=c("None","TasP","Condoms","Seroadaption","TasP &\nCondoms","TasP &\nSeroadaption","Condoms &\nSeroadaption","TasP, Condoms &\nSeroadaption"), cex.names=0.75, ylim=c(0,120))
margins = par()$mar
par(xpd=T, mar=c(1,2,1,2)) 
legend(60, 120, c("0%","1%","5%","10%","15%","20%","25%"), cex=0.9, fill=barcolors, title="% on PrEP", x.intersp=1.5)
par(xpd=F, mar=margins)

#add bars
bars = as.vector(bars)
for (i in 1:length(bars))
{
  arrows(x0=bars[i], y0=(results_infections$Mean[i]-results_infections$SE[i]*1.96), x1=bars[i], y1=(results_infections$Mean[i]+results_infections$SE[i]*1.96), angle=90, length=0.05, code=3)
}
rm(i)

#close file 
dev.off() 

#output to tif for publication 
tiff("Prevented.tif",height=6,width=10,units='in',res=1200) 

#plot prevented by PrEP
barcolors = rev(c(gray.colors(5),"#FFFFFF","#FFFFFF"))
bars = barplot(matrix(data=results_prevented$Mean, nrow=7, ncol=8), col=barcolors, beside=T, ylab="Mean % Infections Prevented", xlab="Additional Prevention Strategy", names.arg=c("None","TasP","Condoms","Seroadaption","TasP &\nCondoms","TasP &\nSeroadaption","Condoms &\nSeroadaption","TasP, Condoms &\nSeroadaption"), cex.names=0.75, ylim=c(0,35))
margins = par()$mar
par(xpd=T, mar=c(1,2,1,2)) 
legend(55, 35, c("1%","5%","10%","15%","20%","25%"), cex=0.8, fill=barcolors[-1], title="% on PrEP", x.intersp=1.5)
par(xpd=F, mar=margins)

#add bars
bars = as.vector(bars)
for (i in 1:length(bars))
{
  arrows(x0=bars[i], y0=(results_prevented$Mean[i]-results_prevented$SE[i]*1.96), x1=bars[i], y1=(results_prevented$Mean[i]+results_prevented$SE[i]*1.96), angle=90, length=0.05, code=3)
}
rm(i)

#close file 
dev.off() 
