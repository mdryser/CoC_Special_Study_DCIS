################################################################################
# Title: aux.fun.R
# Description: Auxiliary functions for cohort verification and survival analysis
#              in the BMJ DCIS outcomes study. Includes custom functions for:
#              - Extracting absolute risks from Kaplan-Meier fits
#              - Disease-specific survival estimation
#              - Modeling expected other-cause mortality using US life tables
#
# Authors: Marc D. Ryser
# Contact: marc.ryser@duke.edu
# Affiliation: Duke University 
#
# Associated publication:
# Ryser et al. "Cancer outcomes in women without upfront surgery for ductal carcinoma in situ (DCIS)"
# BMJ, 2025
#
# License: CC BY-NC
#
# Last updated: 2025-6-17
#
# Notes:
# - Used in: 1_Data_check_FINAL.R and downstream survival analysis scripts
# - Requires mortality rate file: MortalityRatesTables_17Mar2017.csv*
# - Key functions:
#     • Risk.KM(): Computes absolute risk with CIs from Kaplan-Meier estimates
#     • DSS.KM(): Computes disease-specific survival with CIs
#     • surv.fun(): Computes survival under other-cause mortality hazards
#     • OC.surv(): Computes cohort-level expected mortality based on age and year
# *Data from Gangnon, Ronald E., et al. "Contribution of breast cancer to overall mortality for US women." Medical Decision Making 38.1_suppl (2018): 24S-31S.
################################################################################


#-------------------------------------------------------------------------
#-- Cumulative incidence (absolute risk); using Kaplan-Meier estimators
#-------------------------------------------------------------------------

Risk.KM<-function(km_fit, time_points){
  
  # Extract the survival information for each time point for both risk groups
  summary_km <- summary(km_fit, times = time_points)
  
  # Calculate the cumulative risks (1 - survival probability) and convert to percentages
  
  
  if(is.null(summary_km$strata)){
    
    risk_df <- data.frame(
      Time = summary_km$time,
      #N=summary_km$n,
      At.risk = summary_km$n.risk,
      # Risk_group = summary_km$strata,
      Risk = round(100*(1 - summary_km$surv),1),  # Convert to risk percentages
      Lower_CI = round(100*(1 - summary_km$upper),1),  # Convert to risk percentages
      Upper_CI = round(100*(1 - summary_km$lower),1)  # Convert to risk percentages
    )
    
  } else {
    
    risk_df <- data.frame(
      Time = summary_km$time,
      #N=summary_km$n,
      At.risk = summary_km$n.risk,
      Risk_group = summary_km$strata,
      Risk = round(100*(1 - summary_km$surv),1),  # Convert to risk percentages
      Lower_CI = round(100*(1 - summary_km$upper),1),  # Convert to risk percentages
      Upper_CI = round(100*(1 - summary_km$lower),1)  # Convert to risk percentages
    )
    
  }
  
  risk_df <- risk_df %>%
    mutate(CI_95 = paste0(round(Lower_CI, 1), "-", round(Upper_CI, 1))) %>%
    select(-Lower_CI, -Upper_CI) %>%
    mutate()
  
  
  # View the extracted risk values with CIs
  return(risk_df)
  
}


#-----------------------------------------------------------------
# Disease-specific survival using Kaplan-Meier estimators
#-----------------------------------------------------------------

DSS.KM<-function(km_fit, time_points){
  
  # Extract the survival information for each time point for both risk groups
  summary_km <- summary(km_fit, times = time_points)
  
  # Calculate the cumulative risks (1 - survival probability) and convert to percentages
  
  
  if(is.null(summary_km$strata)){
    
    risk_df <- data.frame(
      Time = summary_km$time,
      #N=summary_km$n,
      At.risk = summary_km$n.risk,
      # Risk_group = summary_km$strata,
      Risk = round(100*(summary_km$surv),1),  # Convert to risk percentages
      Lower_CI = round(100*(summary_km$lower),1),  # Convert to risk percentages
      Upper_CI = round(100*(summary_km$upper),1)  # Convert to risk percentages
    )
    
  } else {
    
    risk_df <- data.frame(
      Time = summary_km$time,
      #N=summary_km$n,
      At.risk = summary_km$n.risk,
      Risk_group = summary_km$strata,
      Risk = round(100*(summary_km$surv),1),  # Convert to risk percentages
      Lower_CI = round(100*(summary_km$lower),1),  # Convert to risk percentages
      Upper_CI = round(100*(summary_km$upper),1)  # Convert to risk percentages
    )
    
  }
  
  risk_df <- risk_df %>%
    mutate(CI_95 = paste0(round(Lower_CI, 1), "-", round(Upper_CI, 1))) %>%
    select(-Lower_CI, -Upper_CI) %>%
    mutate()
  
  
  # View the extracted risk values with CIs
  return(risk_df)
  
}


#-----------------------------------------------------------------
# Auxiliary function for general population other cause mortality
#-----------------------------------------------------------------


# Hazard for other cause mortality

lambda_oc<-function(t, OChaz){
  
  ft<- floor(t)
  out=ifelse(ft<=119, OChaz[ft+1],OChaz[120])
  
  return(out)
  
}

lambda_oc<-Vectorize(lambda_oc, "t")


# Other cause mortality survival function 

surv.fun <- function(age0, t, OChaz){
  
  a<-age0
  b<-floor(t)-1
  
  if(b>=a){ 
    
    v=apply(t(a:b),1,lambda_oc, OChaz=OChaz)
    
  } else {
    
    v=0
    
  }
  
  hh<-sum(v) + (t-floor(t))*lambda_oc(t, OChaz=OChaz)
  
  S<-exp(-hh)
  
  return(S)
}

surv.fun<-Vectorize(surv.fun,vectorize.args = "t")



#- Load other cause mortality data
OC.data <- read.csv(file="MortalityRatesTables_17Mar2017.csv")


#--Expected Cohort Survival 
#- time: from diagnosis (years)
#- age.vec: ages of participants at diagnosis (years)
#- year.vec: year of diagnosis 


OC.surv <- function(time, age.vec, year.vec) {
  
  cohort.vec <- year.vec-age.vec
  
  surv.vec <- rep(0, length(cohort.vec))
  for(c in 1:length(age.vec)) {
    
    ind<-which(OC.data$Cohort==max(1950, cohort.vec[c]))
    surv.vec[c]<-surv.fun(age.vec[c], age.vec[c]+time, OC.data$Estimate[ind]/100000)
    
  }
  
  return(1-mean(surv.vec))
  
}

OC.surv<-Vectorize(OC.surv,vectorize.args = "time")
