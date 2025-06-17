################################################################################
# Title: BMJ_2025_Ryser_Hwang_Relative_Risks.R
# Description: Generates supplementary figures and tables reporting relative risks 
#              of ipsilateral invasive breast cancer (iIBC) in the BMJ DCIS outcomes study.
#
# Outputs include:
#   - Supplementary Figure B: Partial residual plot of age at diagnosis vs. iIBC risk
#   - Supplementary Table E:
#       • Part A: Univariate Cox models for static baseline covariates
#       • Part B: Time-dependent Cox models with fixed endocrine therapy (ET) durations
#       • Part C: Time-dependent Cox models with Monte Carlo-imputed ET durations
#
# Author: Marc D. Ryser
# Contact: marc.ryser@duke.edu
# Affiliation: Duke University
#
# Associated publication:
# Ryser et al. "Cancer outcomes in women without upfront surgery for ductal carcinoma in situ (DCIS)"
# BMJ, 2025.
#
# License: CC BY-NC
#
# Last updated: 2025-05-22
#
# Notes:
# - Input data: analytic_cohort.RData (cleaned and pre-processed cohort)
# - No external output files are saved by default; figures and tables are printed to console or memory
# - Key methods:
#     • Multiple imputation for missing baseline covariates (m=20)
#     • Time-dependent Cox modeling with split intervals for ET exposure
#     • Monte Carlo sampling of ET durations across 1000 replicates for robustness
# - Required R packages:
#     haven, dplyr, survival, mice, mitools, ggplot2, ggsci,
#     ggsurvfit, cmprsk, tidycmprsk, survminer
################################################################################


#-  Load the packages
library(haven)
library(dplyr)
library(survival)
library(ggplot2)
library(ggsci)  # For the JCO palette
library(cmprsk)
library(survminer)
library(ggsurvfit)
library(tidycmprsk)
library(mice)       # For multiple imputation
library(survival)   # For Cox proportional hazards model
library(mitools)    # For pooling results from multiple imputations


#- Load the data
load(file = "analytic_cohort.RData")
dat <- data.all %>%
  filter(is.na(slbe_77777777)) %>%  # Remove those where it is unknown if SLBE occurred
  mutate(fumos_O1.yr = fumos_O1 / 12) %>%
  mutate(fumos_O3.yr = fumos_O3 / 12)

rm(data.all)

#---------------------------------------------
#---------------------------------------------
# PREPARE THE MODELING DATA
#---------------------------------------------
#---------------------------------------------

dat.cox <- dat %>% mutate(Grade.bin = ifelse(
  Grade == "1" | Grade == "2",
  "low",
  ifelse(Grade == "3", "high", "UNK")
)) %>%
  mutate(Age.bin = ifelse(Age < 50, "young", ifelse(Age >= 50 &
                                                      Age < 65, "middle", "old"))) %>%
  mutate(Year = as.numeric(Year)) %>%
  mutate(HT1_days_to_tr = as.numeric(HT1_days_to_tr)) %>%
  mutate(HT2_days_to_tr = as.numeric(HT2_days_to_tr)) %>%
  mutate(
    Comedonecrosis = recode(
      Comedonecrosis,
      "1) Present" = "yes",
      "2) Absent" = "no",
      "3) Unknown" = "UNK"
    )
  ) %>%
  mutate(ET.initiate = ifelse(HT1_days_to_tr > 180 |
                                is.na(HT1_days_to_tr) , "no", "yes")) %>%
  mutate(ID = row_number())



#----------------------------------------------------------
#----------------------------------------------------------
# Supplementary Figure B: Relationship of age at diagnosis
# and risk of ipsilateral invasive breast cancer
#----------------------------------------------------------
#----------------------------------------------------------


# Polynomial model of age
model.age <- coxph(Surv(fumos_O1, eventO1_final) ~ poly(Age,2, raw=TRUE) ,
                   data = dat.cox,
                   #subset = Detection!="UNK",
                   ties = "breslow")
termplot(model.age, xlab="Age at diagnosis (years)", ylab="Partial for poly(Age, 2, raw=TRUE)", se = TRUE, col.term = "blue", col.se = "lightblue")

summary(model.age)
extractAIC(model.age)[2]


#----------------------------------------------------------
#----------------------------------------------------------
# Supplementary Table E: Univariate Cox models
# Part A: Static baseline covariates
#----------------------------------------------------------
#----------------------------------------------------------

set.seed(2001)
# Prepare the dataset
dat.cox.mi <- dat.cox %>%
  select(
    fumos_O1,
    eventO1_final,
    Age,
    HR,
    Race,
    Charlson,
    Grade.bin,
    Detection,
    Year,
    Comedonecrosis,
    ET.initiate
  )

# Convert to NA
dat.cox.mi[dat.cox.mi == "UNK"] <- NA

# Turn into factors
dat.cox.mi <- dat.cox.mi %>%
  mutate(
    HR = as.factor(HR),
    Race = as.factor(Race),
    Charlson = as.factor(Charlson),
    Grade.bin = as.factor(Grade.bin),
    Detection = as.factor(Detection),
    Comedonecrosis = as.factor(Comedonecrosis),
    ET.initiate = as.factor(ET.initiate)
  )

# Step 1: Create the predictor matrix
predictor_matrix <- make.predictorMatrix(dat.cox.mi)

# Step 2: Exclude 'fumos_O1' and 'eventO1_final' as predictors for other variables
predictor_matrix[, c("fumos_O1", "eventO1_final")] <- 0

# Step 3: Prevent 'fumos_O1' and 'eventO1_final' from being imputed
predictor_matrix["fumos_O1", ] <- 0  # No imputation for 'status'
predictor_matrix["eventO1_final", ] <- 0    # No imputation for 'time'

# Step 4: Specify the imputation methods
# Set method for 'fumos_O1' and 'eventO1_final' to "" to exclude them from imputation
methods <- make.method(dat.cox.mi)
methods[c("fumos_O1", "eventO1_final")] <- ""


# Step 5: Generate multiple imputations for missing data
imputed_data <- mice(
  dat.cox.mi,
  m = 20,
  predictorMatrix = predictor_matrix,
  method = methods,
  seed = 1234,
  printFlag = TRUE
)

# Step 6: Fit the Cox model on each imputed dataset
cox_models <- with(
  imputed_data,
  coxph(
    Surv(fumos_O1, eventO1_final) ~  
       poly(Age, 2, raw=TRUE) 
      # Race 
      # Charlson 
      # Year 
      #  Detection
      # HR 
      # Grade.bin 
      # Comedonecrosis
    )
)

# Step 7: Pool the results across imputed data sets
pooled_results <- pool(cox_models)

# Step 8: Summarize pooled results
summary(pooled_results, conf.int = TRUE, exponentiate = TRUE)

# Schoenfeld residual test
schoen<-list()
for(m in 1: length(cox_models$analyses)){
  
  schoen[[m]]<-cox.zph(cox_models$analyses[[m]])
  
}

# Visual assessment of Schoenfeld residual plot
plot(schoen[[5]])


#----------------------------------------------------------
#----------------------------------------------------------
# Supplementary Table E: Univariate Cox models
# Part B: Time-dependent endocrine therapy, fixed duration
#----------------------------------------------------------
#----------------------------------------------------------

set.seed(2001)
#-------------------------------------------------------
#-- STEP 1: Do MI on the baseline covariates
#-------------------------------------------------------


dat.cox.mi <- dat.cox %>%
  select(
    ID,
    fumos_O1,
    eventO1_final,
    Age,
    HR,
    Race,
    Charlson,
    Grade.bin,
    Detection,
    Year,
    Comedonecrosis
  ) %>%
  mutate_if(is.character, ~ na_if(., "UNK")) %>%
  mutate(
    HR = as.factor(HR),
    Race = as.factor(Race),
    Charlson = as.factor(Charlson),
    Grade.bin = as.factor(Grade.bin),
    Detection = as.factor(Detection),
    Comedonecrosis = as.factor(Comedonecrosis)
  )


# Step 1: Create the predictor matrix
predictor_matrix <- make.predictorMatrix(dat.cox.mi)

# Step 2: Exclude 'status' and 'time' as predictors for other variables
predictor_matrix[, c("fumos_O1", "eventO1_final", "ID")] <- 0  # Exclude 'status' and 'time' as predictors

# Step 3: Prevent 'status' and 'time' from being imputed
predictor_matrix["fumos_O1", ] <- 0  # No imputation for 'status'
predictor_matrix["eventO1_final", ] <- 0    # No imputation for 'time'
predictor_matrix["ID", ] <- 0    # No imputation for 'ID'


# Step 4: Specify the imputation methods
# Set method for 'status' and 'time' to "" to exclude them from imputation
methods <- make.method(dat.cox.mi)
methods[c("fumos_O1", "eventO1_final", "ID")] <- ""


# Step 5: Generate multiple imputations for missing data
# Specify the method and number of imputations
imputed_data <- mice(
  dat.cox.mi,
  m = 20,
  predictorMatrix = predictor_matrix,
  method = methods,
  seed = 1234,
  printFlag = TRUE
)  # 'pmm' is predictive mean matching


#- Step 6: Extract a list of the imputed data frames
#- Include the non-imputed one (first data set)
all_imputed <- complete(imputed_data, action = "all", include = TRUE)

#-------------------------------------------------------
#-- STEP 2: Attach the time-dependent ET data
#-------------------------------------------------------

all_imputed.ET <- list()

# Choose the mid-points of the ET.duration intervals
gurke <- c(6, 18, 3.5 * 12, 7.5 * 12, 5 * 12)


# Loop over all datasets, including the unimputed one
for (m in 1:length(all_imputed)) {
  
  # Attach the ET data
  dat.temp <- all_imputed[[m]] %>%
    left_join(dat.cox %>% select(ID, HT1_days_to_tr, ET.duration), by = "ID") %>%
    mutate(
      time = fumos_O1,
      event = eventO1_final,
      ET.start = HT1_days_to_tr / 30.44, # convert to month
      ET.duration = recode(
        ET.duration,
        "1) Less than 1 year" = gurke[1],
        "2) 1-2 years" = gurke[2],
        "3) 2-5 years" = gurke[3],
        "4) Greater than 5 years" = gurke[4],
        "5) Duration unknown" = gurke[5],
        .default = -1
      ),
      ET.end = ifelse(!is.na(ET.start), pmin(ET.start +ET.duration, time), NA)
    ) %>%
    mutate(
      tstart = NA,
      tstop = NA,
      ET.status = NA,
      event.status = NA
    )
  
  # Isolate patients without ET during follow-up
  dat.out <- dat.temp %>% filter(is.na(ET.start) | ET.start >= time)
  
  dat.out$tstart = 0
  dat.out$tstop = dat.out$time
  dat.out$event.status = dat.out$event
  dat.out$ET.status = 0
  
  # Isolate patients with ET during follow-up
  ind <- which(dat.temp$ET.start < dat.temp$time)
  
  # record the results for the 409 patients with ET during f/u
  results <- vector("list", length(ind))
  
  for (i in seq_along(ind)) {
    k <- ind[i]
    # Replicate the row three times without `dplyr`
    aha <- dat.temp[k, ]
    aha <- aha[rep(1, 3), ]
    
    # Assign the new values to `aha` directly
    aha$tstart <- c(0, aha$ET.start[1], aha$ET.end[1])
    aha$tstop <- c(aha$ET.start[1], aha$ET.end[1], aha$time[1])
    aha$ET.status <- c(0, 1, 0)
    aha$event.status <- c(0,  ifelse(aha$ET.end[1] == aha$time[1], aha$event[1], 0), aha$event[1])
    
    # Filter out rows where `tstart` equals `tstop`
    aha <- aha[aha$tstart != aha$tstop, ]
    
    # Store in list instead of binding within the loop
    results[[i]] <- aha
  }
  
  # Combine all results at once
  argh <- do.call(rbind, results)
  all_imputed.ET[[m]] <- rbind(dat.out, argh)
  
  
}


# Add two columns for .immp and .id
data_list <- lapply(seq_along(all_imputed.ET), function(i) {
  # Add .imp column as the index - 1 (to label from 0 to M)
  all_imputed.ET[[i]]$.imp <- i - 1
  # Add .id column as the row number within each data frame
  all_imputed.ET[[i]]$.id <- seq_len(nrow(all_imputed.ET[[i]]))
  # Return the modified data frame
  all_imputed.ET[[i]]
})

# Combine the list into a single long data frame
all_imputed_long <- bind_rows(data_list)

# Apply the ET only to those with positive HR 
all_imputed_long <- all_imputed_long %>%
  mutate(ET.status.restricted = ifelse(HR == "positive", ET.status, 0))


# Convert the long data frame into a mids object
all_imputed_mids <- as.mids(all_imputed_long)



#-------------------------------------------------------
#-- STEP 3: Run the Cox models
#-------------------------------------------------------

#- Pooled Cox  analysis
cox_models <- with(
  all_imputed_mids,
  coxph(
    Surv(tstart, tstop, event.status) ~
     #poly(Age, 2, raw=TRUE)+
      ET.status
     # HR +
      #Grade.bin +
     # Detection +
      #Race +
      #Charlson+
      #Year+
      #Comedonecrosis
  )
)
pooled_results <- pool(cox_models)
out<-summary(pooled_results, conf.int = TRUE, exponentiate = TRUE) 


out <- out %>% select(term, estimate, `2.5 %`, `97.5 %`,p.value) %>%
  mutate(estimate = round(estimate,2)) %>%
  mutate(p.value = round(p.value,2)) %>%
  mutate(`2.5 %` = round(`2.5 %`,2)) %>%
  mutate(`97.5 %` = round(`97.5 %`,2)) 






#----------------------------------------------------------
#----------------------------------------------------------
# Supplementary Table E: Univariate Cox models
# Part C: Time-dependent endocrine therapy, variable duration (using Monte Carlo)
#----------------------------------------------------------
#----------------------------------------------------------

set.seed(2001)
#-------------------------------------------------------
#-- STEP 1: Do MI on the baseline covariates
#-------------------------------------------------------

dat.cox.mi <- dat.cox %>%
  select(
    ID,
    fumos_O1,
    eventO1_final,
    Age,
    HR,
    Race,
    Charlson,
    Grade.bin,
    Detection,
    Year,
    Comedonecrosis
  ) %>%
  mutate_if(is.character, ~ na_if(., "UNK")) %>%
  mutate(
    HR = as.factor(HR),
    Race = as.factor(Race),
    Charlson = as.factor(Charlson),
    Grade.bin = as.factor(Grade.bin),
    Detection = as.factor(Detection),
    Comedonecrosis = as.factor(Comedonecrosis)
  )


# Step 1: Create the predictor matrix
predictor_matrix <- make.predictorMatrix(dat.cox.mi)

# Step 2: Exclude 'status' and 'time' as predictors for other variables
predictor_matrix[, c("fumos_O1", "eventO1_final", "ID")] <- 0  # Exclude 'status' and 'time' as predictors

# Step 3: Prevent 'status' and 'time' from being imputed
predictor_matrix["fumos_O1", ] <- 0  # No imputation for 'status'
predictor_matrix["eventO1_final", ] <- 0    # No imputation for 'time'
predictor_matrix["ID", ] <- 0    # No imputation for 'ID'


# Step 4: Specify the imputation methods
# Set method for 'status' and 'time' to "" to exclude them from imputation
methods <- make.method(dat.cox.mi)
methods[c("fumos_O1", "eventO1_final", "ID")] <- ""


# Step 5: Generate multiple imputations for missing data
# Specify the method and number of imputations
imputed_data <- mice(
  dat.cox.mi,
  m = 20,
  predictorMatrix = predictor_matrix,
  method = methods,
  seed = 1234,
  printFlag = TRUE
)  # 'pmm' is predictive mean matching



#- Step 6: Extract a list of the imputed data frames
#- Include the non-imputed one (first data set)
all_imputed <- complete(imputed_data, action = "all", include = TRUE)


#-------------------------------------------------------
#-- STEP 2: Attach the time-dependent ET data
#-------------------------------------------------------

# Save the MC runs
spit.out <- list()


for (Q in 1:1000) {
  print(Q)
  
  # Save the imputed data sets
  all_imputed.ET <- list()
  
  # sample from the intervals
  gurke.1 <- runif(sum(dat.cox$ET.duration == "1) Less than 1 year"), 0, 12)
  gurke.2 <- runif(sum(dat.cox$ET.duration == "2) 1-2 years"), 12, 24)
  gurke.3 <- runif(sum(dat.cox$ET.duration == "3) 2-5 years"), 24, 60)
  gurke.4 <- runif(sum(dat.cox$ET.duration == "4) Greater than 5 years"), 60, 120)
  gurke.5 <- runif(sum(dat.cox$ET.duration == "5) Duration unknown"), 0, 120)

  
  for (m in 1:length(all_imputed)) {
    dat.temp <- all_imputed[[m]] %>%
      left_join(dat.cox %>% select(ID, HT1_days_to_tr, ET.duration), by = "ID") %>%
      mutate(time = fumos_O1,
             event = eventO1_final,
             ET.start = HT1_days_to_tr / 30.44)
    
    #- Insert the ET durations
    dat.temp$ET.duration[dat.temp$ET.duration == "1) Less than 1 year"] <-  gurke.1
    dat.temp$ET.duration[dat.temp$ET.duration == "2) 1-2 years"] <- gurke.2
    dat.temp$ET.duration[dat.temp$ET.duration == "3) 2-5 years"] <- gurke.3
    dat.temp$ET.duration[dat.temp$ET.duration == "4) Greater than 5 years"] <-   gurke.4
    dat.temp$ET.duration[dat.temp$ET.duration == "5) Duration unknown"] <-   gurke.5
    dat.temp$ET.duration <- as.numeric(dat.temp$ET.duration)
    
    dat.temp <- dat.temp %>%
      mutate(ET.end = ifelse(!is.na(ET.start), pmin(ET.start + ET.duration, time), NA)) %>%
      mutate(
        tstart = NA,
        tstop = NA,
        ET.status = NA,
        event.status = NA
      )
    
  
    
    # Isolate patients without ET during follow-up
    dat.out <- dat.temp %>% filter(is.na(ET.start) | ET.start >= time)
    
    dat.out$tstart = 0
    dat.out$tstop = dat.out$time
    dat.out$event.status = dat.out$event
    dat.out$ET.status = 0
    
    # Isolate patients with ET during follow-up
    ind <- which(dat.temp$ET.start < dat.temp$time)
    
    results <- vector("list", length(ind))
    
    for (i in seq_along(ind)) {
      k <- ind[i]
      # Replicate the row three times without `dplyr`
      aha <- dat.temp[k, ]
      aha <- aha[rep(1, 3), ]
      
      # Assign the new values to `aha` directly
      aha$tstart <- c(0, aha$ET.start[1], aha$ET.end[1])
      aha$tstop <- c(aha$ET.start[1], aha$ET.end[1], aha$time[1])
      aha$ET.status <- c(0, 1, 0)
      aha$event.status <- c(0,
                            ifelse(aha$ET.end[1] == aha$time[1], aha$event[1], 0),
                            aha$event[1])
      
      # Filter out rows where `tstart` equals `tstop`
      aha <- aha[aha$tstart != aha$tstop, ]
      
      # Store in list instead of binding within the loop
      results[[i]] <- aha
    }
    
    # Combine all results at once
    argh <- do.call(rbind, results)
    all_imputed.ET[[m]] <- rbind(dat.out, argh)
    
    
  }
  

  
  data_list <- lapply(seq_along(all_imputed.ET), function(i) {
    # Add .imp column as the index - 1 (to label from 0 to M)
    all_imputed.ET[[i]]$.imp <- i - 1
    # Add .id column as the row number within each data frame
    all_imputed.ET[[i]]$.id <- seq_len(nrow(all_imputed.ET[[i]]))
    # Return the modified data frame
    all_imputed.ET[[i]]
  })
  
  # Combine the list into a single long data frame
  all_imputed_long <- bind_rows(data_list) %>%
    mutate(ET.status.restricted = ifelse(HR == "positive", ET.status, 0))
  
  
  # Convert the long data frame into a mids object
  all_imputed_mids <- as.mids(all_imputed_long)
  
  #- pooled analysis
  cox_models <- with(
    all_imputed_mids,
    coxph(
      Surv(tstart, tstop, event.status) ~ ET.status 
    )
  )
  
  # Step 3: Pool the results across imputed datasets
  pooled_results <- pool(cox_models)
  
  # Step 4: Summarize pooled results
  #summary(pooled_results, conf.int = TRUE, exponentiate = TRUE)
  
  spit.out[[Q]] <- summary(pooled_results,
                            conf.int = TRUE,
                            exponentiate = TRUE)
  
}

# Use lapply to extract the 'X' column from each data frame, then bind them together as columns
#- Calculate the proportion of MC runs where time-dependent ET was a significant (<.05) predictor
X_combined <- do.call(cbind, lapply(spit.out, function(df) df$p.value))

rownames(X_combined) <- spit.out[[3]]$term

rowSums(X_combined <= .05) / dim(X_combined)[2]


  