################################################################################
# Title: BMJ_2025_Ryser_Hwang_Absolute_Risks.R
# Description: Generates primary absolute risk figures and tables for the BMJ DCIS outcomes study.
# Outputs include:
#   - Table 1: Patient and tumor characteristics
#   - Figure 2A: Cumulative incidence of ipsilateral invasive breast cancer (iIBC), overall
#   - Figure 2B: Cumulative incidence of iIBC, by risk group
#   - Figure 3A: Disease-specific survival (DSS), overall
#   - Figure 3B: DSS, by risk group
#   - Supplementary Figure A: Time to first ipsilateral surgery by type (LX vs MX)
#   - Supplementary Figure C:
#         - Panel A: Competing risks (breast cancer death vs other-cause death) + general population reference
#         - Panel B: Competing risks stratified by risk group
#
# Authors: Marc D. Ryser
# Contact: marc.ryser@duke.edu
# Affiliation: Duke University
#
# Associated publication:
# Ryser et al. "Cancer outcomes in women without upfront surgery for ductal carcinoma in situ (DCIS)"
# BMJ, 2025.
#
# License: CC BY-NC
#
# Last updated: 2025-06-17
#
# Notes:
# - Input data: analytic_cohort.RData (cleaned cohort); aux.fun.R (utility functions)
# - Output files:
#       - table1.docx
#       - Figure_2A.pdf, Figure_2B.pdf
#       - Figure_3A.pdf, Figure_3B.pdf
#       - Supp_Figure_A.pdf
#       - Supp_Figure_C_panel_A.pdf, Supp_Figure_C_panel_B.pdf
# - Required R packages:
#       haven, dplyr, survival, ggplot2, ggsci, cmprsk, survminer,
#       ggsurvfit, tidycmprsk, gtsummary, flextable
################################################################################




# Packages
library(haven)
library(dplyr)
library(survival)
library(ggplot2)
library(ggsci)  
library(cmprsk)
library(survminer)
library(ggsurvfit)
library(tidycmprsk)
library(gtsummary)


# JCO color palette for the plot
jco_colors <- pal_jco()(4)  # JCO palette with two colors for the two risk groups

# Extract survival probabilities and CIs at specific time points (2, 5, 8 years)
time_points <- c(0, 2, 5, 8)

# Month day conversion
fmo<-30.4375



#-- Aux functions
source('aux.fun.R')

load(file = "analytic_cohort.RData")


#---------------------------------------------
# TABLE 1: Covariate summary
#---------------------------------------------

table1 <- data.all %>%
  select(Age, hispanic, Race,Facility, Comedonecrosis, metro, noHSdiploma, medianHHI, location, insured, Year,  Charlson, Grade, HR, Detection, Risk.group, Biopsy, ET, ET.duration) %>%
  tbl_summary(statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", 
                               all_categorical() ~ "{n} ({p}%)"),
              digits = all_categorical() ~ 3  ) 

table1 %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "table1.docx")  

# Follow-up using inverse KM method

km_fit <- survfit(Surv(fumos_O1.yr, eventO1_final == 0) ~ 1, data = data.all)
km_summary <- summary(km_fit)
median_time <- km_summary$table["median"]
lower_ci <- km_summary$table["0.95LCL"]
upper_ci <- km_summary$table["0.95UCL"]

list(
  Median_Survival_Time = 12*median_time,
  Lower_95_CI = 12*lower_ci,
  Upper_95_CI = 12*upper_ci
)

ggsurvfit(km_fit,
          type="risk")+
  add_risktable() + 
  add_confidence_interval(fill=jco_colors[1])+
    scale_x_continuous(breaks =  c(0,2.5,5,7.5,10))+
  ylab("Cumulative risk of iIBC")+
  xlab("Time since diagnosis (years)")


#---------------------------------------------
# FIGURE 2A: Kaplan-Meier for iIBC 
#---------------------------------------------

# KM risk estimates
km_fit <- survfit(Surv(fumos_O1.yr, eventO1_final == 1) ~ 1, data = data.all)
results <- Risk.KM(km_fit, time_points)
print(results)

# Extract survival estimates
surv_data <- data.frame(
  time = km_fit$time,
  surv = 1 - km_fit$surv # Convert survival probability to cumulative incidence
)

# Filter to exclude times before 0.5 years
surv_data <- surv_data[surv_data$time >= 0.5, ]

# Create the plot
p <- ggsurvfit(km_fit, type = "risk", color = jco_colors[3], alpha=0) +
  add_risktable() + 
  add_confidence_interval(fill = jco_colors[3]) +
  ylim(0, 1) +
  scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10)) +
  ylab("Cumulative incidence of iIBC") +
  xlab("Time since diagnosis (years)") +
  xlim(0, 8) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.margin = margin(20, 0, 30, 0)
  ) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray", linewidth = 1) + # Vertical dashed line at 0.5 years
  geom_step(data = surv_data, aes(x = time, y = surv), color = jco_colors[3], linewidth = .5) # Manually plot curve from 0.5 years

p


ggsave("Figure_2A.pdf", plot = p, device = "pdf", onefile = FALSE,
       width = 10, height = 7, dpi = 300)


#------------------------------------------------
# FIGURE 2B: Kaplan-Meier for iIBC, by Risk group
#------------------------------------------------


# Kaplan-Meier risks by risk.group (excluding unknowns)
km_fit <- survfit(Surv(fumos_O1.yr, eventO1_final == 1) ~ Risk.group, data = data.all, subset = Risk.group != "UNK")
results <- Risk.KM(km_fit, time_points)
print(results)


# p-value by log-rank
surv_pvalue(
  km_fit,
  method = "survdiff")


length.HR<-km_fit$strata[1]
length.LR<-km_fit$strata[2]
 surv_data <- data.frame(
   time=c(km_fit$time[1:length.HR], km_fit$time[(length.HR+1): (length.HR+length.LR)]),
   surv= c(1-km_fit$surv[1:length.HR],1-km_fit$surv[(length.HR+1): (length.HR+length.LR)]),
   group=c(rep("High-risk", km_fit$strata[1]), rep("Low-risk", km_fit$strata[2]))
   
 )

# Filter to exclude times before 0.5 years
surv_data <- surv_data[surv_data$time >= 0.5, ]

options("ggsurvfit.switch-color-linetype" = FALSE)

q<-survfit2(Surv(fumos_O1.yr, eventO1_final == 1) ~ Risk.group, data = data.all, subset = Risk.group != "UNK") %>%
  ggsurvfit(type="risk", alpha=0)+
  add_risktable(risktable_height=.2) + 
  add_confidence_interval()+
  ylim(0,.4)+  
  scale_x_continuous(breaks =  c(0,2.5,5,7.5,10))+xlim(0,8)+
  ylab("Cumulative incidence of iIBC")+
  xlab("Time since diagnosis (years)")+
  scale_color_manual(values = jco_colors) +
  scale_fill_manual(values = jco_colors) +
  theme(
    text = element_text(size = 16),       # Increase base font size
    axis.title = element_text(size = 16),  # Axis titles
    axis.text = element_text(size = 14),   # Axis tick labels
    plot.title = element_text(size = 18),  # Plot title
    legend.title = element_text(size = 14), # Legend title
    legend.text = element_text(size = 12),   # Legend text
    plot.margin = margin(20, 0, 30, 0))+ 
  geom_vline(xintercept = 0.5,
             linetype = "dashed", color = "gray", linewidth = 1)  +  # Vertical dashed line at x = 5
  geom_step(data = surv_data, aes(x = time, y = surv,color=group), linewidth = .5) # Manually plot curve from 0.5 years



q

ggsave("Figure_2B.pdf", plot = q, device = "pdf", onefile = FALSE,
       width = 10, height = 7, dpi = 300)



#---------------------------------------------
# FIGURE 3A: Kaplan-Meier for DSS
#---------------------------------------------


km_fit <- survfit(Surv(fumos_O3.yr, eventO3_final == 1) ~ 1, data = data.all)

# Risk estimates
results <- DSS.KM(km_fit, time_points)
print(results)

# Extract survival estimates
surv_data <- data.frame(
  time = km_fit$time,
  surv = km_fit$surv # Convert survival probability to cumulative incidence
)

# Filter to exclude times before 0.5 years
surv_data <- surv_data[surv_data$time >= 0.5, ]


p2<-ggsurvfit(km_fit,
             color = jco_colors[3], alpha=0)+
  add_risktable() + 
  add_confidence_interval(fill=jco_colors[3])+
  ylim(0,1)+  scale_x_continuous(breaks =  c(0,2.5,5,7.5,10))+
  ylab("Disease-specific survival probability")+
  xlab("Time since diagnosis (years)") +
  xlim(0,8) +
  theme(
    text = element_text(size = 16),       # Increase base font size
    axis.title = element_text(size = 16),  # Axis titles
    axis.text = element_text(size = 14),   # Axis tick labels
    plot.title = element_text(size = 18),  # Plot title
    legend.title = element_text(size = 14), # Legend title
    legend.text = element_text(size = 12),   # Legend text
    plot.margin = margin(20, 0, 30, 0)) + 
  geom_vline(xintercept = 0.5,
             linetype = "dashed", color = "gray", linewidth = 1) +  # Vertical dashed line at x = 5
  geom_step(data = surv_data, aes(x = time, y = surv), color = jco_colors[3], linewidth = .5) # Manually plot curve from 0.5 years




p2

ggsave("Figure_3A.pdf", plot = p2, device = "pdf", onefile = FALSE,
       width = 10, height = 7, dpi = 300)



#---------------------------------------------
# FIGURE 3B: Kaplan-Meier for DSS, by Risk group
#---------------------------------------------

# Kaplan-Meier by risk.group (excluding unknowns)
km_fit <- survfit(Surv(fumos_O3.yr, eventO3_final == 1) ~ Risk.group, data = data.all, subset = Risk.group != "UNK")
results <- DSS.KM(km_fit, time_points)
print(results)

# p-value by log-rank
surv_pvalue(
  km_fit,
  method = "survdiff")





length.HR<-km_fit$strata[1]
length.LR<-km_fit$strata[2]
surv_data <- data.frame(
  time=c(km_fit$time[1:length.HR], km_fit$time[(length.HR+1): (length.HR+length.LR)]),
  surv= c(km_fit$surv[1:length.HR],km_fit$surv[(length.HR+1): (length.HR+length.LR)]),
  group=c(rep("High-risk", km_fit$strata[1]), rep("Low-risk", km_fit$strata[2]))
  
)

# Filter to exclude times before 0.5 years
surv_data <- surv_data[surv_data$time >= 0.5, ]



options("ggsurvfit.switch-color-linetype" = FALSE)

q2<-survfit2(Surv(fumos_O3.yr, eventO3_final == 1) ~ Risk.group, data = data.all, subset = Risk.group != "UNK") %>%
  ggsurvfit(type="survival", alpha=0)+
  add_risktable(risktable_height=.2) + 
  add_confidence_interval()+
  ylim(.6,1)+  
  scale_x_continuous(breaks =  c(0,2.5,5,7.5,10))+xlim(0,8)+
  ylab("Disease-specific survival probability")+
  xlab("Time since diagnosis (years)")+
  scale_color_manual(values = jco_colors) +
  scale_fill_manual(values = jco_colors) +
  theme(
    text = element_text(size = 16),       # Increase base font size
    axis.title = element_text(size = 16),  # Axis titles
    axis.text = element_text(size = 14),   # Axis tick labels
    plot.title = element_text(size = 18),  # Plot title
    legend.title = element_text(size = 14), # Legend title
    legend.text = element_text(size = 12),   # Legend text
    plot.margin = margin(20, 0, 30, 0))+ 
  geom_vline(xintercept = 0.5,
             linetype = "dashed", color = "gray", linewidth = 1) + # Vertical dashed line at x = 5
  geom_step(data = surv_data, aes(x = time, y = surv,color=group), linewidth = .5) # Manually plot curve from 0.5 years


q2

ggsave("Figure_3B.pdf", plot = q2, device = "pdf", onefile = FALSE,
       width = 10, height = 7, dpi = 300)



#-----------------------------------------------------------
# SUPPLEMENTARY FIGURE A: Time to first ipsilateral surgery
#-----------------------------------------------------------

dat.surg <- data.all %>% 
  mutate(First.ipsi.Surgery.type=recode(First.ipsi.Surgery.type, "1) Lumpectomy/Quadrantectomy/Partial mastectomy"="LX", "2) Mastectomy"="MX")) %>%
  mutate(First.ipsi.Surgery.type=as.factor(First.ipsi.Surgery.type)) %>%
  filter(!is.na(First.ipsi.Surgery.type)) %>%
  mutate(First.ipsi.Surgery.months = First.ipsi.surgery.days / fmo) %>% # Convert days to months
  select(First.ipsi.Surgery.months, First.ipsi.Surgery.type)


# Create a histogram with different colors for the two types
t<-ggplot(dat.surg, aes(x = First.ipsi.Surgery.months, fill = First.ipsi.Surgery.type)) +
  geom_histogram(binwidth = 2, color = "white", position = "stack", alpha = 0.8) +  # Adjust 'bins' as needed
  scale_fill_manual(values = c("LX" = "#0073C2FF", "MX" = "#EFC000FF")) +  # Customize colors for LX and MX
  labs(
    #title = "Time from Diagnosis to Surgery by Surgery Type",
    x = "Time from diagnosis to surgery (months)",
    y = "Frequency (n)",
    fill = "Surgery Type"
  ) +
  scale_x_continuous(
    limits = c(0, 120),                   # Set axis limits here
    breaks = seq(0, 120, by = 6)          # Tick every 6 months
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  ) 
t



ggsave("Supp_Figure_A.pdf", plot = t, device = "pdf", onefile = FALSE,
       width = 10, height = 5, dpi = 300)




#----------------------------------------------------------------
# SUPPLEMENTARY FIGURE C, panel A: Competing risk analysis of DSS and OS
#----------------------------------------------------------------

#- Expected other cause survival in general population

x<- seq(0, 10, .1) # 10-year horizon

#- Assuming we start at Age+.5 (midpoint), we add another .5 because we condition
#- on survival for 6 months, so the clock starts at Age+1
#- But Age+1 is at the 6-month mark on the x-axis, so we need to shift the whole thing to the right

OC.predict<-data.frame(time=x+.5, rate=OC.surv(x, data.all$Age+1, as.numeric(data.all$Year)))

ggplot(data=OC.predict, aes(x=time, y=rate)) + geom_point()

#- Competing risks in our cohort

dat.comp <- data.all %>%
  mutate(surv_type = as.factor(recode(O3_type, "Censored at Date of Last Contact" = "Censored", "Breast Cancer Death" = "BC death", "Unrelated Death" = "OC death"))) %>%
  mutate(surv_type = factor(surv_type, levels = c("Censored", "BC death", "OC death"))) 

cuminc_fit <- tidycmprsk::cuminc(Surv(fumos_O3.yr, surv_type) ~ 1, dat.comp)

surv_data <- cuminc_fit %>%
  tidycmprsk::tidy() %>%  
  select(time, estimate, outcome)

# Filter to exclude times before 0.5 years
surv_data <- surv_data[surv_data$time >= 0.5, ]

# Competing risks analysis for DSS and OS, stratified by Risk.group
options("ggsurvfit.switch-color-linetype" = TRUE)

r<-tidycmprsk::cuminc(Surv(fumos_O3.yr, surv_type) ~ 1, dat.comp) %>%
  ggcuminc(outcome = c("BC death", "OC death"), alpha=0) +
  scale_color_manual(values = c(jco_colors[1], jco_colors[2])) +
  scale_fill_manual(values = c(c(jco_colors[1], jco_colors[2])))+
  add_risktable(risktable_group="strata") + 
  add_confidence_interval()+
  ylim(0,.4) + xlim(0,8)+ 
  xlab("Time since diagnosis (years)") +  # Change x-axis label here
  ylab("Cumulative incidence") +
  theme(
    text = element_text(size = 16),       # Increase base font size
    axis.title = element_text(size = 16),  # Axis titles
    axis.text = element_text(size = 14),   # Axis tick labels
    plot.title = element_text(size = 18),  # Plot title
    legend.title = element_text(size = 14), # Legend title
    legend.text = element_text(size = 12),   # Legend text
    plot.margin = margin(20, 0, 30, 0)) + 
  geom_vline(xintercept = 0.5,
             linetype = "dashed", color = "gray", linewidth = 1) +  # Vertical dashed line at x = 5
  geom_step(data = surv_data, aes(x = time, y = estimate,color=outcome), linewidth = .5) # Manually plot curve from 0.5 years

  
#- Add the predicted survival in general population
r<- r+   geom_line(data = OC.predict, aes(x = time, y=rate), color = "black")  # Line for df2
r

#- Save the figure
ggsave("Supp_Figure_C_panel_A.pdf", plot = r, device = "pdf", onefile = FALSE,
       width = 10, height = 7, dpi = 300)



# Competing risks analysis using cuminc()
cif_fit <- cmprsk::cuminc(ftime = dat.comp$fumos_O3.yr, fstatus = dat.comp$surv_type, cencode="Censored")

cif_timepoints <- timepoints(cif_fit, times = time_points)

cif_out.PE <- round(100*cif_timepoints$est,1)
cif_out.lower <- round(100*(cif_timepoints$est-1.96*sqrt(cif_timepoints$var)),1)
cif_out.upper <- round(100*(cif_timepoints$est+1.96*sqrt(cif_timepoints$var)),1)





#---------------------------------------------
# SUPPLEMENTARY FIGURE C, panel B: DSS and OS
# Low-risk vs high-risk
#---------------------------------------------

# Filter out "UNK" values in Risk.group
dat_filtered <- dat.comp %>%
  filter(Risk.group %in% c("Low-risk", "High-risk")) 

# Competing risks analysis for DSS and OS, stratified by Risk.group

s<-tidycmprsk::cuminc(Surv(fumos_O3.yr, surv_type) ~ Risk.group, dat_filtered) %>%
  ggcuminc(outcome = c("BC death", "OC death")) +
  add_risktable() + 
 # add_confidence_interval() +
  scale_color_manual(values = c(jco_colors[1], jco_colors[2])) +
  scale_fill_manual(values = c(c(jco_colors[1], jco_colors[2])))+
  ylim(0,.4)+
  xlim(0,8)+ 
  xlab("Time since diagnosis (years)") +  # Change x-axis label here
  ylab("Cumulative incidence")+
  theme(
    text = element_text(size = 16),       # Increase base font size
    axis.title = element_text(size = 16),  # Axis titles
    axis.text = element_text(size = 14),   # Axis tick labels
    plot.title = element_text(size = 18),  # Plot title
    legend.title = element_text(size = 14), # Legend title
    legend.text = element_text(size = 12),   # Legend text
    plot.margin = margin(20, 0, 30, 0))+ 
  geom_vline(xintercept = 0.5,
             linetype = "dashed", color = "gray", linewidth = 1)   # Vertical dashed line at x = 5


s 

ggsave("Supp_Figure_C_pal_B.pdf", plot = s, device = "pdf", onefile = FALSE,
       width = 10, height = 7, dpi = 300)
                  
 # Competing risks analysis using cuminc()
cif_fit <- cmprsk::cuminc(ftime = dat_filtered$fumos_O3.yr, fstatus = dat_filtered$surv_type, group = dat_filtered$Risk.group, cencode="Censored")

# Use timepoints() to extract estimates and variances at a specific time (e.g., time = 96)
cif_timepoints <- timepoints(cif_fit, times = c(2,5,8))

cif_out.PE <- round(100*cif_timepoints$est,1)
cif_out.lower <- round(100*(cif_timepoints$est-1.96*sqrt(cif_timepoints$var)),1)
cif_out.upper <- round(100*(cif_timepoints$est+1.96*sqrt(cif_timepoints$var)),1)

# Significance test

aha<-tidycmprsk::cuminc(Surv(fumos_O3.yr, surv_type) ~ Risk.group, dat_filtered)
aha
# Table

tidycmprsk::cuminc(Surv(fumos_O3.yr, surv_type) ~ Risk.group, dat_filtered)  %>%
  tbl_cuminc(outcomes= c("BC death", "OC death"), 
             times = c(2,5,8),
             label_header = "**Year {time}**",
             estimate_fun = gtsummary::label_style_number(scale = 100, digits = 1)) %>%
  add_nevent() %>%
  add_n() 


