## Causal Machine Learning Course Project by Yu Nu (yn292@cornell.edu)
## R codes for data replication purpose 

## Aggregate ATT via staggered synthetic DiD ##

# Install R package by Arkhangelsky et al. (2021) 
devtools::install_github("synth-inference/synthdid")   
library(synthdid)  
library(data.table)
library(plm)

# Load and prepare data 
Mdat = fread("Winnow_waste_data.csv")
Mdat <- Mdat[, logweight := log(weight_g+1)]
Mdat <- Mdat[, logwasteticket := log(waste_per_ticket+1)]
Mdat[, time_to_treat := period - treat_time_var]  
Mdat <- Mdat[time_to_treat <= 180]  # across half a year post adoption

site_list <- unique(Mdat$target_site)
result_matrix_sdid = matrix(ncol = 6, nrow = length(site_list))
result_matrix_sc   = matrix(ncol = 6, nrow = length(site_list))
result_matrix_did  = matrix(ncol = 6, nrow = length(site_list))


for (i in 1:length(site_list)){
  print(i)
  Mdat_s <- subset(Mdat, (target_site == site_list[i])) 
  Mdat_s <- Mdat_s[, c("site_label", "period", "treated_or_not", "logweight"), with = FALSE]
  Mdat_s_balanced <- make.pbalanced(Mdat_s, balance.type = "fill", index = c('site_label', 'period') )
  
  data <- Mdat_s_balanced
  untreated = 0
  unit = 1
  time_var = 2
  treatment_var = 3
  outcome_var = 4
  
  # Creating initial treatment time variable
  state_initial_treat = matrix(nrow=length(unique(data[, unit])), ncol=2)
  state_initial_treat[, 1] = unique(data[, unit])
  
  for (j in 1:nrow(state_initial_treat)){
    state_initial_treat[j, 2] = min(data[which(data[, unit] == state_initial_treat[j, 1] & data[, treatment_var] == 1), time_var])
    state_initial_treat[j, 2] = ifelse(state_initial_treat[j, 2] == "Inf", 0, state_initial_treat[j, 2])
  }
  
  state_initial_treat = as.data.frame(state_initial_treat)
  colnames(state_initial_treat)[1] = "unit"
  colnames(state_initial_treat)[2] = "initial_treat_period"
  colnames(data)[unit] = "unit"
  colnames(data)[time_var] = "period" 
  col_name_vec = colnames(data)
  
  data = merge(data, state_initial_treat, by = "unit")
  time_var = as.numeric(which(colnames(data) == col_name_vec[time_var]))
  outcome_var = as.numeric(which(colnames(data) == col_name_vec[outcome_var]))
  treatment_var = as.numeric(which(colnames(data) == col_name_vec[treatment_var]))
  data$initial_treat_var = as.numeric(data$initial_treat_period)
  initial_treat_var = as.numeric(which(colnames(data) == "initial_treat_var"))
  data[, treatment_var] = ifelse(data[, treatment_var] == untreated, 0, data[, treatment_var])
  
  treated_periods = unique(data[which(data[, initial_treat_var] > 0), initial_treat_var])
  result_matrix_sdid[i, 1] = treated_periods 
  result_matrix_sc[i, 1]   = treated_periods 
  result_matrix_did[i, 1]  = treated_periods
  
  data = na.omit(data)
  data <- make.pbalanced(data, balance.type = "shared.times", index = c("unit", "period") )
  
  # Synthetic DiD 
  subbed = data[which(data[, initial_treat_var] == treated_periods|data[, initial_treat_var] == untreated), ]
  colnames(subbed)[ncol(subbed)] = "post_treat"
  subbed[, ncol(subbed)+1] = 0
  subbed$post_treat = ifelse(subbed[, time_var] >= subbed[, initial_treat_var] & subbed[, initial_treat_var]!=0, 1, 0)
  subbed_test <- subset(subbed,(initial_treat_period==treated_periods))
  
  if (min(subbed_test$treated_or_not)==0){
    
    setup = synthdid::panel.matrices(subbed, unit=unit, time=time_var, 
                                     outcome = outcome_var, treatment = "post_treat")
    if (setup$T0 >= 2) { 
      tau_hat   = synthdid_estimate(setup$Y, setup$N0, setup$T0)
      tau_sc    = sc_estimate(setup$Y, setup$N0, setup$T0)
      tau_did   = did_estimate(setup$Y, setup$N0, setup$T0)
      estimates = list(tau_did, tau_sc, tau_hat)
      names(estimates) = c('Diff-in-Diff', 'Synthetic Control', 'Synthetic Diff-in-Diff')
    }
    
    result_matrix_sdid[i, 2] = as.numeric(tau_hat)
    result_matrix_sdid[i, 3] = sqrt(vcov(tau_hat, method = "placebo"))
    result_matrix_sdid[i, 4] = result_matrix_sdid[i, 2] - (1.96*result_matrix_sdid[i, 3])
    result_matrix_sdid[i, 5] = result_matrix_sdid[i, 2] + (1.96*result_matrix_sdid[i, 3])
    result_matrix_sdid[i, 6] = nrow(subbed[which(subbed$post_treat == 1), ])
    
    result_matrix_sc[i, 2]   = as.numeric(tau_sc)
    result_matrix_sc[i, 3]   = sqrt(vcov(tau_sc, method = "placebo"))
    result_matrix_sc[i, 4]   = result_matrix_sc[i, 2] - (1.96*result_matrix_sc[i, 3])
    result_matrix_sc[i, 5]   = result_matrix_sc[i, 2] + (1.96*result_matrix_sc[i, 3])
    result_matrix_sc[i, 6]   = nrow(subbed[which(subbed$post_treat == 1), ])
    
    result_matrix_did[i, 2]  = as.numeric(tau_did)
    result_matrix_did[i, 3]  = sqrt(vcov(tau_did, method = "placebo"))
    result_matrix_did[i, 4]  = result_matrix_did[i, 2] - (1.96*result_matrix_did[i, 3])
    result_matrix_did[i, 5]  = result_matrix_did[i, 2] + (1.96*result_matrix_did[i, 3])
    result_matrix_did[i, 6]  = nrow(subbed[which(subbed$post_treat == 1), ])
    
  }}

# Output result matrix via synthetic DiD and calculate ATT
result_matrix_sdid = as.data.frame(result_matrix_sdid)
colnames(result_matrix_sdid)[1] = "treatment_period"
colnames(result_matrix_sdid)[2] = "tau_hat"
colnames(result_matrix_sdid)[3] = "standard_error"
colnames(result_matrix_sdid)[4] = "lower_95_CI"
colnames(result_matrix_sdid)[5] = "upper_95_CI"
result_matrix_sdid$weight = result_matrix_sdid[, 6]/sum(result_matrix_sdid[, 6])
result_matrix_sdid[, 6] = NULL

overall_estimator_sdid = sum(result_matrix_sdid[, 6] * result_matrix_sdid[, 2])
B = result_matrix_sdid[, 2]
x = data[, 3]
y = data[, 4]
x = as.matrix(x)
y = as.matrix(y)
B = as.matrix(B)
influence = matrix(nrow=nrow(x), ncol=nrow(B))

# Empirical equivalents for influence functions
for(r in 1:nrow(B)){
  for(q in 1:nrow(x)){
    influence[q, r] = ((x[q] - mean(x))/var(x))*((y[q] - mean(y))*B[r]*(x[q] - mean(x)))}
}

V = rnorm(nrow(data), 0, 1)
variance = (1/((nrow(x)^2)))*(t(influence) %*% influence)
weights_sdid = result_matrix_sdid$weight
variance_sdid = t(weights_sdid) %*% variance %*% weights_sdid
overall_se_sdid = sqrt(sqrt(variance_sdid^2))        
overall_t_sdid = (overall_estimator_sdid/(overall_se_sdid/sqrt(nrow(data))))
overall_p_sdid = 2*pt(q = sqrt(((overall_estimator_sdid/overall_se_sdid)^2)), 
                      df = (nrow(data)-length(unique(data[, time_var]))-length(unique(data[, unit]))), 
                      lower.tail=FALSE)
lower_95_CI_sdid = overall_estimator_sdid - (1.96*overall_se_sdid)
upper_95_CI_sdid = overall_estimator_sdid + (1.96*overall_se_sdid)
print(paste0("ATT= ", overall_estimator_sdid, 
             ifelse(overall_p_sdid < .05 & overall_p_sdid > .01,"**", 
                    ifelse(overall_p_sdid < .1 & overall_p_sdid > .05,"*",
                           ifelse(overall_p_sdid < .01,"***","")))))
print(paste0("SE= ", overall_se_sdid))
print(paste0("p= ", overall_p_sdid))
print(paste0("95% CI: (", lower_95_CI_sdid, ", ", upper_95_CI_sdid,")"))
print("*p<0.1; **p<0.5, ***p<0.01")


# Output result matrix via synthetic control and calculate ATT
result_matrix_sc = as.data.frame(result_matrix_sc)
colnames(result_matrix_sc)[1] = "treatment_period"
colnames(result_matrix_sc)[2] = "tau_sc"
colnames(result_matrix_sc)[3] = "standard_error"
colnames(result_matrix_sc)[4] = "lower_95_CI"
colnames(result_matrix_sc)[5] = "upper_95_CI"
result_matrix_sc$weight = result_matrix_sc[, 6]/sum(result_matrix_sc[, 6])
result_matrix_sc[, 6] = NULL

overall_estimator_sc = sum(result_matrix_sc[, 6] * result_matrix_sc[, 2])
B_sc = result_matrix_sc[, 2]
x = data[, 3]
y = data[, 4]
x = as.matrix(x)
y = as.matrix(y)
B_sc = as.matrix(B_sc)
influence_sc = matrix(nrow=nrow(x), ncol=nrow(B_sc))

for(r in 1:nrow(B_sc)){
  for(q in 1:nrow(x)){
    influence_sc[q, r] = ((x[q] - mean(x))/var(x))*((y[q] - mean(y))*B_sc[r]*(x[q] - mean(x)))}
}

V = rnorm(nrow(data), 0, 1)
variance_sc = (1/((nrow(x)^2)))*(t(influence_sc) %*% influence_sc)
weights_sc = result_matrix_sc$weight
variance_sc   = t(weights_sc) %*% variance_sc %*% weights_sc
overall_se_sc = sqrt(sqrt(variance_sc^2))       
overall_t_sc  = (overall_estimator_sc/(overall_se_sc/sqrt(nrow(data))))
overall_p_sc  = 2*pt(q=sqrt(((overall_estimator_sc/overall_se_sc)^2)), 
                     df=(nrow(data)-length(unique(data[, time_var]))-length(unique(data[, unit]))), 
                     lower.tail=FALSE)
lower_95_CI_sc = overall_estimator_sc - (1.96*overall_se_sc)
upper_95_CI_sc = overall_estimator_sc + (1.96*overall_se_sc)
print(paste0("ATT= ", overall_estimator_sc, 
             ifelse(overall_p_sc < .05 & overall_p_sc > .01,"**",
                    ifelse(overall_p_sc < .1 & overall_p_sc > .05,"*",
                           ifelse(overall_p_sc < .01,"***","")))))
print(paste0("SE= ", overall_se_sc))
print(paste0("p= ", overall_p_sc))
print(paste0("95% CI: (", lower_95_CI_sc, ", ", upper_95_CI_sc,")"))
print("*p<0.1; **p<0.5, ***p<0.01")


# Output result matrix via DiD and calculate ATT 
result_matrix_did = as.data.frame(result_matrix_did)
colnames(result_matrix_did)[1] = "treatment_period"
colnames(result_matrix_did)[2] = "tau_did"
colnames(result_matrix_did)[3] = "standard_error"
colnames(result_matrix_did)[4] = "lower_95_CI"
colnames(result_matrix_did)[5] = "upper_95_CI"
result_matrix_did$weight = result_matrix_did[, 6]/sum(result_matrix_did[, 6])
result_matrix_did[, 6] = NULL

overall_estimator_did = sum(result_matrix_did[, 6] * result_matrix_did[, 2])
B_did = result_matrix_did[, 2]
x = data[, 3]
y = data[, 4]
x = as.matrix(x)
y = as.matrix(y)
B_did = as.matrix(B_did)
influence_did = matrix(nrow=nrow(x), ncol=nrow(B_did))

for(r in 1:nrow(B_did)){
  for(q in 1:nrow(x)){
    influence_did[q, r] = ((x[q] - mean(x))/var(x))*((y[q] - mean(y))*B_did[r]*(x[q] - mean(x)))}
}

V = rnorm(nrow(data), 0, 1)
variance_did = (1/((nrow(x)^2)))*(t(influence_did) %*% influence_did)
weights_did = result_matrix_did$weight
variance_did = t(weights_did) %*% variance_did %*% weights_did
overall_se_did = sqrt(sqrt(variance_did^2))        
overall_t_did = (overall_estimator_did/(overall_se_did/sqrt(nrow(data))))
overall_p_did = 2*pt(q = sqrt(((overall_estimator_did/overall_se_did)^2)), 
                     df = (nrow(data)-length(unique(data[, time_var]))-length(unique(data[, unit]))), 
                     lower.tail=FALSE)
lower_95_CI_did = overall_estimator_did - (1.96*overall_se_did)
upper_95_CI_did = overall_estimator_did + (1.96*overall_se_did)
print(paste0("ATT= ", overall_estimator_did, 
             ifelse(overall_p_did < .05 & overall_p_did > .01,"**", 
                    ifelse(overall_p_did < .1 & overall_p_did > .05,"*",
                           ifelse(overall_p_did < .01,"***","")))))
print(paste0("SE= ", overall_se_did))
print(paste0("p= ", overall_p_did))
print(paste0("95% CI: (", lower_95_CI_did, ", ", upper_95_CI_did,")"))
print("*p<0.1; **p<0.5, ***p<0.01")


# Plotting outcome var over time for a certain treated unit
synthdid_plot(estimates, se.method='placebo')


