####################################################################################################

## Project:        CADSET 2.0 - COPD epigenetics across the lifespan
## File name:      run_robust_regressions.R
## Description:    script to test differentially methylated position using robust regression
## Authors:        Julieta Viglino
##                 Sandra Casas
## Version:    	   2023-11
## R version used: R 4.2.2 (for the design and testing of this script)

####################################################################################################


# IMPORTANT: THIS SCRIPT SHOULD ONLY BE IF YOUR COHORT INCLUDES ADOLESCENTS/ADULTS
## IF ITS A ONLY CHILDREN COHORT, DONT NEED TO RUN THIS SCRIPT

#### *********************** ####
#### LOAD DATA AND LIBRARIES ####
#### *********************** ####

# if any of the following libraries are not installed, please install them through "install.packages()" 
suppressMessages({
  library(dplyr)
  library(MASS)
  library(parallel)
  library(sandwich)
  library(lmtest)
})

# NOTE 1: if during your cohort-specific methylation analysis batch was corrected while performing the models
# (i.e using negative control PCs or batch as model covariate), this should also be included in the formula!
# Please add your specific batch variables to models formulas is applicable

# NOTE 2: Modify the cell populations in the formulas if you have different ones


# load data:
# load clinical_data that was already created in the script create_clinical_df.R
clinical_df <-  read.csv("../data/clinical_df.csv") %>%
   # keep only the adolescents/adults from the cohort [if you only have children in your cohort dont need to run this script!!! ]
    subset(age > 12) %>%
    mutate(sex = as.factor(sex),
           smoking_status = factor(smoking_status, levels = c("NS", "FS", "CS")), # modify accordingly if there are less ctageories in your cohort
           age_category = as.factor(age_category),
           pack_years = as.numeric(pack_years))

# load methylation data
# note: user should complete with the path were you methylation data is. note that probes should be rows and patients columns
methyl_data <- readRDS("../data/methylation_matrix.RDS")
methyl_data <- as.matrix(methyl_data)
dim(methyl_data)

# order methylation matrix so that patients are in the same order than the clinical df and transpose methylation mat
methyl_data <- methyl_data[,as.character(clinical_df$ID)] %>% t()

# make sure the order of the IDs match
rownames(clinical_df) <- clinical_df$ID
table(rownames(clinical_df) == rownames(methyl_data))


#### **************************************************************************** ####
#### ROBUST REGRESSION MODEL WITH FEV AS RESPONSE: WITH BLOOD CELLS AS COVARIATES ####
#### **************************************************************************** ####

# set cores_ parameter for parallelization. parallel::detectCores() uses all cores available but can be modified by user
cores_ <- detectCores()

# split methylation data set into multiple smaller datasets 
methyl_s <- split.data.frame(t(methyl_data), 1:ncol(methyl_data) %% 2500) # you can change this number depending on your system computational resources

# empty lists to store results
cpgs <- list()
effect <- list()
std_error <- list() 
p_value <- list()
gc()

# set a progress bar to display in screen -> it will show the progress of models estimation
pb <- txtProgressBar(min = 0, max = length(methyl_s), style = 3)


start.time <- Sys.time()

for(idx_split in seq_along(methyl_s)){
  
  methyl_s_t <- methyl_s[[idx_split]] %>% t()
  df <- bind_cols(clinical_df, methyl_s_t)
  cpgs_name <- colnames(methyl_s_t)
  
  # construct formulas -> modify formula to include PCs/Batch effects or change cell populations if appropiate!!!
  formulas_ <- mclapply(cpgs_name, 
                       function(i) formula(paste("fevml ~", i, " + age + age_squared + height_meters + height_meters_squared + sex + smoking_status + pack_years + CD8T + CD4T + NK + Bcell + Neu")),
                       mc.cores = cores_)
  
  # perform rlm and estimate se and pvalues
  rlm_fit <- mclapply(formulas_, function(f) rlm(f, data = df), mc.cores = cores_)
  coef_test <- mclapply(rlm_fit, function(fit_) coeftest(fit_, vcoc = vcocHC(fit_)), mc.cores = cores_)
  
  # get estimates 
  cpgs.tmp <- mclapply(coef_test, function(i) rownames(i)[2] , mc.cores = cores_) %>% unlist()
  effect.tmp <-  mclapply(coef_test, function(i) i[2, "Estimate"] , mc.cores = cores_) %>% unlist()
  std_error.tmp <-  mclapply(coef_test, function(i) i[2, "Std. Error"] , mc.cores = cores_) %>% unlist()
  p_value.tmp <-   mclapply(coef_test, function(i) i[2, "Pr(>|z|)"] , mc.cores = cores_) %>% unlist()
  
  # append estimates of current model to results list
  cpgs <- append(cpgs, cpgs.tmp)
  effect <- append(effect, effect.tmp)
  std_error <- append(std_error, std_error.tmp)
  p_value <- append(p_value, p_value.tmp)

  setTxtProgressBar(pb, idx_split)
  
}

end.time <- Sys.time()
(time.taken <- end.time - start.time)


# combine results into a data frame and save it as csv
res <- data.frame(probe = unlist(cpgs), 
                 effect_rlm = unlist(effect), 
                 std_error_rlm = unlist(std_error),
                 p_value_rlm = unlist(p_value))

write.csv(res, file = "../results/rlm/rlm_results_fevml_with_cell_pop.csv", row.names = F)



#### ******************************************************************************* ####
#### ROBUST REGRESSION MODEL WITH FEV AS RESPONSE: WITHOUT BLOOD CELLS AS COVARIATES ####
#### ******************************************************************************* ####

# remove all r objects genereated during previous model estimation
# we do want to keep the main datasets
objects_to_keep <- c("clinical_df", "methyl_data", "methyl_s", "cores_")
rm(list = setdiff(ls(), objects_to_keep))

# follow same steps than above to obatin the model

# empty lists to store results
cpgs <- list()
effect <- list()
std_error <- list() 
p_value <- list()
gc()

# set a progress bar to display in screen -> it will show the progress of models estimation
pb <- txtProgressBar(min = 0, max = length(methyl_s), style = 3)


start.time <- Sys.time()

for(idx_split in seq_along(methyl_s)){
  
  methyl_s_t <- methyl_s[[idx_split]] %>% t()
  df <- bind_cols(clinical_df, methyl_s_t)
  cpgs_name <- colnames(methyl_s_t)
  
  # construct formulas -> modify formula to include PCs/Batch effects!!!
  formulas_ <- mclapply(cpgs_name, 
                       function(i) formula(paste("fevml ~", i, " + age + age_squared + height_meters + height_meters_squared + sex + smoking_status + pack_years")),
                       mc.cores = cores_)
  
  # perform rlm and estimate se and pvalues
  rlm_fit <- mclapply(formulas_, function(f) rlm(f, data = df), mc.cores = cores_)
  coef_test <- mclapply(rlm_fit, function(fit_) coeftest(fit_, vcoc = vcocHC(fit_)), mc.cores = cores_)
  
  # get estimates 
  cpgs.tmp <- mclapply(coef_test, function(i) rownames(i)[2] , mc.cores = cores_) %>% unlist()
  effect.tmp <-  mclapply(coef_test, function(i) i[2, "Estimate"] , mc.cores = cores_) %>% unlist() 
  std_error.tmp <-  mclapply(coef_test, function(i) i[2, "Std. Error"] , mc.cores = cores_) %>% unlist()
  p_value.tmp <- mclapply(coef_test, function(i) i[2, "Pr(>|z|)"] , mc.cores = cores_) %>% unlist()
  
  # append estimates of current model to results list
  cpgs <- append(cpgs, cpgs.tmp)
  effect <- append(effect, effect.tmp)
  std_error <- append(std_error, std_error.tmp)
  p_value <- append(p_value, p_value.tmp)
  
  setTxtProgressBar(pb, idx_split)
  
}

end.time <- Sys.time()
(time.taken <- end.time - start.time)


# combine results into a data frame and save it as csv
res <- data.frame(probe = unlist(cpgs), 
                 effect_rlm = unlist(effect), 
                 std_error_rlm = unlist(std_error),
                 p_value_rlm = unlist(p_value))

write.csv(res, file = "../results/rlm/rlm_results_fevml_without_cell_pop.csv", row.names = F)

