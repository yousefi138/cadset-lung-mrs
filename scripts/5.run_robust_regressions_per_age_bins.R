####################################################################################################

## Project:        CADSET 2.0 - COPD epigenetics across the lifespan
## File name:      run_robust_regressions_per_age_bins.R
## Description:    script to test differentially methylated position per age group using robust regression
## Authors:        Sandra Casas
##                 Julieta Viglino
## Version:    	   2023-11
## R version used: R 4.2.2 (for the design and testing of this script)

####################################################################################################

rm(list = ls())
gc()

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
    mutate(sex = factor(sex, levels = c("F", "M")),
           smoking_status = factor(smoking_status, levels = c("NS", "FS", "CS")), # modify accordingly if there are less ctageories in your cohort
           maternal_smoking = factor(maternal_smoking, levels = c(0,1)), # comment if not applicable (if you dont have children in your cohort)
           age_category = as.factor(age_category),
           pack_years = as.numeric(pack_years))
rownames(clinical_df) <- clinical_df$ID

# load methylation data
# note: user should complete with the path were you methylation data is. note that probes should be rows and patients columns
methyl_data <- readRDS("../data/methylation_matrix.RDS")
methyl_data <- as.matrix(methyl_data)
# order methylation matrix so that patients are in the same order than the clinical df
methyl_data <- methyl_data[,as.character(clinical_df$ID)] %>% t()

# make sure the order of the IDs match
table(rownames(clinical_df) == rownames(methyl_data))


#### **************************************************************************** ####
#### ROBUST REGRESSION MODEL WITH FEV AS RESPONSE: WITH BLOOD CELLS AS COVARIATES ####
#### **************************************************************************** ####

# set cores_ parameter for parallelization. parallel::detectCores() uses all cores available but can be modified by user
cores_ <- parallel::detectCores()

# split methylation data set into multiple smaller datasets 
methyl_s <- split.data.frame(t(methyl_data), 1:ncol(methyl_data) %% 2000)


table(clinical_df$age_category)
age_ranges <- unique(clinical_df$age_category)

for (age_range in age_ranges){

  # empty lists to store results
  cpgs <- list()
  effect <- list()
  std_error <- list()
  p_value <- list()
  gc()
  

  start.time <- Sys.time()
  pb <- txtProgressBar(min = 0, max = length(methyl_s), style = 3) # progress bar associated to current age bin

  # subset data frame by age category of current iteration
  clinical_df_tmp <- clinical_df %>% subset(age_category == age_range)
  
  # if there are less than 30 patients in current age category, dont compute any model, n too small
  if (nrow(clinical_df_tmp) < 30) {
    next
  }

  for(idx_split in seq_along(methyl_s)){
    
    methyl_s_t <- methyl_s[[idx_split]] %>% t()
    methyl_s_t <- methyl_s_t[clinical_df_tmp$ID,]

    df <- bind_cols(clinical_df_tmp, methyl_s_t)
    cpgs_name <- colnames(methyl_s_t)
    
    # construct formulas
    # formulas depend on age group: for children/babies maternal_smoking is included
    # modify cell populations/include batch variables if appropiate
    if ((age_range == "coord_blood") | (age_range == "children")) {
      formulas_ <- mclapply(cpgs_name,
                         function(i) formula(paste("fevml ~", i, " + age + age_squared + height_meters + height_meters_squared + sex + maternal_smoking + CD8T + CD4T + NK + Bcell + Neu")),
                         mc.cores = cores_)
    }
    else{
      formulas_ <- mclapply(cpgs_name,
                   function(i) formula(paste("fevml ~", i, " + age + age_squared + height_meters + height_meters_squared + sex + smoking_status + pack_years + CD8T + CD4T + NK + Bcell + Neu")),
                   mc.cores = cores_)
    }

    # perform rlm and estimate se and pvalues
    rlm_fit <- mclapply(formulas_, function(f) rlm(f, data = df), mc.cores = cores_)
    coef_test <- mclapply(rlm_fit, function(fit_) coeftest(fit_, vcoc = vcocHC(fit_)), mc.cores = cores_)
    
    # get estimates 
    cpgs.tmp <- mclapply(coef_test, function(i) rownames(i)[2] , mc.cores = cores_) %>% unlist()
    effect.tmp <- mclapply(coef_test, function(i) i[2, "Estimate"] , mc.cores = cores_) %>% unlist()
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

  write.csv(res, file = paste("../results/rlm/rlm_results_fevml_with_cell_pop_", age_range, ".csv", sep = ""), row.names = F)
  gc()

}





#### ******************************************************************************* ####
#### ROBUST REGRESSION MODEL WITH FEV AS RESPONSE: WITHOUT BLOOD CELLS AS COVARIATES ####
#### ******************************************************************************* ####

# remove all r objects genereated during previous model estimation
# we do want to keep the main datasets
objects_to_keep <- c("clinical_df", "methyl_data", "methyl_s", "cores_", "age_ranges")
rm(list = setdiff(ls(), objects_to_keep))

# set cores_ parameter for parallelization. parallel::detectCores() uses all cores available but can be modified by user
cores_ <- parallel::detectCores()


for (age_range in age_ranges){

  # empty lists to store results
  cpgs <- list()
  effect <- list()
  std_error <- list()
  p_value <- list()
  gc()

  start.time <- Sys.time()
  pb <- txtProgressBar(min = 0, max = length(methyl_s), style = 3) # progress bar associated to current age bin

  clinical_df_tmp <- clinical_df %>% subset(age_category == age_range)
  if (nrow(clinical_df_tmp) < 30) {
    next
  }

  for(idx_split in seq_along(methyl_s)){
    
    methyl_s_t <- methyl_s[[idx_split]] %>% t()
    methyl_s_t <- methyl_s_t[as.character(clinical_df_tmp$ID),]

    df <- bind_cols(clinical_df_tmp, methyl_s_t)
    cpgs_name <- colnames(methyl_s_t)
    
    # construct formulas
    # formulas depend on age group: for children/babies maternal_smoking is included
    # include batch variables if appropiate
    if ((age_range == "coord_blood") | (age_range == "children")) {
      formulas_ <- mclapply(cpgs_name,
                         function(i) formula(paste("fevml ~", i, " + age + age_squared + height_meters + height_meters_squared + sex + maternal_smoking")),
                         mc.cores = cores_)
    }
    else {
      formulas_ <- mclapply(cpgs_name,
                   function(i) formula(paste("fevml ~", i, " + age + age_squared + height_meters + height_meters_squared + sex + smoking_status + pack_years")),
                   mc.cores = cores_)
    }

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

  
  write.csv(res, file = paste("../results/rlm/rlm_results_fevml_without_cell_pop_", age_range, ".csv", sep = ""), row.names = F)
  gc()

}


