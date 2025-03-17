####################################################################################################

## Project:        CADSET 2.0 - COPD epigenetics across the lifespan
## File name:      analyze_methylation_scores.R
## Description:    script to analyze methlation risk scores performance once obtained
## Authors:        Julieta Viglino
##                 Sandra Casas
## Version:    	   2023-11
## R version used: R 4.2.2 (for the design and testing of this script)

####################################################################################################


#### *********************************** ####
#### LOAD DATA & DEFINE USEFULL FUNCTION ####
#### *********************************** ####

rm(list = ls())


# load libraries
library(dplyr)
library(pROC)

# load calculated MRSs and clinical df that was created in **************
clinical_df <- read.csv("../data/clinical_df.csv")  %>%
  mutate(sex = factor(sex, levels = c("F", "M")),
         smoking_status = factor(smoking_status, levels = c("NS", "FS", "CS")), # modify accordingly if there are less ctageories in your cohort
         maternal_smoking = factor(maternal_smoking, levels = c(0, 1)),
         age_category = as.factor(age_category),
         pack_years = as.numeric(pack_years))

scores <- read.csv("../data/mrs_results.csv")

# left join clinical data and methylation scores
clinical_df <- clinical_df %>% left_join(scores, by = "ID")


# define vectors that will be used later to perform analysis
pulm_fn_vars <- c("fevml", "fvcml", "fev_fvc_ratio", "fev_perc")
calc_scores <- colnames(scores)[!"ID" == colnames(scores)]

#### ***************************************** ####
#### GET ASSOCIATION OF MRS WITH CLINICAL VARS ####
#### ***************************************** ####

# 1) obtain spearman correlation between calculated scores and clinical features
# note: modify lines below depending on your specific dataset (remove pack_years if your cohort is only children and change your cell populations if appropiate
cor_ <- cor(clinical_df[,c(calc_scores, pulm_fn_vars, "age", "age_squared", "pack_years",
                           "height_meters", "height_meters_squared",
                           "CD8T", "CD4T", "NK", "Bcell", "Neu")],
           use = "pairwise.complete.obs",
           method = "spearman")

write.csv(cor_, file = "../results/mrs/insights/scores_correlation_tot.csv", row.names = TRUE)
rm(cor_)

# 2) adjust a linear model of <lung function metric> ~ MRS and evaluate adjustment

# create function to perform linear regressions of methylation risk scores and pulmonary function
get_score_adj <- function(lung_function, score){
  lm.fit <- lm(lung_function ~ score)
  res <- list(
    "estimate" =  summary(lm.fit)$coefficients["score", "Estimate"],
    "p_value" = summary(lm.fit)$coefficients["score", "Pr(>|t|)"],
    "std_error" = summary(lm.fit)$coefficients["score", "Std. Error"],
    "min_conf" = confint(lm.fit, level = 0.95)["score", "2.5 %"],
    "max_conf" = confint(lm.fit, level = 0.95)["score", "97.5 %"],
    "r.squared" = summary(lm.fit)$r.squared,
    "adj.r.squared" = summary(lm.fit)$adj.r.squared)

  res <- data.frame(unlist(res))
  return(res)
  }

results.df <- data.frame(row.names = c("estimate", "p_value", "std_error", "min_conf", "max_conf", "r.squared", "adj.r.squared"))
for (i in pulm_fn_vars){ 
  for (k in calc_scores){
    tmp.name <- paste(i, k, sep = "_")
    res.tmp.df <- get_score_adj(lung_function = clinical_df[, i], score = clinical_df[, k]) %>%
                    rename(!!quo_name(tmp.name) :=  1)
    
    results.df <- bind_cols(results.df, res.tmp.df)
    rm(tmp.name, res.tmp.df)
  }
}

write.csv(results.df, file = "../results/mrs/insights/scores_lm_tot.csv", row.names = TRUE)
rm(results.df, i, k)

# 3) separate score in deciles and get mean and sd values of lung function metrics in each decil

# function to assign scores deciles and calculate the mean and sd of fevml, fev1%, fvcml and fev/fvc in each decil
get_lung_function_avg_per_score_decil <- function(score, fev_perc, fevml, fvcml, fev_fvc_ratio){
  dec_score <- ntile(score, 10) # compute deciles of score
  mean.dec <- data.frame(fevml = fevml,
                         fev_perc = fev_perc,
                         fev_fvc_ratio = fev_fvc_ratio,
                         fvcml = fvcml,
                         dec_score = dec_score) %>%
    group_by(dec_score) %>% 
    summarise(n = n(),
              mean_fevml = mean(fevml, na.rm = T),
              sd_fevml = sd(fevml, na.rm = T),
              mean_fev_perc = mean(fev_perc, na.rm = T),
              sd_fev_perc = sd(fev_perc, na.rm = T),
              mean_fev_fvc_ratio = mean(fev_fvc_ratio, na.rm = T),
              sd_fev_fvc_ratio = sd(fev_fvc_ratio, na.rm = T),
              mean_fvcml = mean(fvcml, na.rm = T),
              sd_fvcml = sd(fvcml, na.rm = T),
              missing_fev_fvc_ratio = sum(is.na(fev_fvc_ratio)),
              missing_fvcml = sum(is.na(fvcml)),
              missing_fev_perc = sum(is.na(fev_perc)))
  return(mean.dec)
}
pulm_fn_vars <- c("fevml", "fvcml", "fev_fvc_ratio", "fev_perc")
# generate data frame to export results
decils_results_df <- data.frame()
for (k in calc_scores){
  decils_results_tmp <- get_lung_function_avg_per_score_decil(
                            score = clinical_df[,k],
                            fev_perc = clinical_df[,"fev_perc"],
                            fevml = clinical_df[,"fevml"],
                            fvcml = clinical_df[,"fvcml"],
                            fev_fvc_ratio = clinical_df[,"fev_fvc_ratio"]
                                ) %>%
    mutate(score = k)

  decils_results_df <- bind_rows(decils_results_df, decils_results_tmp)
  rm(decils_results_tmp)
}

write.csv(decils_results_df, paste("../results/mrs/insights/stats_per_score_dec_tot.csv", sep = ""), row.names = FALSE)
rm(decils_results_df)



#### ******************************************************* ####
#### GET ASSOCIATION OF MRS WITH CLINICAL VARS PER AGE RANGE ####
#### ******************************************************* ####

# remove all r objects genereated previously: only keep definitions and things to keep
objects_to_keep <-  c("get_score_adj", "get_lung_function_avg_per_score_decil", "pulm_fn_vars", "calc_scores", "clinical_df")
rm(list = setdiff(ls(), objects_to_keep))

age_ranges <- unique(clinical_df$age_category)

# compute different associations between mrs and lung function for each age range of the cohort
# in each iteration, several associations are computed and results statistics exportes as csv
for (age_range in age_ranges){


  cli.df.age.tmp <-  clinical_df %>%
    subset(age_category == age_range)
  
  # if there are less than 30 subjects in the given age group, continue to the next iteration - n too small to get stats!
  if (nrow(cli.df.age.tmp) < 30) {
    next  
  }
  
  # CORRELATION
  # calculate spearman correlation between scores and clinical features
  # note: modify lines below depending on your specific dataset (remove pack_years if your cohort is only children and change your cell populations if appropiate
  cor_ <- cor(cli.df.age.tmp[,c(calc_scores, pulm_fn_vars, "age", "pack_years", "height_meters", "CD8T", "CD4T", "NK", "Bcell", "Neu")],
              use = "pairwise.complete.obs",
             method = "spearman")
  write.csv(cor_, file = paste("../results/mrs/insights/scores_correlation_", age_range, ".csv", sep = ""), row.names = TRUE)
  
  
  # LM ADJUSTMENT
  # adjust linear regression per age group
  results.df <- data.frame(row.names = c("estimate","p_value", "min_conf","max_conf","std_error","r.squared","adj.r.squared"))
  for (i in pulm_fn_vars){
    for (k in calc_scores){
      tmp.name <- paste(i, k, sep = "_")
      res.tmp.df <- get_score_adj(lung_function =  cli.df.age.tmp[, i], score = cli.df.age.tmp[, k]) %>%
        rename(!!quo_name(tmp.name) :=  1)
      
      results.df <- bind_cols(results.df, res.tmp.df)
    }
  }
  write.csv(results.df, paste("../results/mrs/insights/scores_lm_", age_range, ".csv", sep = ""), row.names = TRUE)
  rm(results.df)


  # SCORES DECILS ANALYSIS
  # separate score in deciles and get proportions lung function averages in each decil
  decils_results_df <- data.frame()
  for (k in calc_scores){
   decils_results_tmp <- get_lung_function_avg_per_score_decil(
                            score = cli.df.age.tmp[,k],
                            fev_perc = cli.df.age.tmp[,"fev_perc"],
                            fevml = cli.df.age.tmp[,"fevml"],
                            fvcml = cli.df.age.tmp[,"fvcml"],
                            fev_fvc_ratio = cli.df.age.tmp[,"fev_fvc_ratio"]
                                ) %>%
    mutate(score = k)

   decils_results_df <- bind_rows(decils_results_df, decils_results_tmp)
   rm(decils_results_tmp)
  }
  
  write.csv(decils_results_df, paste("../results/mrs/insights/stats_per_score_dec_", age_range, ".csv", sep = ""), row.names = FALSE)
  rm(decils_results_df)


}

#### ******************************************************* ####
#### ***** EVALUATE IF MRS CAN INCREASE COPD PREDICTION **** ####
#### ******************************************************* ####
# aim: computete different models (binomial glm) and evaluate if adding the MRS as predictor variable can improves COPD prediction
# models will be performed for whole cohort and per age bins

# This part should only be run if:
## 1. your cohort is a mix of copd patients and controls (if no copd patients are present in your cohort, dont run)
(table(clinical_df$copd))
## 2. your cohort includes adults! (if its only a children cohort, dont need to run)
(table(clinical_df$age_category))


# remove all r objects genereated previously: only keep definitions and things to keep
objects_to_keep <-  c("get_score_adj", "get_copds_and_fev_per_score_decil", "pulm_fn_vars", "calc_scores", "clinical_df", "age_ranges")
rm(list = setdiff(ls(), objects_to_keep))

# important: keep only subjects that dont have NAs in copd variable [patients with NAs in the other relevant variables were already filtered in script one]
# and exclude childs
clinical_df_mod <- clinical_df %>% filter(!is.na(copd)) %>% subset(age_category != c("coord_blood", "children"))

# base models to compute:
# formulas for adolescents/adults:
f_m0 <- formula("copd ~ smoking_status + pack_years")
f_m1 <- formula("copd ~ smoking_status + pack_years + sex + age")
# adding raw MRS:
f_m2 <- formula("copd ~ smoking_status + pack_years + sex + age + MRS_fev_pop")
f_m3 <- formula("copd ~ smoking_status + pack_years + sex + age + MRS_fev_clin")
f_m4 <- formula("copd ~ smoking_status + pack_years + sex + age + MRS_sev")
# adding std MRS:
f_m2b <- formula("copd ~ smoking_status + pack_years + sex + age + std_MRS_fev_pop")
f_m3b <- formula("copd ~ smoking_status + pack_years + sex + age + std_MRS_fev_clin")
f_m4b <- formula("copd ~ smoking_status + pack_years + sex + age + std_MRS_sev")

# list all the models formulas:
f_list <- list("f_m0" = f_m0, "f_m1" = f_m1, "f_m2" = f_m2, "f_m3" = f_m3, "f_m4" = f_m4,  "f_m2b" = f_m2b, "f_m3b" = f_m3b, "f_m4b" = f_m4b)

# function to compute GLM of COPDs vs Controls based on previous formulas and obtain AUC 
get_risk_mod <- function(df, formula) {
  mod <- glm(formula, df, family = "binomial")
  fit_roc <- roc(response = df[,"copd"],
                 predictor = mod$fitted.values,
                 plot = F, percent = F, ci = T)
  
  res_auc <- c(fit_roc$ci[1], fit_roc$ci[2], fit_roc$ci[3]) # auc lower ci limit, auc, auc upper ci limit
  return(res_auc)
}



# 1) run models for all cohort (without subsetting age bins)
auc_res_tmp <- data.frame(row.names = c("auc.low", "auc", "auc.up"))

clinical_df_mod_adults <- clinical_df_mod %>% subset(age > 12)

for (f_n in names(f_list)){ # iterate over models formulas
  auc_res_tmp[, f_n] <- get_risk_mod(df = clinical_df_mod_adults, formula = f_list[[f_n]]) # check in R terminal messages that control = 0, case = 1
}

write.csv(auc_res_tmp, paste("../results/mrs/insights/copd_models_auc_tot.csv", sep = ""), row.names = TRUE)
rm(auc_res_tmp)


# 2) run models per age bins, exlucidng coord_blood and babies
for (age_range in age_ranges) {
  cli.df.age.tmp <-  clinical_df_mod %>% subset(age_category == age_range)
  if ((nrow(cli.df.age.tmp) < 30) | (age_range %in% c("coord_blood", "children"))) {
    next  
  }
  # iterate over formulas and get auc and auc confidence intervals
  auc_res_tmp <- data.frame(row.names = c("auc.low", "auc", "auc.up"))
  for (f_n in names(f_list)){
    auc_res_tmp[, f_n] <- get_risk_mod(df = cli.df.age.tmp, formula = f_list[[f_n]])
  }

  write.csv(auc_res_tmp, paste("../results/mrs/insights/copd_models_auc_", age_range, ".csv", sep = ""), row.names = TRUE)

}




#### *********************************************************** ####
#### ***** EVALUATE IF MRS CAN INCREASE FEV1/FVC PREDICTION **** ####
#### *********************************************************** ####

# aim: evaluate if adding the MRS as predictor variable to a base model including age, sex and smoking status can improves FEV1/FVC prediction
# model will be performed for whole cohort and per age bins

# remove all r objects genereated previously: only keep definitions and things to keep
objects_to_keep <-  c("get_score_adj", "get_copds_and_fev_per_score_decil", "pulm_fn_vars", "calc_scores", "clinical_df", "age_ranges")
rm(list = setdiff(ls(), objects_to_keep))

# important: keep only subjects that dont have NAs IN fev_fvc_ratio variable [patients with NAs in the other relevant variables were already filtered in script one]
clinical_df_mod <- clinical_df %>% filter(!is.na(fev_fvc_ratio))

# base models to compute for adolescents/adults:
f_m0 <- formula("fev_fvc_ratio ~ smoking_status + pack_years")
f_m1 <- formula("fev_fvc_ratio ~ smoking_status + pack_years + sex + age")
f_m2 <- formula("fev_fvc_ratio ~ smoking_status + pack_years + sex + age + MRS_fev_pop")
f_m3 <- formula("fev_fvc_ratio ~ smoking_status + pack_years + sex + age + MRS_fev_clin")
f_m4 <- formula("fev_fvc_ratio ~ smoking_status + pack_years + sex + age + MRS_sev")
f_m2b <- formula("fev_fvc_ratio ~ smoking_status + pack_years + sex + age + std_MRS_fev_pop")
f_m3b <- formula("fev_fvc_ratio ~ smoking_status + pack_years + sex + age + std_MRS_fev_clin")
f_m4b <- formula("fev_fvc_ratio ~ smoking_status + pack_years + sex + age + std_MRS_sev")
# list all the models formulas:
f_list <- list("f_m0" = f_m0, "f_m1" = f_m1, "f_m2" = f_m2, "f_m3" = f_m3, "f_m4" = f_m4,  "f_m2b" = f_m2b, "f_m3b" = f_m3b, "f_m4b" = f_m4b)

# base models to compute for childrens:
f_m0_c <- formula("fev_fvc_ratio ~ maternal_smoking")
f_m1_c <- formula("fev_fvc_ratio ~ maternal_smoking + sex + age")
f_m2_c <- formula("fev_fvc_ratio ~ maternal_smoking + sex + age + MRS_fev_pop")
f_m3_c <- formula("fev_fvc_ratio ~ maternal_smoking + sex + age + MRS_fev_clin")
f_m4_c <- formula("fev_fvc_ratio ~ maternal_smoking + sex + age + MRS_sev")
f_m2b_c <- formula("fev_fvc_ratio ~ maternal_smoking + sex + age + std_MRS_fev_pop")
f_m3b_c <- formula("fev_fvc_ratio ~ maternal_smoking + sex + age + std_MRS_fev_clin")
f_m4b_c <- formula("fev_fvc_ratio ~ maternal_smoking + sex + age + std_MRS_sev")

# list all the models formulas:
f_list <- list("f_m0" = f_m0, "f_m1" = f_m1, "f_m2" = f_m2, "f_m3" = f_m3, "f_m4" = f_m4,  "f_m2b" = f_m2b, "f_m3b" = f_m3b, "f_m4b" = f_m4b)
f_list_c <- list("f_m0_c" = f_m0_c, "f_m1_c" = f_m1_c, "f_m2_c" = f_m2_c, "f_m3_c" = f_m3_c, "f_m4_c" = f_m4_c,  "f_m2b_c" = f_m2b_c, "f_m3b_c" = f_m3b_c, "f_m4b_c" = f_m4b_c)


# function to compute linear models and obtain prediction mean squared error (MSE)
get_mse_mod <- function(df, formula) {
  mod <- lm(formula, df)
  mse <- mean(mod$residuals^2)
  return(mse)
}

# 1) run models for all adolescents/adults in cohort
mse_res_tmp <- data.frame(row.names = c("MSE"))
clinical_df_mod_adults <- clinical_df_mod %>% subset(age > 12)
for (f_n in names(f_list)){ # iterate over models formulas
  mse_res_tmp[, f_n] <- get_mse_mod(df = clinical_df_mod_adults, formula = f_list[[f_n]]) # check in R terminal messages that control = 0, case = 1
}
write.csv(mse_res_tmp, paste("../results/mrs/insights/fev_fvc_models_mse_tot.csv", sep = ""), row.names = TRUE)
rm(mse_res_tmp)


# 2) run models per age bins (including chiuldren,  adolescents, and adults)
for (age_range in age_ranges) {
  cli.df.age.tmp <-  clinical_df_mod %>% subset(age_category == age_range)
  if (nrow(cli.df.age.tmp) < 30) {
    next
  }
  # iterate over formulas and get auc and auc confidence intervals
  mse_res_tmp <- data.frame(row.names = c("MSE"))

  # define forumlas that will be used based on age_range
  # if children, maternal smoking will be used !!
  if (age_range  %in% c("children", "coord_blood")) {
    for (f_n in names(f_list_c)){
      mse_res_tmp[, f_n] <- get_mse_mod(df = cli.df.age.tmp, formula = f_list_c[[f_n]])
    }
  }
  # if adults, self-smoking will be used
  else {
    for (f_n in names(f_list)){
      mse_res_tmp[, f_n] <- get_mse_mod(df = cli.df.age.tmp, formula = f_list[[f_n]])
    }
  }

  write.csv(mse_res_tmp, paste("../results/mrs/insights/fev_fvc_models_mse_", age_range, ".csv", sep = ""), row.names = TRUE)
  rm(mse_res_tmp)
}



#### *********************************************************** ####
#### ***** EVALUATE IF MRS CAN INCREASE FEV1% PREDICTION **** ####
#### *********************************************************** ####

# aim: evaluate if adding the MRS as predictor variable to a base model including age, sex and smoking status can improves FEV1 prediction
# model will be performed for whole cohort and per age bins

# remove all r objects genereated previously: only keep definitions and things to keep
objects_to_keep <-  c("get_score_adj", "get_copds_and_fev_per_score_decil", "pulm_fn_vars", "calc_scores", "clinical_df", "age_ranges")
rm(list = setdiff(ls(), objects_to_keep))

# important: keep only subjects that dont have NAs IN fev_fvc_ratio variable [patients with NAs in the other relevant variables were already filtered in script one]
clinical_df_mod <- clinical_df %>% filter(!is.na(fev_perc))

# base models to compute for adolescents/adults:
f_m0 <- formula("fev_perc ~ smoking_status + pack_years")
f_m1 <- formula("fev_perc ~ smoking_status + pack_years + sex + age")
f_m2 <- formula("fev_perc ~ smoking_status + pack_years + sex + age + MRS_fev_pop")
f_m3 <- formula("fev_perc ~ smoking_status + pack_years + sex + age + MRS_fev_clin")
f_m4 <- formula("fev_perc ~ smoking_status + pack_years + sex + age + MRS_sev")
f_m2b <- formula("fev_perc ~ smoking_status + pack_years + sex + age + std_MRS_fev_pop")
f_m3b <- formula("fev_perc ~ smoking_status + pack_years + sex + age + std_MRS_fev_clin")
f_m4b <- formula("fev_perc ~ smoking_status + pack_years + sex + age + std_MRS_sev")
# list all the models formulas:
f_list <- list("f_m0" = f_m0, "f_m1" = f_m1, "f_m2" = f_m2, "f_m3" = f_m3, "f_m4" = f_m4,  "f_m2b" = f_m2b, "f_m3b" = f_m3b, "f_m4b" = f_m4b)

# base models to compute for childrens:
f_m0_c <- formula("fev_perc ~ maternal_smoking")
f_m1_c <- formula("fev_perc ~ maternal_smoking + sex + age")
f_m2_c <- formula("fev_perc ~ maternal_smoking + sex + age + MRS_fev_pop")
f_m3_c <- formula("fev_perc ~ maternal_smoking + sex + age + MRS_fev_clin")
f_m4_c <- formula("fev_perc ~ maternal_smoking + sex + age + MRS_sev")
f_m2b_c <- formula("fev_perc ~ maternal_smoking + sex + age + std_MRS_fev_pop")
f_m3b_c <- formula("fev_perc ~ maternal_smoking + sex + age + std_MRS_fev_clin")
f_m4b_c <- formula("fev_perc ~ maternal_smoking + sex + age + std_MRS_sev")

# list all the models formulas:
f_list <- list("f_m0" = f_m0, "f_m1" = f_m1, "f_m2" = f_m2, "f_m3" = f_m3, "f_m4" = f_m4,  "f_m2b" = f_m2b, "f_m3b" = f_m3b, "f_m4b" = f_m4b)
f_list_c <- list("f_m0_c" = f_m0_c, "f_m1_c" = f_m1_c, "f_m2_c" = f_m2_c, "f_m3_c" = f_m3_c, "f_m4_c" = f_m4_c,  "f_m2b_c" = f_m2b_c, "f_m3b_c" = f_m3b_c, "f_m4b_c" = f_m4b_c)


# function to compute linear models and obtain prediction mean squared error (MSE)
get_mse_mod <- function(df, formula) {
  mod <- lm(formula, df)
  mse <- mean(mod$residuals^2)
  return(mse)
}

# 1) run models for all adolescents/adults in cohort
mse_res_tmp <- data.frame(row.names = c("MSE"))
clinical_df_mod_adults <- clinical_df_mod %>% subset(age > 12)
for (f_n in names(f_list)){ # iterate over models formulas
  mse_res_tmp[, f_n] <- get_mse_mod(df = clinical_df_mod_adults, formula = f_list[[f_n]]) # check in R terminal messages that control = 0, case = 1
}
write.csv(mse_res_tmp, paste("../results/mrs/insights/fev_perc_models_mse_tot.csv", sep = ""), row.names = TRUE)
rm(mse_res_tmp)


# 2) run models per age bins (including chiuldren,  adolescents, and adults)
for (age_range in age_ranges) {
  cli.df.age.tmp <-  clinical_df_mod %>% subset(age_category == age_range)
  if (nrow(cli.df.age.tmp) < 30) {
    next
  }
  # iterate over formulas and get auc and auc confidence intervals
  mse_res_tmp <- data.frame(row.names = c("MSE"))
  
  # define forumlas that will be used based on age_range
  # if children, maternal smoking will be used !!
  if (age_range  %in% c("children", "coord_blood")) {
    for (f_n in names(f_list_c)){
      mse_res_tmp[, f_n] <- get_mse_mod(df = cli.df.age.tmp, formula = f_list_c[[f_n]])
    }
  }
  # if adults, self-smoking will be used
  else {
    for (f_n in names(f_list)){
      mse_res_tmp[, f_n] <- get_mse_mod(df = cli.df.age.tmp, formula = f_list[[f_n]])
    }
  }
  
  write.csv(mse_res_tmp, paste("../results/mrs/insights/fev_perc_models_mse_", age_range, ".csv", sep = ""), row.names = TRUE)
  rm(mse_res_tmp)
}


#### *********************************************************** ####
#### ***** EVALUATE IF MRS CAN INCREASE FEV1ml PREDICTION **** ####
#### *********************************************************** ####

# aim: evaluate if adding the MRS as predictor variable to a base model including age, sex and smoking status can improves FEV1 prediction
# model will be performed for whole cohort and per age bins

# remove all r objects genereated previously: only keep definitions and things to keep
objects_to_keep <-  c("get_score_adj", "get_copds_and_fev_per_score_decil", "pulm_fn_vars", "calc_scores", "clinical_df", "age_ranges")
rm(list = setdiff(ls(), objects_to_keep))

# important: keep only subjects that dont have NAs IN fev_fvc_ratio variable [patients with NAs in the other relevant variables were already filtered in script one]
clinical_df_mod <- clinical_df %>% filter(!is.na(fevml))

# base models to compute for adolescents/adults:
f_m0 <- formula("fevml ~ smoking_status + pack_years")
f_m1 <- formula("fevml ~ smoking_status + pack_years + sex + age")
f_m2 <- formula("fevml ~ smoking_status + pack_years + sex + age + MRS_fev_pop")
f_m3 <- formula("fevml ~ smoking_status + pack_years + sex + age + MRS_fev_clin")
f_m4 <- formula("fevml ~ smoking_status + pack_years + sex + age + MRS_sev")
f_m2b <- formula("fevml ~ smoking_status + pack_years + sex + age + std_MRS_fev_pop")
f_m3b <- formula("fevml ~ smoking_status + pack_years + sex + age + std_MRS_fev_clin")
f_m4b <- formula("fevml ~ smoking_status + pack_years + sex + age + std_MRS_sev")
# list all the models formulas:
f_list <- list("f_m0" = f_m0, "f_m1" = f_m1, "f_m2" = f_m2, "f_m3" = f_m3, "f_m4" = f_m4,  "f_m2b" = f_m2b, "f_m3b" = f_m3b, "f_m4b" = f_m4b)

# base models to compute for childrens:
f_m0_c <- formula("fevml ~ maternal_smoking")
f_m1_c <- formula("fevml ~ maternal_smoking + sex + age")
f_m2_c <- formula("fevml ~ maternal_smoking + sex + age + MRS_fev_pop")
f_m3_c <- formula("fevml ~ maternal_smoking + sex + age + MRS_fev_clin")
f_m4_c <- formula("fevml ~ maternal_smoking + sex + age + MRS_sev")
f_m2b_c <- formula("fevml ~ maternal_smoking + sex + age + std_MRS_fev_pop")
f_m3b_c <- formula("fevml ~ maternal_smoking + sex + age + std_MRS_fev_clin")
f_m4b_c <- formula("fevml ~ maternal_smoking + sex + age + std_MRS_sev")

# list all the models formulas:
f_list <- list("f_m0" = f_m0, "f_m1" = f_m1, "f_m2" = f_m2, "f_m3" = f_m3, "f_m4" = f_m4,  "f_m2b" = f_m2b, "f_m3b" = f_m3b, "f_m4b" = f_m4b)
f_list_c <- list("f_m0_c" = f_m0_c, "f_m1_c" = f_m1_c, "f_m2_c" = f_m2_c, "f_m3_c" = f_m3_c, "f_m4_c" = f_m4_c,  "f_m2b_c" = f_m2b_c, "f_m3b_c" = f_m3b_c, "f_m4b_c" = f_m4b_c)


# function to compute linear models and obtain prediction mean squared error (MSE)
get_mse_mod <- function(df, formula) {
  mod <- lm(formula, df)
  mse <- mean(mod$residuals^2)
  return(mse)
}

# 1) run models for all adolescents/adults in cohort
mse_res_tmp <- data.frame(row.names = c("MSE"))
clinical_df_mod_adults <- clinical_df_mod %>% subset(age > 12)
for (f_n in names(f_list)){ # iterate over models formulas
  mse_res_tmp[, f_n] <- get_mse_mod(df = clinical_df_mod_adults, formula = f_list[[f_n]]) # check in R terminal messages that control = 0, case = 1
}
write.csv(mse_res_tmp, paste("../results/mrs/insights/fevml_models_mse_tot.csv", sep = ""), row.names = TRUE)
rm(mse_res_tmp)


# 2) run models per age bins (including chiuldren,  adolescents, and adults)
for (age_range in age_ranges) {
  cli.df.age.tmp <-  clinical_df_mod %>% subset(age_category == age_range)
  if (nrow(cli.df.age.tmp) < 30) {
    next
  }
  # iterate over formulas and get auc and auc confidence intervals
  mse_res_tmp <- data.frame(row.names = c("MSE"))
  
  # define forumlas that will be used based on age_range
  # if children, maternal smoking will be used !!
  if (age_range  %in% c("children", "coord_blood")) {
    for (f_n in names(f_list_c)){
      mse_res_tmp[, f_n] <- get_mse_mod(df = cli.df.age.tmp, formula = f_list_c[[f_n]])
    }
  }
  # if adults, self-smoking will be used
  else {
    for (f_n in names(f_list)){
      mse_res_tmp[, f_n] <- get_mse_mod(df = cli.df.age.tmp, formula = f_list[[f_n]])
    }
  }
  
  write.csv(mse_res_tmp, paste("../results/mrs/insights/fevml_models_mse_", age_range, ".csv", sep = ""), row.names = TRUE)
  rm(mse_res_tmp)
}



#### *********************************************************** ####
#### ***** EVALUATE IF MRS CAN INCREASE FEVCml PREDICTION **** ####
#### *********************************************************** ####

# aim: evaluate if adding the MRS as predictor variable to a base model including age, sex and smoking status can improves FVC prediction
# model will be performed for whole cohort and per age bins

# remove all r objects genereated previously: only keep definitions and things to keep
objects_to_keep <-  c("get_score_adj", "get_copds_and_fev_per_score_decil", "pulm_fn_vars", "calc_scores", "clinical_df", "age_ranges")
rm(list = setdiff(ls(), objects_to_keep))

# important: keep only subjects that dont have NAs IN fev_fvc_ratio variable [patients with NAs in the other relevant variables were already filtered in script one]
clinical_df_mod <- clinical_df %>% filter(!is.na(fvcml))

# base models to compute for adolescents/adults:
f_m0 <- formula("fvcml ~ smoking_status + pack_years")
f_m1 <- formula("fvcml ~ smoking_status + pack_years + sex + age")
f_m2 <- formula("fvcml ~ smoking_status + pack_years + sex + age + MRS_fev_pop")
f_m3 <- formula("fvcml ~ smoking_status + pack_years + sex + age + MRS_fev_clin")
f_m4 <- formula("fvcml ~ smoking_status + pack_years + sex + age + MRS_sev")
f_m2b <- formula("fvcml ~ smoking_status + pack_years + sex + age + std_MRS_fev_pop")
f_m3b <- formula("fvcml ~ smoking_status + pack_years + sex + age + std_MRS_fev_clin")
f_m4b <- formula("fvcml ~ smoking_status + pack_years + sex + age + std_MRS_sev")
# list all the models formulas:
f_list <- list("f_m0" = f_m0, "f_m1" = f_m1, "f_m2" = f_m2, "f_m3" = f_m3, "f_m4" = f_m4,  "f_m2b" = f_m2b, "f_m3b" = f_m3b, "f_m4b" = f_m4b)

# base models to compute for childrens:
f_m0_c <- formula("fvcml ~ maternal_smoking")
f_m1_c <- formula("fvcml ~ maternal_smoking + sex + age")
f_m2_c <- formula("fvcml ~ maternal_smoking + sex + age + MRS_fev_pop")
f_m3_c <- formula("fvcml ~ maternal_smoking + sex + age + MRS_fev_clin")
f_m4_c <- formula("fvcml ~ maternal_smoking + sex + age + MRS_sev")
f_m2b_c <- formula("fvcml ~ maternal_smoking + sex + age + std_MRS_fev_pop")
f_m3b_c <- formula("fvcml ~ maternal_smoking + sex + age + std_MRS_fev_clin")
f_m4b_c <- formula("fvcml ~ maternal_smoking + sex + age + std_MRS_sev")

# list all the models formulas:
f_list <- list("f_m0" = f_m0, "f_m1" = f_m1, "f_m2" = f_m2, "f_m3" = f_m3, "f_m4" = f_m4,  "f_m2b" = f_m2b, "f_m3b" = f_m3b, "f_m4b" = f_m4b)
f_list_c <- list("f_m0_c" = f_m0_c, "f_m1_c" = f_m1_c, "f_m2_c" = f_m2_c, "f_m3_c" = f_m3_c, "f_m4_c" = f_m4_c,  "f_m2b_c" = f_m2b_c, "f_m3b_c" = f_m3b_c, "f_m4b_c" = f_m4b_c)


# function to compute linear models and obtain prediction mean squared error (MSE)
get_mse_mod <- function(df, formula) {
  mod <- lm(formula, df)
  mse <- mean(mod$residuals^2)
  return(mse)
}

# 1) run models for all adolescents/adults in cohort
mse_res_tmp <- data.frame(row.names = c("MSE"))
clinical_df_mod_adults <- clinical_df_mod %>% subset(age > 12)
for (f_n in names(f_list)){ # iterate over models formulas
  mse_res_tmp[, f_n] <- get_mse_mod(df = clinical_df_mod_adults, formula = f_list[[f_n]]) # check in R terminal messages that control = 0, case = 1
}
write.csv(mse_res_tmp, paste("../results/mrs/insights/fvcml_models_mse_tot.csv", sep = ""), row.names = TRUE)
rm(mse_res_tmp)


# 2) run models per age bins (including chiuldren,  adolescents, and adults)
for (age_range in age_ranges) {
  cli.df.age.tmp <-  clinical_df_mod %>% subset(age_category == age_range)
  if (nrow(cli.df.age.tmp) < 30) {
    next
  }
  # iterate over formulas and get auc and auc confidence intervals
  mse_res_tmp <- data.frame(row.names = c("MSE"))
  
  # define forumlas that will be used based on age_range
  # if children, maternal smoking will be used !!
  if (age_range  %in% c("children", "coord_blood")) {
    for (f_n in names(f_list)){
      mse_res_tmp[, f_n] <- get_mse_mod(df = cli.df.age.tmp, formula = f_list_c[[f_n]])
    }
  }
  # if adults, self-smoking will be used
  else {
    for (f_n in names(f_list)){
      mse_res_tmp[, f_n] <- get_mse_mod(df = cli.df.age.tmp, formula = f_list[[f_n]])
    }
  }
  
  write.csv(mse_res_tmp, paste("../results/mrs/insights/fvcml_models_mse_", age_range, ".csv", sep = ""), row.names = TRUE)
  rm(mse_res_tmp)
}





#### ******************************************************* ####
#### ***** EVALUATE IF MRS CAN INCREASE SEVERITY PREDICTION **** ####
#### ******************************************************* ####
# aim: computete different models (binomial glm) and evaluate if adding the MRS as predictor variable can improves SEVERITY prediction
# models will be performed for whole cohort and per age bins
clinical_df_sev = clinical_df %>%
  subset(age > 12) %>%
  mutate(severity = case_when(copd == 0 ~ NA,
                              copd == 1 & fev_perc>= 50 ~ '0',
                              copd == 1 & fev_perc< 50 ~ '1')) %>%
  mutate(severity = as.factor(severity))


# This part should only be run if:
## 1. your cohort is a mix of copd patients and controls (if no copd patients are present in your cohort, dont run)
(table(clinical_df_sev$severity))
## 2. your cohort includes adults! (if its only a children cohort, dont need to run)
(table(clinical_df_sev$age_category))


# remove all r objects genereated previously: only keep definitions and things to keep
objects_to_keep <-  c("get_score_adj", "get_copds_and_fev_per_score_decil", "pulm_fn_vars", "calc_scores", "clinical_df_sev", "age_ranges")
rm(list = setdiff(ls(), objects_to_keep))

# important: keep only subjects that dont have NAs in copd variable [patients with NAs in the other relevant variables were already filtered in script one]
# and exclude childs
clinical_df_mod <- clinical_df_sev %>% filter(!is.na(severity)) %>% subset(age_category != c("coord_blood", "children"))

# base models to compute:
# formulas for adolescents/adults:
f_m0 <- formula("severity ~ smoking_status + pack_years")
f_m1 <- formula("severity ~ smoking_status + pack_years + sex + age")
# adding raw MRS:
f_m2 <- formula("severity ~ smoking_status + pack_years + sex + age + MRS_fev_pop")
f_m3 <- formula("severity ~ smoking_status + pack_years + sex + age + MRS_fev_clin")
f_m4 <- formula("severity ~ smoking_status + pack_years + sex + age + MRS_sev")
# adding std MRS:
f_m2b <- formula("severity ~ smoking_status + pack_years + sex + age + std_MRS_fev_pop")
f_m3b <- formula("severity ~ smoking_status + pack_years + sex + age + std_MRS_fev_clin")
f_m4b <- formula("severity ~ smoking_status + pack_years + sex + age + std_MRS_sev")

# list all the models formulas:
f_list <- list("f_m0" = f_m0, "f_m1" = f_m1, "f_m2" = f_m2, "f_m3" = f_m3, "f_m4" = f_m4,  "f_m2b" = f_m2b, "f_m3b" = f_m3b, "f_m4b" = f_m4b)

# function to compute GLM of COPDs vs Controls based on previous formulas and obtain AUC 
get_risk_mod <- function(df, formula) {
  mod <- glm(formula, df, family = "binomial")
  fit_roc <- roc(response = df[,"severity"],
                 predictor = mod$fitted.values,
                 plot = F, percent = F, ci = T)
  
  res_auc <- c(fit_roc$ci[1], fit_roc$ci[2], fit_roc$ci[3]) # auc lower ci limit, auc, auc upper ci limit
  return(res_auc)
}



# 1) run models for all cohort (without subsetting age bins)
auc_res_tmp <- data.frame(row.names = c("auc.low", "auc", "auc.up"))

clinical_df_mod_adults <- clinical_df_mod %>% subset(age > 12)

for (f_n in names(f_list)){ # iterate over models formulas
  auc_res_tmp[, f_n] <- get_risk_mod(df = clinical_df_mod_adults, formula = f_list[[f_n]]) # check in R terminal messages that control = 0, case = 1
}

write.csv(auc_res_tmp, paste("../results/mrs/insights/severity_models_auc_tot.csv", sep = ""), row.names = TRUE)
rm(auc_res_tmp)


# 2) run models per age bins, exlucidng coord_blood and babies
for (age_range in age_ranges) {
  cli.df.age.tmp <-  clinical_df_mod %>% subset(age_category == age_range)
  if ((nrow(cli.df.age.tmp) < 30) | (age_range %in% c("coord_blood", "children"))) {
    next  
  }
  # iterate over formulas and get auc and auc confidence intervals
  auc_res_tmp <- data.frame(row.names = c("auc.low", "auc", "auc.up"))
  for (f_n in names(f_list)){
    auc_res_tmp[, f_n] <- get_risk_mod(df = cli.df.age.tmp, formula = f_list[[f_n]])
  }
  
  write.csv(auc_res_tmp, paste("../results/mrs/insights/severity_models_auc_", age_range, ".csv", sep = ""), row.names = TRUE)
  
}




