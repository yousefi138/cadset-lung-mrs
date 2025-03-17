################################################################################################################

## Project:        CADSET 2.0 - COPD epigenetics across the lifespan
## File name:      create_clinical_df.R
## Description:    code to get summary statistics of the cohort clinical characteristics
## Authors:        Julieta Viglino
##                 Sandra Casas
## Version:    	   2023-11
## R version used: R 4.2.2 (for the design and testing of this script)

################################################################################################################

library(dplyr)
library(reshape2)

# load clinical_data that was already created in the script create_clinical_df.R
clinical_df <-  read.csv("../data/clinical_df.csv") %>%
    mutate(sex = as.factor(sex),
           copd = as.factor(copd),
           ID = as.character(ID),
           smoking_status = as.factor(smoking_status), # delete if you only have children in your cohort
           pack_years = as.numeric(pack_years),    # delete if you only have children in your cohort
           maternal_smoking = factor(maternal_smoking, levels = c(0,1)), # delete if not applicable (if you dont have children in your cohort)
           age_category = as.factor(age_category))     # delete if not applicable (if you dont have children in your cohort)



factor_vars <- colnames(clinical_df)[sapply(clinical_df, is.factor)]
numeric_vars <- colnames(clinical_df)[sapply(clinical_df, is.numeric)]


# _________________________________________________________________________________________________
# Get number of patients per age category 
n_per_age_cat <- table(clinical_df$age_category) %>% 
          as.data.frame() %>%
          dplyr::rename("age_category" = "Var1",
                        "n" = "Freq")

write.csv(n_per_age_cat, file = "../results/summary-stats/n_per_age_category.csv", row.names = T)




# _________________________________________________________________________________________________
# Get number of NAs per variable

nas_summary <-  colSums(is.na(clinical_df)) %>% as.data.frame()
write.csv(nas_summary, file = "../results/summary-stats/na_count.csv", row.names = T)


# _________________________________________________________________________________________________
# Including all cohort

summary_num_all <- clinical_df %>%
       dplyr::select(-age_category, -age_squared, -height_meters_squared) %>%
       melt(., id.vars  = c("ID")) %>% 
       subset(variable %in% numeric_vars) %>%
       mutate(variable = as.factor(as.character(variable)),
              value = as.numeric(value)) %>%
       filter(!is.na(value)) %>%
       group_by(variable) %>%
       reframe(., 
            mean = mean(value),
            std_dev = sd(value),
            min = min(value),
            q25 = quantile(value, 0.25),
            median = median(value),
            q75 = quantile(value, 0.75),
            max = max(value),
            n = n()) %>%
       mutate(age_category = "all-cohort") %>%
       mutate(copd = "patients&controls")

# the only categorical variables should be sex (F - M), tobacco status (current - never - former),  copd (0 - 1) and maternal_smoking (0 - 1) - if applicable
summary_cat_all <- clinical_df %>%
       dplyr::select(-age_category, -age_squared, -height_meters_squared) %>%
       melt(., id.vars  = c("ID")) %>% 
       subset(variable %in% factor_vars) %>%
       filter(!is.na(value)) %>%
       mutate(variable = as.factor(as.character(variable)),
              value = as.character(value)) %>%
       group_by(variable, value) %>%
       summarise(., n = n()) %>%
       group_by(variable) %>%
       mutate(perc = 100 * n/sum(n)) %>%
       mutate(age_category = "all-cohort") %>%
       mutate(copd = "patients&controls")
            

# __________________________________________________________________________________________________
# Separating by age bins

summary_num_agebins <- 
  clinical_df %>%
  dplyr::select(-age_squared, -height_meters_squared) %>%
  melt(., id = c("ID", "age_category")) %>%
  subset(variable %in% numeric_vars) %>%
  mutate(variable = as.factor(as.character(variable)),
         value = as.numeric(value)) %>%
  filter(!is.na(value)) %>%
  group_by(age_category, variable) %>%
  reframe(., 
          mean = mean(value),
          std_dev = sd(value),
          min = min(value),
          q25 = quantile(value, 0.25),
          median = median(value),
          q75 = quantile(value, 0.75),
          max = max(value),
          n = n()) %>%
  mutate(copd = "patients&controls")


# the only categorical variables should be sex (F - M) and tabaco status (current - never - former)
summary_cat_agebins <- clinical_df %>%
       dplyr::select(-age_squared, -height_meters_squared) %>%
       melt(., id = c("ID", "age_category")) %>%
       subset(variable %in% factor_vars) %>%
       mutate(variable = as.factor(as.character(variable)),
              value = as.character(value)) %>%
       filter(!is.na(value)) %>%
       group_by(age_category, variable, value) %>%
       summarise(., n = n()) %>%
       group_by(age_category, variable) %>%
       mutate(perc = 100 * n/sum(n)) %>%
       mutate(copd = "patients&controls")


# __________________________________________________________________________________________________
# Separating by age bins and COPDs

summary_num_agebins_copds <-
  clinical_df %>%
  dplyr::select(-age_squared, -height_meters_squared) %>%
  melt(., id = c("ID", "age_category", "copd")) %>%
  subset(variable %in% numeric_vars) %>%
  mutate(variable = as.factor(as.character(variable)),
         value = as.numeric(value)) %>%
  filter(!is.na(value)) %>%
  group_by(age_category, copd, variable) %>%
  reframe(.,
          mean = mean(value),
          std_dev = sd(value),
          min = min(value),
          q25 = quantile(value, 0.25),
          median = median(value),
          q75 = quantile(value, 0.75),
          max = max(value),
          n = n()) %>%
  mutate(copd = as.factor(copd))



# the only categorical variables should be sex (F - M) and tabaco status (current - never - former)
summary_cat_agebins_copds <- clinical_df %>%
       dplyr::select(-age_squared, -height_meters_squared) %>%
       melt(., id = c("ID", "age_category", "copd")) %>%
       subset(variable %in% factor_vars) %>%
       mutate(variable = as.factor(as.character(variable)),
              value = as.character(value)) %>%
       filter(!is.na(value)) %>%
       group_by(age_category, copd, variable, value) %>%
       summarise(., n = n()) %>%
       group_by(age_category, copd, variable) %>%
       mutate(perc = 100 * n/sum(n)) %>%
       mutate(copd = as.factor(copd))



# __________________________________________________________________________________________________
# Save results

numeric_res <- bind_rows(summary_num_all, summary_num_agebins, summary_num_agebins_copds)
categ_res <- bind_rows(summary_cat_all, summary_cat_agebins, summary_cat_agebins_copds)

write.csv(numeric_res, file = "../results/summary-stats/cohort_characteristics_num.csv", row.names = F)
write.csv(categ_res, file = "../results/summary-stats/cohort_characteristics_cat.csv", row.names = F)



# _________________________________________________________________________________________________

# Check how severe are COPD patients based on GOLD
# Get number of COPD patients for each GOLD category group

n_copds_per_age_cat <-  clinical_df %>%
       subset(copd == 1) %>%
       # define GOLD categories for COPD patients
        mutate(GOLD = NA,
               GOLD = ifelse(fev_perc >= 80, "GOLD-1", GOLD),
               GOLD = ifelse((fev_perc >= 50) & (fev_perc < 80), "GOLD-2", GOLD),
               GOLD = ifelse((fev_perc >= 30) & (fev_perc < 50), "GOLD-3", GOLD),
               GOLD = ifelse(fev_perc < 30, "GOLD-4", GOLD)) %>%
       group_by(age_category, GOLD) %>%
       summarise(n = n())

write.csv(n_copds_per_age_cat, file = "../results/summary-stats/n_copds_per_age_cat.csv", row.names = T)


