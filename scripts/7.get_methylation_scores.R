
####################################################################################################

## Project:        CADSET 2.0 - COPD epigenetics across the lifespan
## File name:      get_methylation_risk_scores.R
## Description:    script to obtain methylation risk scores for each subject
##                 uses three different scores: one associated to FEVlt in population based studies
##                 and two scores associated to FEVlt and severity in a COPD clinical cohort.
## Authors:        Julieta Viglino
##                 Sandra Casas
## Version:    	   2023-11
## R version used: R 4.2.2 (for the design and testing of this script)

####################################################################################################

rm(list = ls())
library(dplyr)

# load cpgs associated to each of the methylation risk scores (three different csvs files)
# these files were provided by us in the mrs/ subdirectory from the shared working directory folder
mrs_fev_pop <- read.csv("mrs/mrs_population_fevlt.csv", row.names = 1) %>% mutate("probe" = rownames(.))
mrs_fev_clin <- read.csv("mrs/mrs_copd_fevlt.csv", row.names = 1) %>% mutate("probe" = rownames(.))
mrs_severity <- read.csv("mrs/mrs_copd_severity.csv", row.names = 1) %>% mutate("probe" = rownames(.))


# load methylation data (beta values)
methyl_data <-  readRDS("../data/methylation_matrix.RDS")

### function to calculate methylation risk scores

# this function calculate MRS using both standarized and non standarized coefficients
# it returns a dataframe with patients ids and the three mrs
# in case any cpg probe from the MRS is not found in the methylation matrix, those cpgs names are also returned.

getMRS <-  function(methyl_data, mrs){
  
  keep_cpgs <-  intersect(rownames(methyl_data), rownames(mrs))
  cpgs_used <-  data.frame(cpgs_mrs = rownames(mrs)) %>%
    mutate(used = ifelse(cpgs_mrs %in% keep_cpgs, 1, 0))

  methyl_data <-  methyl_data[keep_cpgs,]
  mrs <-  mrs[keep_cpgs,]
  
  res <-  t(methyl_data) %*% mrs[,"effect"] %>%
    as.data.frame() %>% 
    mutate(ID = rownames(.)) %>%
    dplyr::rename("MRS" = 1)

  res <- res %>%
    mutate(std_MRS = (MRS-mean(MRS)) / sd(MRS)) %>%
    select(ID, MRS, std_MRS)
  

  cor_ <-  cor(t(methyl_data), method = "spearman", use = "pairwise.complete.obs")

  return(list(results = res,
              cpgs_used = cpgs_used,
              cpgs_cor = cor_))
}



res_mrs_fev_pop <-  getMRS(methyl_data, mrs_fev_pop)
res_mrs_fev_clin <-  getMRS(methyl_data, mrs_fev_clin)
res_mrs_sev_clin <-  getMRS(methyl_data, mrs_severity)

mrs_results <-  res_mrs_fev_pop$results %>% rename_at(2:3, ~ paste(., "fev_pop", sep = "_")) %>%
        left_join(., (res_mrs_fev_clin$results %>% rename_at(2:3, ~ paste(., "fev_clin", sep = "_")))) %>%
        left_join(., (res_mrs_sev_clin$results %>% rename_at(2:3, ~ paste(., "sev", sep = "_"))))
                        

write.csv(mrs_results, file = "../data/mrs_results.csv", row.names=F)

write.csv(res_mrs_fev_pop$cpgs_used, file = "../results/mrs/cpgs_used_in_fev_population_mrs.csv", row.names=F)
write.csv(res_mrs_fev_clin$cpgs_used, file = "../results/mrs/cpgs_used_in_fev_clinical_mrs.csv", row.names=F)
write.csv(res_mrs_sev_clin$cpgs_used, file = "../results/mrs/cpgs_used_in_severity_mrs.csv", row.names=F)

write.csv(res_mrs_fev_pop$cpgs_cor, file = "../results/mrs/cpgs_cor_fev_population_mrs.csv")
write.csv(res_mrs_fev_clin$cpgs_cor, file = "../results/mrs/cpgs_cor_fev_clinical_mrs.csv")
write.csv(res_mrs_sev_clin$cpgs_cor, file = "../results/mrs/cpgs_cor_severity_mrs.csv")


