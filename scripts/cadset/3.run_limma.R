####################################################################################################

## Project:        CADSET 2.0 - COPD epigenetics across the lifespan
## File name:      run_limma.R
## Description:    script to test differentially methylated position using limma tool
## Authors:        Julieta Viglino
##                 Sandra Casas
## Version:    	   2023-11
## R version used: R 4.2.2 (for the design and testing of this script)

####################################################################################################

# IMPORTANT NOTE: The scripts provided are general and require the user to load specific studies database and rename variables names accordingly
# Please, first use the option 'install.packages' to install the libraries if any have not previously installed

rm(list=ls())
gc()


library("limma")
library("dplyr")

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

# load methylation data (beta values)
methyl_data <- readRDS("../data/methylation_matrix.RDS")
# arrange methylation matrix so that patients are in the same order than the clinical df
methyl_data <- methyl_data[,as.character(clinical_df$ID)]

# make sure the order of the IDs match
rownames(clinical_df) <- clinical_df$ID
table(rownames(clinical_df) == colnames(methyl_data))


# ---------------------------------------------------------------------------------------------------------------------------------- #

#### 1) RUN LIMMA MODEL WITHOUT CELL POPULATIONS USING ALL COHORT (ADOLESCENTS/ADULTS) ####
### if your cohort is from children only, no need to run this part - you can go directly to step 3 ##
# IMPORTANT: if applicable, add technical variation (sample_plate or PCs)

# subset adults
clinical_df_adults <- clinical_df %>% subset(age > 12)
methyl_data_adults <- methyl_data[,as.character(clinical_df_adults$ID)]

formula_ <- " ~ fevml + age + age_squared + sex + height_meters + height_meters_squared + smoking_status + pack_years" # add your specific batch variables if the beta matrix is not batch adjusted.

design_matrix <- model.matrix(as.formula(formula_), data = clinical_df_adults)

set.seed(1)
fit <- lmFit(methyl_data_adults, design_matrix)
fit <- eBayes(fit)
res_fit <- data.frame(
      probe = rownames(fit$coefficients),
      mean_methyl = fit$Amean,
      fevml_coef = fit$coefficients[, "fevml"],
      fevml_std_error =  (sqrt(fit$s2.post) * fit$stdev.unscaled)[, "fevml"],
      fevml_p_val = fit$p.value[, "fevml"]
      ) 

write.csv(res_fit, "../results/limma/limma_results_without_cell_populations.csv", row.names = F)

# ---------------------------------------------------------------------------------------------------------------------------------- #

#### 2) RUN LIMMA MODEL WITH CELL POPULTIONS USING ALL COHORT (ADOLESCENTS/ADULTS) ####
### if your cohort is from children only, no need to run this part - you can go directly to step 3 ##

# IMPORTANT: if applicable, add technical variation (sample_plate or PCs) to formulas and modify cell populations appropriately
formula_cp <- " ~ fevml + age + age_squared + sex + height_meters + height_meters_squared + smoking_status + pack_years +  CD8T + CD4T + NK + Bcell + Neu " # add your specific batch variables if the beta matrix is not batch adjusted.
design_matrix_cp <- model.matrix(as.formula(formula_cp), data = clinical_df_adults)

set.seed(1)
fit_cp <- lmFit(methyl_data_adults, design_matrix_cp)
fit_cp <- eBayes(fit_cp)
res_cp <- data.frame(
  probe = rownames(fit_cp$coefficients),
  mean_methyl = fit_cp$Amean,
  fevml_coef = fit_cp$coefficients[, "fevml"],
  fevml_std_error =  (sqrt(fit_cp$s2.post) * fit_cp$stdev.unscaled)[, "fevml"],
  fevml_p_val = fit_cp$p.value[, "fevml"]
)
write.csv(res_cp, "../results/limma/limma_results_with_cell_populations.csv", row.names = F)

# ---------------------------------------------------------------------------------------------------------------------------------- #

#### 3) RUN LIMMA MODEL PER AGE BINS  ####
# including patients from all age bins (coord blood, children, adults, etc)
# note: include technical variation PCs/Batch effects or modify cell populations in formulas if appropiate!

age_ranges <- unique(clinical_df$age_category)

for (age_range in age_ranges){

  # if there are less than 30 patients in a given age-bin, limma model is not run for that age-bin.
  clinical_df_tmp <- clinical_df %>% subset(age_category == age_range)
  if (nrow(clinical_df_tmp) < 30) {
    next
  }

  # Define model formulas
  # IMPORTANT: if applicable, add technical variation (sample_plate or PCs) to formulas below and modify cell populations accordingly

  # formulas for adolescents/adults (they include self smoking_status & pack_years variable):
  formula_ <- " ~ fevml + age + age_squared + sex + height_meters + height_meters_squared + smoking_status + pack_years" # add your specific batch variables if the beta matrix is not batch adjusted.
  formula_cp <- " ~ fevml + age + age_squared + sex + height_meters + height_meters_squared + smoking_status + pack_years +  CD8T + CD4T + NK + Bcell + Neu " # add your specific batch variables if the beta matrix is not batch adjusted.

  # formulas for children or babies (thay include maternal_smoking as smoking-related variable):
  if ((age_range == "coord_blood") | (age_range == "children")) {
    formula_ <- " ~ fevml + age + age_squared + sex + height_meters + height_meters_squared + maternal_smoking" # add your specific batch variables if the beta matrix is not batch adjusted.
    formula_cp <- " ~ fevml + age + age_squared + sex + height_meters + height_meters_squared + maternal_smoking + CD8T + CD4T + NK + Bcell + Neu " # add your specific batch variables if the beta matrix is not batch adjusted.
  }

  methyl_data_age <- methyl_data[,as.character(clinical_df_tmp$ID)]
  
  # WITHOUT CELL POPULATIONS
  design_matrix_tmp <- model.matrix(as.formula(formula_), data = clinical_df_tmp)
  set.seed(1)
  fit_tmp <- lmFit(methyl_data_age, design_matrix_tmp)
  fit_tmp <- eBayes(fit_tmp)
  res_tmp <- data.frame(
    probe = rownames(fit_tmp$coefficients),
    mean_methyl = fit_tmp$Amean,
    fevml_coef = fit_tmp$coefficients[, "fevml"],
    fevml_std_error =  (sqrt(fit_tmp$s2.post) * fit_tmp$stdev.unscaled)[, "fevml"],
    fevml_p_val = fit_tmp$p.value[, "fevml"]
  )
  write.csv(res_tmp, paste("../results/limma/limma_results_without_cell_populations_", age_range, ".csv", sep = ""), row.names = F)

  # WITH CELL POPULATIONS
  design_matrix_tmp_cp <- model.matrix(as.formula(formula_cp), data = clinical_df_tmp)
  set.seed(1)
  fit_cp_tmp <- lmFit(methyl_data_age, design_matrix_tmp_cp)
  fit_cp_tmp <- eBayes(fit_cp_tmp)
  res_cp_tmp <- data.frame(
    probe = rownames(fit_cp_tmp$coefficients),
    mean_methyl = fit_cp_tmp$Amean,
    fevml_coef = fit_cp_tmp$coefficients[, "fevml"],
    fevml_std_error =  (sqrt(fit_cp_tmp$s2.post) * fit_cp_tmp$stdev.unscaled)[, "fevml"],
    fevml_p_val = fit_cp_tmp$p.value[, "fevml"]
  )
  write.csv(res_cp_tmp, paste("../results/limma/limma_results_with_cell_populations_", age_range, ".csv", sep = ""), row.names = F)

  }
  
  
  
    
  
  
  
  
