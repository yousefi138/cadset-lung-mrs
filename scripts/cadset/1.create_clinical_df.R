################################################################################################################

## Project:        CADSET 2.0 - COPD epigenetics across the lifespan
## File name:      create_clinical_df.R
## Description:    code to create a dataframe with all the clinical data that will be used in the project
## Authors:        Julieta Viglino
##                 Sandra Casas
## Version:    	   2023-11
## R version used: R 4.2.2 (for the design and testing of this script)

################################################################################################################

## IMPORTANT NOTE: The scripts provided are general and require the user to load specific studies database and rename variables names accordingly

# install.packages("dplyr")
library(dplyr)


rm(list=ls())
gc()

# read clinical dataset
clinical_df_path <- "" # note: User should provide clinical_df_path (i.e: "/home/user/Documents/databases/clinical_data.csv")
clinical_df <- read.csv(clinical_df_path, sep = ",") # add separator if needed

# read cell populations dataset if not already included in clinica dataset
cell_pop_path <- "" # note: user should provide cell populations dataset path.
cell_pop_df <- read.csv(cell_pop_path, sep = ",")

# if applicable, merge clinical dataset and cell populations dataset by ID
# note: make sure both clinical_df and cell_pop_df have a column named ID matching samples ID.
clinical_df <- clinical_df %>% left_join(cell_pop_df, by = "ID")

# load methylation data (*** beta values ***)
# note: user should complete with the path were you methylation data is. note that probes should be rows and patients columns
# imporant: if you have M values, please convert them to Beta values
methyl_path <- "" # note: user should provide methylation matrix path.
methyl_data <- readRDS(methyl_path)

# keep patients for which methyl assay was performed
# before runnning make sure that the ID coding is the same in clinical_df and methylation matrix : i.e patient X123 vs X-123
clinical_df <- clinical_df %>% subset(ID %in% colnames(methyl_data))  
dim(clinical_df)


####  PREPROCESS CLINICAL DATASET  #####

# create clinical data set with standarized variable names
# user will have to rename and probably recode some variables

## Important Notes:
# note 1: if pre and pro bronquidilatador measurments are available, include post
# note 2: make sure fev unit is ml and height unit is meters (otherwise see lines 94-100)
# note 3: if during your specific cohort data modeling, technical variation effects were included in models (as PCs or as categorical technical variable),  include those in the clinical data table
# note 4: there should be an ID variable named ID
# note 5: tabaco-use will be added into analysis in two separate ways depending in patients in a child (<12 yearls-old) or adolescent/adult (>12 years old)
### for children, maternal smoking will be used and should indicate if the mother was a smoker during pregnancy (variable: maternal_smoking)
### for adolescents/adults, self smoking will be used (variables: smoking_status -former, current, never smoker- and pack-years smoked)

clinical_df <-  clinical_df %>%

  # all the variables below should be present in your clinical_df (otherwise please let us know)
  # if some variables are not present in your dataset, create them using mutate if possible (i.e mutate(age_squared = age^2))

  # rename variables names accordingly: ("new_var_name" = "old_var_name" > note that new_var_name is already provided)
  dplyr::rename("fevml" = "fevml", # i.e : "fevlt" = "FEV1LT"  # note: in ml
         "fvcml" = "fvcml",  # note: in ml
         "fev_fvc_ratio" = "FEV1_FVC", # range 0-1
         "fev_perc" = "FEV1", # FEV1% in percentage
         "age" = "age",
         "sex" = "sex", # sex should be "F" when females or "M" if males, use line 99 to modify if sex is not encoded like F/M
         "height_meters" = "height", # note: in meters

         # tabacco exposure for adults or adolescents (age greater than 12 years):
         "smoking_status" = "smoking", # former (FS), current (CS), never (NS), use line 103-109 to modify if it is not encoded like this
         "pack_years" = "pack_year", # never smokers are expected to have pack years = 0 (if never smokers have pack_years = NA, impute to pack_years = 0)

         # smoking exposure for childs (age 12 years or less):
         "maternal_smoking" = "", # maternal smoking during pregnancy: 1 if smoker or 0 if non-smoker, not applicable if no childs are present in cohort

         # note that maternal_smoking should be NA for adolescents/adults, while smoking_status and pack_years should be NA for children (12 years or less)

         # deconvoluted or lab measured cell counts (whatever populations were used in original analysis)
         # these will be cohort specific: if some of these cells types are not present, or you have others not listed here, modify accordingly:
         "CD8T" = "CD8T",
         "CD4T" = "CD4T",
         "NK" = "NK",
         "Bcell" = "Bcell",
         "Neu" = "Neu",
         "Mono" = "Mono") %>%

  # !! IMPORTANT: ADD BATCH VARIABLES (PCs, Sample_Plate, Sample_Slide) if these needs to be used as covariate when performing linear models.

  # IF APPLICABLE:
  # mutate(fevml = fevml*1000) %>% # if your FEV variable was in litres instead of ml
  # mutate(fvcml = fvcml*1000) %>% # if your FVC variable was in litres instead of ml
  # mutate(fev_fvc_ratio = fev_fvc_ratio/100) %>%  # if fev_fvc_ratio was in 0-100 range instead of 0-1 range
  # mutate(height_meters = height_meters /100) %>% # if meters was in cm
  # mutate(sex = ifelse(sex == 1, "F", ifelse(sex == 0, "M", NA))) %>% # note: modify accordingly to have M and F as coding for sex variable

  mutate(age_squared = age^2,
         height_meters_squared = height_meters^2)  %>%

  # modify smoking_status coding depending on your original data to have NS (never smokers), FS (former smokers) or CS (current smokers)
  # note that all categories (NS, FS, CS) might not be present in your dataset
  # mutate(smoking_status = case_when(smoking_status == 1 ~  "NS",
  #                                   smoking_status == 2 ~  "FS",
  #                                   smoking_status == 3 ~  "CS"))  %>%

  # categorize age in age ranges
  mutate(age_category = NA,
         age_category = ifelse(age <= 0, "coord_blood", age_category),
         age_category = ifelse((age > 0) & (age <= 12), "children", age_category),
         age_category = ifelse((age > 12) & (age <= 20), "adolescence", age_category),
         age_category = ifelse((age > 20) & (age <= 40), "young-adulthood", age_category),
         age_category = ifelse((age > 40) & (age <= 60), "adulthood-1", age_category),
         age_category = ifelse((age > 60) & (age <= 80), "adulthood-2", age_category),
         age_category = ifelse(age > 80, "seniors", age_category)) %>%


  mutate(sex = as.factor(sex),
         smoking_status = as.factor(smoking_status),
         maternal_smoking = as.factor(maternal_smoking), # Comment line if childs are not present in cohort (childs: 12 years-old or less)
         age_category = as.factor(age_category),
         pack_years = as.numeric(pack_years)) %>%

  # create copd variable - based on ratio value
  mutate(copd = as.factor(ifelse(fev_fvc_ratio < 0.7, 1, 0))) %>% # check that fev_fvc_ratio is in 0-1 scale and not 0-100

  # select only relevant variables         
  dplyr::select(ID, fevml, fev_fvc_ratio, fvcml, fev_perc, copd,
         age, age_squared, age_category, sex,
         CD8T, CD4T, NK, Bcell, Neu, Mono, # replace by your specific cell population
         height_meters, height_meters_squared,
         maternal_smoking,   # only if childs are present in cohort (12 years-old or less), otherwise comment line
         smoking_status, pack_years)


# get total number of subjects per age bin before filtering patients with NA in key variables:
totNsubjects <- clinical_df %>%
  group_by(age_category) %>%
  summarise(n_subjects = n())
write.csv(totNsubjects,  file = "../results/summary-stats/n_per_age_category_before_filtering.csv", row.names = T)



# we will filter data set in order to keep only subjects that do not have NA in the variables that will be used for linear models
clinical_df <- clinical_df %>%
  filter(!is.na(fevml) & !is.na(age) & !is.na(age_squared) &
         !is.na(sex) & !is.na(height_meters) & !is.na(height_meters_squared) &
         # modify line below depending on the specific cell populations types that you have
         !is.na(CD8T) & !is.na(CD4T) & !is.na(NK) & !is.na(Bcell) & !is.na(Neu) & !is.na(Mono))  %>%

  # smoking status will also be used in analysis, but filtering of missing data will be conditioned to subjects being childs or adolescents/adults:
  # if patient is a baby/child (age <= 12), "maternal_smoking" should not have NAs:
  filter((age <= 12 & !is.na(maternal_smoking)) | age > 12) %>%
  # if patients is older than 12, "smoking_status" and "pack_years" should not have NAs:
  filter((age > 12 & !is.na(smoking_status) & !is.na(pack_years)) | age <= 12)


dim(clinical_df)


# after filtering subjects that have NA in relevant features, make sure you also subset the methylation matrix:
methyl_data <- methyl_data[, as.character(clinical_df$ID)]



saveRDS(methyl_data, "../data/methylation_matrix.RDS")
write.csv(clinical_df, "../data/clinical_df.csv", row.names = F)  # NOT TO SHARE, used only for running the rest of the scripts
