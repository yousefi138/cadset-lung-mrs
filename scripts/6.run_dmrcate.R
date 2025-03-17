####################################################################################################

## Project:        CADSET 2.0 - COPD epigenetics across the lifespan
## File name:      run_dmrcate.R
## Description:    script to test differentially methylated regions using DMRcate tool
## Authors:        Julieta Viglino
##                 Sandra Casas
## Version:    	   2023-11
## R version used: R 4.2.2 (for the design and testing of this script)

####################################################################################################

# IMPORTANT: THIS SCRIPT SHOULD ONLY BE IF YOUR COHORT INCLUDES ADOLESCENTS/ADULTS
## IF ITS A ONLY CHILDREN COHORT (AGE < 12 YEARS), DONT NEED TO RUN THIS SCRIPT

library("DMRcate") # install through Bioconductor
library("minfi") # install through Bioconductor
library("dplyr")


rm(list=ls())
gc()

# NOTE 1: if during your cohort-specific methylation analysis batch was corrected while performing the models
# (i.e using negative control PCs or batch as model covariate), this should also be included in the formula!
# Please add your specific batch variables to models formulas is applicable

# NOTE 2: Modify the cell populations in the formulas if you have different ones


# load clinical_data that was already created in the script create_clinical_df.R
clinical_df <-  read.csv("../data/clinical_df.csv") %>%
  # to be run only for adolescents/adults, children will be excluded:
  subset(age > 12) %>%
  mutate(sex = factor(sex, levels = c("F", "M")),
         smoking_status = factor(smoking_status, levels = c("NS", "FS", "CS")),
         age_category = as.factor(age_category),
         pack_years = as.numeric(pack_years))
rownames(clinical_df) <- clinical_df$ID

# load methylation data
# note: user should complete with the path were you methylation data is. note that probes should be rows and patients columns
methyl_data <- readRDS("../data/methylation_matrix.RDS")
# arrange methylation matrix so that patients are in the same order than the clinical df
methyl_data <- methyl_data[,as.character(clinical_df$ID)]

# make sure the order of the IDs match
table(rownames(clinical_df) == colnames(methyl_data))

# obtain M values from Beta values:
methyl_data_M <- logit2(methyl_data)
methyl_data_M <- as.matrix(methyl_data_M)

# ---------------------------------------------------------------------------------------------------------------------------------- #

#### RUN MODEL WITHOUT CELL POPULATIONS USING ALL COHORT ####

# include batch variables if applicable (i.e PC1 + PC2 + PC3 ...)
formula_ <- " ~ fevml + age + age_squared + sex + height_meters + height_meters_squared + smoking_status + pack_years" 
design_matrix <- model.matrix(as.formula(formula_), data = clinical_df)
set.seed(1)
tryCatch(
  
  { 
    myannotation_fevml <- cpg.annotate(methyl_data_M,
                                          datatype = "array", 
                                          arraytype = "EPIC",  # "EPIC" or "450K" choose depending on the type of array
                                          what = "M",
                                          analysis.type="differential", 
                                          design=design_matrix, 
                                          coef = "fevml")
    dmrcoutput_fevml <- dmrcate(myannotation_fevml, lambda=1000, C=2)
    dmrcoutput_fevml.ranges <- extractRanges(dmrcoutput_fevml, genome = "hg19") %>% as_tibble() 
    write.csv(dmrcoutput_fevml.ranges, "../results/dmr/dmr_fevml_results_without_cell_populations.csv", row.names = F)
    rm(dmrcoutput_fevml, myannotation_fevml, design_matrix, dmrcoutput_fevml.ranges)
  },
  
  error = function(e){ print(e) },
  warning  = function(w){ print(w) }
)

gc()

# ---------------------------------------------------------------------------------------------------------------------------------- #

#### RUN DMR CATE MODEL WITH CELL POPULTIONS USING ALL COHORT ####

# include batch variables and modify cell populations if applicable:
formula_cp <- " ~ fevml + age + age_squared + sex + height_meters + height_meters_squared + smoking_status + pack_years +  CD8T + CD4T + NK + Bcell + Neu " 
design_matrix_cp <- model.matrix(as.formula(formula_cp), data = clinical_df)

set.seed(1)
tryCatch(
  
  { 
    myannotation_fevml_cp <- cpg.annotate(methyl_data_M,
                                          datatype = "array", 
                                          arraytype = "EPIC", # "EPIC" or "450K" choose depending on the type of array
                                          what = "M",
                                          analysis.type="differential", 
                                          design=design_matrix_cp, 
                                          coef = "fevml")
    dmrcoutput_fevml_cp <- dmrcate(myannotation_fevml_cp, lambda=1000, C=2)
    dmrcoutput_fevml_cp.ranges <- extractRanges(dmrcoutput_fevml_cp, genome = "hg19") %>% as_tibble() 
    write.csv(dmrcoutput_fevml_cp.ranges, "../results/dmr/dmr_fevml_results_with_cell_populations.csv", row.names = F)
    rm(dmrcoutput_fevml_cp, myannotation_fevml_cp, design_matrix_cp, dmrcoutput_fevml_cp.ranges)
  },
  
  error = function(e){ print(e); message("Could not obtain DMR regions")},
  warning  = function(w){ print(w) }
)


gc()
