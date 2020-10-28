setwd('/Users/raphaelsinclair/Desktop/MSc Health Data Analytics - IC/HDA/SPH028 - Computational Epidemiology/Project/')

# Packages
library(dplyr)
library(lme4)

# Load data
covars <- readRDS('Data/Covariates.rds')
transcripts <- readRDS('Data/Transcripts.rds')     # Very large
exposures <- readRDS('Data/Exposures.rds')

###### Data preparation ####################

# Remove variables with too much missing data
missing_list <- NULL
for (i in 1:ncol(transcripts)) {
  if ( (sum(is.na(transcripts[,i]))/length(transcripts[,i])) > 0.95) {
    missing_list <- c(missing_list, i)
  }
}

transcripts <- transcripts[, -missing_list]

# ------------------
# SUBSETTING DATA. FILE TOO LARGE. COMMENT OUT FOR HPC
s = sample(ncol(transcripts), size = 10)
transcripts <- transcripts[,s]
# ------------------


# Log transform exposure and OMICs (where necessary)
exposures[, 3:ncol(exposures)] <- log(exposures[, 3:ncol(exposures)])


# Merge covariates and exposure dataframes then remove duplicate columns
covar_exp_merge <- merge(covars, exposures, by = 0)
rownames(covar_exp_merge) <- covar_exp_merge[,1]
covar_exp_merge <- subset(covar_exp_merge, select = -c(subjectidp.y, id.y, Row.names))
covar_exp_merge <- covar_exp_merge %>% rename(subjectidp = subjectidp.x, id = id.x)

# Filter merged dataset to match all observations in OMIC data
samples <- rownames(transcripts)
covar_exp_merge <- filter(covar_exp_merge, subjectidp %in% samples)
transcripts <- filter(transcripts, rownames(transcripts) %in% covar_exp_merge$subjectidp)


###### Obtain vector of differences in ICCs for ID ####################
## Use linear mixed models to isolate the variation within IDs (variability within time points)

# Run linear mixed models 
tol = 10^(-4)    # Set a tolerance parameter
icc_diffs <- rep(0, ncol(transcripts))     # Create item to store differences in ICCs
names(icc_diffs) <- colnames(transcripts)

fixed_eff <- paste0(c('age', 'gender', 'bmi'), collapse = '+')     # List of fixed effects
rand_eff_list <- c('(1 | id)', '(1 | isolation)', '(1| labeling)', '(1| hybridization)')     # List of random effects (extracted in mixed models)
rand_eff <- paste0(rand_eff_list, collapse = '+')

for (column in 1:ncol(transcripts)) {
  y <- paste0('transcripts[, ', column, '] ~ ')
  f1 <- paste0(y, rand_eff, '+', fixed_eff)     # formula WITHOUT exposure
  f2 <- paste0(f1, '+', 'pncmedian')     # formula WITH exposure
  
  model1 = lmer(as.formula(f1), data = covar_exp_merge)
  vcov1 = as.data.frame(VarCorr(model1))$vcov
  names(vcov1) = as.data.frame(VarCorr(model1))[, 1]
  
  model2 = lmer(as.formula(f2), data = covar_exp_merge)
  vcov2 = as.data.frame(VarCorr(model2))$vcov
  names(vcov2) = as.data.frame(VarCorr(model2))[, 1]
  
  icc_diffs[column] = (vcov1[1]/sum(vcov1)) - (vcov2[1]/sum(vcov2))
}
  
  
###### Create null distribution of differences in ICCs ######
## Permute exposure to break association with exposure and create distribution of null ICCs

temp <- rep(0,ncol(transcripts))     # temp vector to store 'null' ICC differences
null_icc_diff_matrix <- NULL     # Initiate matrix to store 'null' ICC differences over iterations
iterations = 10

# Run linear mixed models on permuted data
for (i in 1:iterations) {
  permuted <- sample(covar_exp_merge$pncmedian)     # Permute exposure 
  covar_exp_merge$permuted <- permuted
  
  for (column in 1:ncol(transcripts)) {
    
    y <- paste0('transcripts[, ', column, '] ~ ')
    f1 <- paste0(y, rand_eff, '+', fixed_eff)
    f2 <- paste0(f1, '+', 'permuted')
    
    model1 = lmer(as.formula(f1), data = covar_exp_merge)
    vcov1 = as.data.frame(VarCorr(model1))$vcov
    names(vcov1) = as.data.frame(VarCorr(model1))[, 1]
    
    model2 = lmer(as.formula(f2), data = covar_exp_merge)
    vcov2 = as.data.frame(VarCorr(model2))$vcov
    names(vcov2) = as.data.frame(VarCorr(model2))[, 1]
    
    
    temp[column] = (vcov1[1]/sum(vcov1)) - (vcov2[1]/sum(vcov2))
  }
  
  null_icc_diff_matrix <- cbind(null_icc_diff_matrix, temp)
}

# Assign row/column names to ICC matrix
rownames(null_icc_diff_matrix) <- colnames(transcripts)
colnames(null_icc_diff_matrix) <- paste0('iter ', seq(ncol(null_icc_diff_matrix)))


###### Calculate p-value for observed ICCs ######
## Null hypothesis is that there is NO association between OMIC stability and exposure.
## Need to compare observed differences in ICCs to null distribution of differences

# Rank null ICC and obtain p-values for observed ICCs
exp_diff_trans <- rep(0, length(icc_diffs))
names(exp_diff_trans) <- names(icc_diffs)

for (column in 1:length(icc_diffs)) {
  exp_diff_trans[column] = sum(null_icc_diff_matrix[column,] >= icc_diffs[column])/length(null_icc_diff_matrix[column,])
}


###### Save p-values ######
ifelse(dir.exists("Results"),"",dir.create("Results", showWarnings = FALSE))
saveRDS(exp_diff_trans, file = 'Results/exp_diff_trans.rds')
