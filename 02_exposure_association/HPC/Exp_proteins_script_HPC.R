setwd('/Users/raphaelsinclair/Desktop/MSc Health Data Analytics - IC/HDA/SPH028 - Computational Epidemiology/Project/')

# Packages
library(dplyr)
library(lme4)

# Load data
covars <- readRDS('Data/Covariates.rds')
proteins <- readRDS('Data/Proteins.rds')
exposures <- readRDS('Data/Exposures.rds')

###### Data preparation ####################

# Log transform exposure and OMICs (where necessary)
exposures[, 3:ncol(exposures)] <- log(exposures[, 3:ncol(exposures)])
proteins <- log(proteins)

# Merge covariates and exposure dataframes then remove duplicate columns
covar_exp_merge <- merge(covars, exposures, by = 0)
rownames(covar_exp_merge) <- covar_exp_merge[,1]
covar_exp_merge <- subset(covar_exp_merge, select = -c(subjectidp.y, id.y, Row.names))
covar_exp_merge <- covar_exp_merge %>% rename(subjectidp = subjectidp.x, id = id.x)

# Filter merged dataset to match all observations in OMIC data
samples <- rownames(proteins)
covar_exp_merge <- filter(covar_exp_merge, subjectidp %in% samples)
proteins <- filter(proteins, rownames(proteins) %in% covar_exp_merge$subjectidp)


###### Obtain vector of differences in ICCs for ID ####################
## Use linear mixed models to isolate the variation within IDs (variability within time points)

# Run linear mixed models 
tol = 10^(-4)    # Set a tolerance parameter
icc_diffs <- rep(0, ncol(proteins))     # Create item to store differences in ICCs
names(icc_diffs) <- colnames(proteins)

fixed_eff <- paste0(c('age', 'gender', 'bmi'), collapse = '+')     # List of fixed effects
rand_eff_list <- c('(1 | id)', '(1 | plate)')     # List of random effects (extracted in mixed models)
rand_eff <- paste0(rand_eff_list, collapse = '+')

for (column in 1:ncol(proteins)) {
  y <- paste0('proteins[, ', column, '] ~ ')
  f1 <- paste0(y, rand_eff, '+', fixed_eff)     # formula WITHOUT exposure
  f2 <- paste0(f1, '+', 'pncmedian')     # formula WITH exposure
  
  model1 = lmer(as.formula(f1), data = covar_exp_merge)
  vcov1 = as.data.frame(VarCorr(model1))$vcov
  names(vcov1) = as.data.frame(VarCorr(model1))[, 1]
  
  model2 = lmer(as.formula(f2), data = covar_exp_merge)
  vcov2 = as.data.frame(VarCorr(model2))$vcov
  names(vcov2) = as.data.frame(VarCorr(model2))[, 1]
  
  if (sum(vcov1[2:(length(vcov1)-1)] < tol) > 0) {
    vcov1 <- vcov1[-length(vcov1)]
    new_rand_eff <- rand_eff_list
    new_rand_eff <- new_rand_eff[vcov1 > tol]
    new_rand_eff <- paste0(new_rand_eff, collapse = '+')
    f = paste0(y, new_rand_eff, '+', fixed_eff)
    
    model1 = lmer(as.formula(f), data = covar_exp_merge)
    vcov1 = as.data.frame(VarCorr(model1))$vcov
    names(vcov1) = as.data.frame(VarCorr(model1))[, 1]
  }
  
  if (sum(vcov2[2:(length(vcov2)-1)] < tol) > 0) {
    vcov2 <- vcov2[-length(vcov2)]
    new_rand_eff <- rand_eff_list
    new_rand_eff <- new_rand_eff[vcov2 > tol]
    new_rand_eff <- paste0(new_rand_eff, collapse = "+")
    f = paste0(y, new_rand_eff, '+', fixed_eff)
    
    model2 = lmer(as.formula(f), data = covar_exp_merge)
    vcov2 = as.data.frame(VarCorr(model2))$vcov
    names(vcov2) = as.data.frame(VarCorr(model2))[, 1]
  }
  
  icc_diffs[column] = (vcov1[1]/sum(vcov1)) - (vcov2[1]/sum(vcov2))
}
  
  
###### Create null distribution of differences in ICCs ######
## Permute exposure to break association with exposure and create distribution of null ICCs

temp <- rep(0,ncol(proteins))     # temp vector to store 'null' ICC differences
null_icc_diff_matrix <- NULL     # Initiate matrix to store 'null' ICC differences over iterations
iterations = 500

# Run linear mixed models on permuted data
for (i in 1:iterations) {
  permuted <- sample(covar_exp_merge$pncmedian)     # Permute exposure 
  covar_exp_merge$permuted <- permuted
  
  for (column in 1:ncol(proteins)) {
    
    y <- paste0('proteins[, ', column, '] ~ ')
    f1 <- paste0(y, rand_eff, '+', fixed_eff)
    f2 <- paste0(f1, '+', 'permuted')
    
    model1 = lmer(as.formula(f1), data = covar_exp_merge)
    vcov1 = as.data.frame(VarCorr(model1))$vcov
    names(vcov1) = as.data.frame(VarCorr(model1))[, 1]
    
    model2 = lmer(as.formula(f2), data = covar_exp_merge)
    vcov2 = as.data.frame(VarCorr(model2))$vcov
    names(vcov2) = as.data.frame(VarCorr(model2))[, 1]
    
    if (sum(vcov1[2:(length(vcov1)-1)] < tol) > 0) {
      vcov1 <- vcov1[-length(vcov1)]
      new_rand_eff <- rand_eff_list
      new_rand_eff <- new_rand_eff[vcov1 > tol]
      new_rand_eff <- paste0(new_rand_eff, collapse = '+')
      f = paste0(y, new_rand_eff, '+', fixed_eff)
      
      model1 = lmer(as.formula(f), data = covar_exp_merge)
      vcov1 = as.data.frame(VarCorr(model1))$vcov
      names(vcov1) = as.data.frame(VarCorr(model1))[, 1]
    }
    
    if (sum(vcov2[2:(length(vcov2)-1)] < tol) > 0) {
      vcov2 <- vcov2[-length(vcov2)]
      new_rand_eff <- rand_eff_list
      new_rand_eff <- new_rand_eff[vcov2 > tol]
      new_rand_eff <- paste0(new_rand_eff, collapse = "+")
      f = paste0(y, new_rand_eff, '+', fixed_eff)
      
      model2 = lmer(as.formula(f), data = covar_exp_merge)
      vcov2 = as.data.frame(VarCorr(model2))$vcov
      names(vcov2) = as.data.frame(VarCorr(model2))[, 1]
    }
    
    temp[column] = (vcov1[1]/sum(vcov1)) - (vcov2[1]/sum(vcov2))
  }
  
  null_icc_diff_matrix <- cbind(null_icc_diff_matrix, temp)
}

# Assign row/column names to ICC matrix
rownames(null_icc_diff_matrix) <- colnames(proteins)
colnames(null_icc_diff_matrix) <- paste0('iter ', seq(ncol(null_icc_diff_matrix)))


###### Calculate p-value for observed ICCs ######
## Null hypothesis is that there is NO association between OMIC stability and exposure.
## Need to compare observed differences in ICCs to null distribution of differences

# Rank null ICC and obtain p-values for observed ICCs
exp_diff_pro <- rep(0, length(icc_diffs))
names(exp_diff_pro) <- names(icc_diffs)

for (column in 1:length(icc_diffs)) {
  exp_diff_pro[column] = sum(null_icc_diff_matrix[column,] >= icc_diffs[column])/length(null_icc_diff_matrix[column,])
}


###### Save p-values ######
ifelse(dir.exists("Results"),"",dir.create("Results", showWarnings = FALSE))
saveRDS(exp_diff_pro, file = 'Results/exp_diff_pro.rds')
