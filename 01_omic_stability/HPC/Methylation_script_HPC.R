# ###### Initialisation ####################
# 
# # Local setup
# setwd('/Users/raphaelsinclair/Desktop/MSc Health Data Analytics - IC/HDA/SPH028 - Computational Epidemiology/Project/')
# 
# library(dplyr)
# library(lme4)
# library(parallel)
# 
# covars <- readRDS('Data/Covariates.rds')
# methylation <- readRDS('Data/methylation.rds')     # Very large
# 
# nchunks <- 10
# no_jobs <- 20 

# =========================== #

# Server setup
work_dir<-"~/ce_test/Project/"
setwd(work_dir)
library(lme4,lib.loc = "~/anaconda3/lib/R/library")
library(parallel,lib.loc = "~/anaconda3/lib/R/library")
library(dplyr,lib.loc = "~/anaconda3/lib/R/library")

# Load data
covars <- readRDS('Data/Covariates.rds')
methylation <- readRDS('Data/Methylation.rds')     # Very large

args = commandArgs(trailingOnly=TRUE)
nchunks = as.numeric(args[1])
node = as.numeric(args[2])


no_jobs <- 25     # NOTE: THIS MUST BE THE SAME AS THE JOB NUMBER IN THE PBS SCRIPT

# =========================== #


###### Data preparation ####################

# Remove variables with too much missing data
threshold = 0.95
missing_list <- NULL
for (i in 1:ncol(methylation)) {
  if ( (sum(is.na(methylation[,i]))/length(methylation[,i])) > threshold) {
    missing_list <- c(missing_list, i)
  }
}
if(!is.null(missing_list)) {
  methylation <- methylation[, -missing_list]}


# # ------------------
# # SUBSETTING DATA. FILE TOO LARGE. COMMENT OUT FOR HPC
# s = sample(ncol(methylation), size = 1000)
# methylation <- methylation[,s]
# # ------------------

# Match OMICs data with covariates
index_methyl <- match(rownames(methylation), covars$subjectidp)
covars_methyl <- covars[index_methyl,]
# all(rownames(covars_methyl) == rownames(methylation))   # Quick check that they match


###### Parallelisation preparation ####################

# Tolerance and iterations
tol = 10^(-4)
iterations = 500

# Split columns (biomarkers) using nchunks
ids = as.character(cut(1:ncol(methylation), breaks = nchunks, labels = 1:nchunks))

# Functions to be parallelised (foo1: observed ICCs; foo2: null ICC distribution)
foo1=function(X) {
  
  model1 = lmer(methylation[, column] ~  (1 | id) + (1 | chip) + (1 | position) + age + gender + bmi, data = covars_methyl)
  vcov = as.data.frame(VarCorr(model1))$vcov
  names(vcov) = as.data.frame(VarCorr(model1))[, 1]
  ## For chip > tol, position > tol
  
  # Re-run for singular fits (locating and removing source of singular fit)
  if (vcov[2] < tol) {
    
    if(vcov[3] < tol) {
      model1 = lmer(methylation[, column] ~  (1 | id) + age + gender + bmi, data = covars_methyl)
      vcov = as.data.frame(VarCorr(model1))$vcov
      names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      ## For chip < tol, position < tol
    }
    
    else {
      model1 = lmer(methylation[, column] ~  (1 | id) + (1 | position) + age + gender + bmi, data = covars_methyl)
      vcov = as.data.frame(VarCorr(model1))$vcov
      names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      ## For chip < tol, position > tol
    }
  }
  
  else {
    if (vcov[3] < tol) {
      model1 = lmer(methylation[, column] ~  (1 | id) + (1 | chip) + age + gender + bmi, data = covars_methyl)
      vcov = as.data.frame(VarCorr(model1))$vcov
      names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      ## For chip > tol, position < tol
    }
  }
  
  return(vcov[1]/sum(vcov, na.rm = TRUE))
}


foo2=function(X,iterations) {
  
  tmp <- NULL     # Store vector of iterations
  
  for (i in 1:iterations) {
    
    permuted <- sample(covars_methyl$id)
    covars_methyl$permuted <- permuted
    
    model1 = lmer(methylation[, column] ~  (1 | permuted) + (1 | chip) + (1 | position) + age + gender + bmi, data = covars_methyl)
    vcov = as.data.frame(VarCorr(model1))$vcov
    names(vcov) = as.data.frame(VarCorr(model1))[, 1]
    ## For chip > tol, position > tol
    
    # Re-run for singular fits (locating and removing source of singular fit)
    if (vcov[2] < tol) {
      
      if(vcov[3] < tol) {
        model1 = lmer(methylation[, column] ~  (1 | permuted) + age + gender + bmi, data = covars_methyl)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        ## For chip < tol, position < tol
      }
      
      else {
        model1 = lmer(methylation[, column] ~  (1 | permuted) + (1 | position) + age + gender + bmi, data = covars_methyl)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        ## For chip < tol, position > tol
      }
    }
    
    else {
      if (vcov[3] < tol) {
        model1 = lmer(methylation[, column] ~  (1 | permuted) + (1 | chip) + age + gender + bmi, data = covars_methyl)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        ## For chip > tol, position < tol
      }
    }
    
    
    # Collect matrix of ICCs by biomarker for each iteration
    tmp <- c(tmp, vcov[1]/sum(vcov, na.rm = TRUE))
  }
  return(tmp)
}


###### Parallelisation ######

# Set the number of iterations for a single job to run by splitting the total number of iterations into chunks
iter.chunks = as.character(cut(1:iterations, breaks = no_jobs, labels = 1:no_jobs))
iter = length(iter.chunks[iter.chunks == node])     # Number of iterations for each chunk

# Run parallelisation
t0=Sys.time()
no_cores=detectCores()
cl <- makeCluster(no_cores-1) 
clusterExport(cl, c("methylation", "covars", "covars_methyl", "foo1", "foo2", "ids", "node", "tol", "iterations", "iter"))
clusterEvalQ(cl, library(lme4, dplyr))


icc_methyl=as.matrix(unlist(parLapply(cl=cl,1:nchunks,fun=function(k) {
  X_chunk=methylation[,ids==k, drop=FALSE]
  return(apply(X_chunk, 2, FUN= foo1))
})))

null_icc_matrix=parLapply(cl=cl,1:nchunks, fun=function(k) {
  X_chunk = methylation[,ids==k, drop=FALSE]
  return(apply(X_chunk, 2, FUN =function(x) {foo2(x,iter)}))
})


stopCluster(cl)
t1=Sys.time()
print(t1-t0)

# Reformat nested lists into matrix
null_icc_matrix=t(do.call(cbind,null_icc_matrix))
null_icc_matrix=data.matrix(null_icc_matrix)
colnames(null_icc_matrix) <- NULL


###### Save observed ICCs and array jobs ######
ifelse(dir.exists("Array_jobs"),"",dir.create("Array_jobs", showWarnings = FALSE))
saveRDS(icc_methyl, file = 'icc_methyl.rds')
saveRDS(null_icc_matrix, file = paste0('Array_jobs/null_icc_matrix','_methyl_', node,'.rds'))

# NOTE: P-VALUES NEED TO BE EVALUATED IN A SEPARATE SCRIPT

