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
# transcripts <- readRDS('Data/Transcripts.rds')     # Very large
# 
# nchunks <- 10
# no_jobs <- 20 

# =========================== #

# Server setup
work_dir<-"/rds/general/user/mmk4218/home/Test"
setwd(work_dir)
library(lme4,lib.loc ="/rds/general/user/mmk4218/home/anaconda3/lib/R/library")
library(parallel,lib.loc ="/rds/general/user/mmk4218/home/anaconda3/lib/R/library")
library(dplyr,lib.loc ="/rds/general/user/mmk4218/home/anaconda3/lib/R/library")

# Load data
covars <- readRDS('Covariates.rds')
transcripts <- readRDS('Transcripts.rds')     # Very large

args=commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
node=as.numeric(args[2])

no_jobs <- 20     # NOTE: THIS MUST BE THE SAME AS THE JOB NUMBER IN THE PBS SCRIPT

# =========================== #


###### Data preparation ####################

# Remove variables with too much missing data
threshold = 0.95
missing_list <- NULL
for (i in 1:ncol(transcripts)) {
  if ( (sum(is.na(transcripts[,i]))/length(transcripts[,i])) > threshold) {
    missing_list <- c(missing_list, i)
  }
}
if(!is.null(missing_list)) {
transcripts <- transcripts[, -missing_list]}


# # ------------------
# # SUBSETTING DATA. FILE TOO LARGE. COMMENT OUT FOR HPC
# s = sample(ncol(transcripts), size = 1000)
# transcripts <- transcripts[,s]
# # ------------------

# Match OMICs data with covariates
index_trans <- match(rownames(transcripts), covars$subjectidp)
covars_trans <- covars[index_trans,]
# all(rownames(covars_trans) == rownames(transcripts))   # Quick check that they match


###### Parallelisation preparation ####################

# Tolerance and iterations
tol = 10^(-4)
iterations = 1000

# Split columns (biomarkers) using nchunks
ids = as.character(cut(1:ncol(transcripts), breaks = nchunks, labels = 1:nchunks))

# Functions to be parallelised (foo1: observed ICCs; foo2: null ICC distribution)
foo1=function(X) {
  
  model1 = lmer(X ~  (1 | id) + (1 | isolation) + (1| labeling) + (1| hybridization) + age + gender + bmi, data = covars_trans)
  vcov = as.data.frame(VarCorr(model1))$vcov
  names(vcov) = as.data.frame(VarCorr(model1))[, 1]
  ## For isolation, labelling, hybridization > tol
  
  # Re-run for singular fits (locating and removing source of singular fit)
  
  if (vcov[2] < tol) {
    
    if(vcov[3] < tol) {
      
      if(vcov[4] < tol) {
        model1 = lmer(X ~  (1 | id) + age + gender + bmi, data = covars_trans)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      }
      ## For isolation, labelling, hybridization < tol
      
      else{
        model1 = lmer(X ~  (1 | id) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      }
      ## For isolation, labelling < tol | hybridization > tol
      
    }
    
    else {
      if(vcov[4] < tol) {
        model1 = lmer(X ~  (1 | id) + (1| labeling) + age + gender + bmi, data = covars_trans)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      }
      
      else {
        model1 = lmer(X ~  (1 | id) + (1| labeling) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      }
    }
  }
  
  
  else {
    if(vcov[3] < tol) {
      
      if(vcov[4] < tol) {
        model1 = lmer(X ~  (1 | id) + (1 | isolation) + age + gender + bmi, data = covars_trans)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      }
      
      else {
        model1 = lmer(X ~  (1 | id) + (1 | isolation) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      }
    }
    
    else {
      if(vcov[4] < tol) {
        model1 = lmer(X ~  (1 | id) + (1 | isolation) + (1| labeling) + age + gender + bmi, data = covars_trans)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      }
      
      else {
        model1 = lmer(X ~  (1 | id) + (1 | isolation) + (1| labeling) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      }
    }
  }
  
  return(vcov[1]/sum(vcov, na.rm = TRUE))
}


foo2=function(X,iterations) {
  
  tmp <- NULL     # Store vector of iterations
  
  for (i in 1:iterations) {
    
    permuted <- sample(covars_trans$id)
    covars_trans$permuted <- permuted
    
    model1 = lmer(X ~  (1 | permuted) + (1 | isolation) + (1| labeling) + (1| hybridization) + age + gender + bmi, data = covars_trans)
    vcov = as.data.frame(VarCorr(model1))$vcov
    names(vcov) = as.data.frame(VarCorr(model1))[, 1]
    ## For isolation, labelling, hybridization > tol
    
    # Re-run for singular fits (locating and removing source of singular fit)
    
    if (vcov[2] < tol) {
      
      if(vcov[3] < tol) {
        
        if(vcov[4] < tol) {
          model1 = lmer(X ~  (1 | permuted) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
        ## For isolation, labelling, hybridization < tol
        
        else{
          model1 = lmer(X ~  (1 | permuted) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
        ## For isolation, labelling < tol | hybridization > tol
        
      } else {
        if(vcov[4] < tol) {
          model1 = lmer(X ~  (1 | permuted) + (1| labeling) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
        
        else {
          model1 = lmer(X ~  (1 | permuted) + (1| labeling) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
      }
    } else {
      if(vcov[3] < tol) {
        
        if(vcov[4] < tol) {
          model1 = lmer(X ~  (1 | permuted) + (1 | isolation) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }  else {
          model1 = lmer(X ~  (1 | permuted) + (1 | isolation) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
      } else {
        if(vcov[4] < tol) {
          model1 = lmer(X ~  (1 | permuted) + (1 | isolation) + (1| labeling) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        } else {
          model1 = lmer(X ~  (1 | permuted) + (1 | isolation) + (1| labeling) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
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
clusterExport(cl, c("transcripts", "covars", "covars_trans", "foo1", "foo2", "ids", "node", "tol", "iterations", "iter"))
clusterEvalQ(cl, library(lme4, dplyr))


icc_trans=as.matrix(unlist(parLapply(cl=cl,1:nchunks,fun=function(k) {
  X_chunk=transcripts[,ids==k, drop=FALSE]
  return(apply(X_chunk, 2, FUN= foo1))
})))

null_icc_matrix=parLapply(cl=cl,1:nchunks, fun=function(k) {
  X_chunk = transcripts[,ids==k, drop=FALSE]
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
saveRDS(icc_trans, file = 'icc_trans.rds')
saveRDS(null_icc_matrix, file = paste0('Array_jobs/null_icc_matrix','_trans_', node,'.rds'))

# NOTE: P-VALUES NEED TO BE EVALUATED IN A SEPARATE SCRIPT

