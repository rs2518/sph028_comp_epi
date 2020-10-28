# Server setup
# work_dir<-"/rds/general/user/mmk4218/home/Test"
# setwd(work_dir)


# no_jobs <- 20      # Check this is the same as the HPC script and PBS script

# Read observed ICCs
icc_trans <- readRDS('Array_jobs/icc_trans.rds')

# Loop results from all nodes and bind into one matrix
null_icc_matrix <- NULL
for (node in (1:no_jobs)) {
  tmp_data <- readRDS(paste0('Array_jobs/null_icc_matrix','_trans_', node,'.rds'))
  null_icc_matrix <- cbind(null_icc_matrix, tmp_data)
}

###### Calculate p-value for observed ICCs ######

# Rank null ICC and obtain p-values for observed ICCs
stab_trans <- rep(0, length(icc_trans))
names(stab_trans) <- dimnames(null_icc_matrix)[[1]]

for (column in 1:length(icc_trans)) {
  stab_trans[column] = sum(null_icc_matrix[column,] >= icc_trans[column])/length(null_icc_matrix[column,])
  
}

###### Save p-values ######

# Set directory to location where 'Results' folder can be found
# dir <- "< PUT DIRECTORY HERE >"
# setwd(dir)

ifelse(dir.exists("Results"),"",dir.create("Results", showWarnings = FALSE))
saveRDS(stab_trans, file = 'Results/stab_trans.rds')



