# Server setup
work_dir<-"~/ce_test/Project/"
setwd(work_dir)

no_jobs <- 25      # Check this is the same as the HPC script and PBS script

# Read observed ICCs
icc_methyl <- readRDS('Array_jobs/icc_methyl.rds')

# Loop results from all nodes and bind into one matrix
null_icc_matrix <- NULL
for (node in (1:no_jobs)) {
  tmp_data <- readRDS(paste0('Array_jobs/null_icc_matrix','_methyl_', node,'.rds'))
  null_icc_matrix <- cbind(null_icc_matrix, tmp_data)
}

###### Calculate p-value for observed ICCs ######

# Rank null ICC and obtain p-values for observed ICCs
stab_methyl <- rep(0, length(icc_methyl))
names(stab_methyl) <- dimnames(null_icc_matrix)[[1]]

for (column in 1:length(icc_methyl)) {
  stab_methyl[column] = sum(null_icc_matrix[column,] >= icc_methyl[column])/length(null_icc_matrix[column,])
  
}

###### Save p-values ######

# Set directory to location where 'Results' folder can be found
# dir <- "< PUT DIRECTORY HERE >"
# setwd(dir)

ifelse(dir.exists("Results"),"",dir.create("Results", showWarnings = FALSE))
saveRDS(stab_methyl, file = 'Results/stab_methyl.rds')



