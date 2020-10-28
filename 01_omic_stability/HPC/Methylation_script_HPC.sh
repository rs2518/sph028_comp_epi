#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=24:mem=48gb
#PBS -N NEW_CODE
#PBS -J 1-25

cd ~/ce_test/Project
module load anaconda3/personal

node=$PBS_ARRAY_INDEX # model index
nchunks=24

Rscript Methylation_script_HPC.R $nchunks $node
