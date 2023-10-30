#!/usr/bin/env bash
#SBATCH -n 2 -c 20
#SBATCH --job-name=do-lasso
#SBATCH --output=%x-%a-%A.SLURMout
#SBATCH -a 1-5

#The above commands (i.e., the "#SBATCH ..." lines) reflect that initial development was done on a server with the SLURM workload manager: https://slurm.schedmd.com/
#The included commands ask for 
#SBATCH -n 2 -c 20  -> requesting 2 nodes with 20 cores
#SBATCH --job-name=do-lasso  -> defining the output file name head
#SBATCH --output=%x-%a-%A.SLURMout --> defining the rest of the output file format. In this case it corresponds to "file head name"-"array ID"-"job ID".SLURMout
#SBATCH -a 1-5 -> defining slurm array variables (used here for cross-validation)

# because this code is developed in python it assumes that you have a specific python virutal environment that you want to use.
# original development used GCC 8.3.0, OpenMPI 3.1.4, and the intel math kernel library (imkl) 2019.5.281, as well as python 3.8.3
# pysnptools (https://github.com/fastlmm/PySnpTools) is used to load variant level genetic data 
module purge
module load GCC/8.3.0 OpenMPI/3.1.4 imkl/2019.5.281
module load Python/3.8.3
source 'FULLPATH-TO-PYSNPTOOLS/pysnptools/bin/activate'
#slurm array index for performing cross validation
k=$SLURM_ARRAY_TASK_ID
echo $k

#gather input variables: trait name, number of SNPs per block, and block label
traitname=$1
snps=$3
block=$4

OUTDIR='FULLPATH-TO-OUTPATH'/$traitname/
#make output directory if it does not exist
mkdir -p $OUTDIR

genoPATH='FULLPATH-TO-BED-MATRIX'

ML='LASSO'
echo $ML

python ml-single.py --geno-path $genoPATH \
    --trait $traitname \
    --cv-fold $k \
    --ml-type $ML\
    --working-path $OUTDIR \
    --snp-size $snps \
    --blk $block 
