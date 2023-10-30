#!/usr/bin/env bash

trait=$1
workPATH="PATH TO WORKING DIRECTORY"
# SNPSIZE file should contain a row for the number of SNPs per block
# initial tests used: 10, 23, 50, 100, 227, 500, 1000, 2273, 5000, 10000, 22726 snps per block
mapfile -t S < ${workPATH}SNPSIZE-FILENAME

#memory in GB. These amounts were chosen based on the above default SNP sizes
mems=(4 4 4 8 8 16 16 32 64 128 256)
#run time limits. These correspond to empirical testing for the default SNP sizes
tims=(00:59:00 00:59:00 00:59:00 00:59:00 00:59:00 00:59:00 00:59:00 00:59:00 02:59:00 5:59:00 11:59:00)
#loop over the the snp sizes and loop over the blocks (in this case the 22 autosomes) and submit a lasso job for each
for s in {0..10};do
    for c in {1..22}; do
        sbatch --mem=${mems[$s]}G --time=${tims[$s]} lasso.sh $trait $m ${S[$s]} $c $gwas
    done
done
wait