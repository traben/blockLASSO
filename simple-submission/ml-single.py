import numpy as np
import pandas as pd
import scipy as sp
import math as math
import time as time
import argparse
import textwrap
import os
from sklearn import linear_model
import sklearn as skl
from pysnptools.snpreader import Bed
import pipeline_utilities as pu
import multiprocessing as mp



def runML(genPATH,trait,index,MLTYPE,workPATH,snpSIZE,blk,**kwargs):
    print(trait)
    print(workPATH)
    print('cv fold: '+str(index))
    # number of steps in the ML path
    nstep = 190
    print(str(nstep)+' steps in ML path')

    # lamratio defines the ratio of the smallest to largest LASSO penalty coefficient (i.e., it determines the length of the path)
    # typically a 0.001<lamratio<0.01 works for most traits
    lamratio = 0.01
    print('lambda ratio: '+str(lamratio))
    blk = int(blk)
    
    
    print('load paths')
    #input paths

    #gwas file should include at least SNP names and p-values
    gwasPATH = 'FULL-PATH-TO-GWAS'

    #train set should include ids for individuals used in training
    trainPATH = 'FULL-PATH-TO-TRAINING-SET'
    # output paths
    # the scikit-learn lasso_path has 3 main outputs (https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.lasso_path.html)
    lamPATH = 'FULL-PATH-TO-lambda-values-OUTPUT'
    betaPATH = 'FULL-PATH-TO-beta-values-OUTPUT'
    gapPATH = 'FULL-PATH-TO-dualgap-values-OUTPUT'
    
    print('load fam/bim/phen/gwas')    
    G = Bed(genPATH,count_A1=False)
    fam = pd.read_csv(genPATH+".fam",header=None,sep=' ')
    bim = pd.read_csv(genPATH+".bim",header=None,sep='\t')

    # phenotype file must have at least individual ID (IID) and phenotype value (PHENO)
    phen = pd.read_csv(phenPATH,header=None,sep='\s+',names=['FID','IID','PHENO'])
    # GWAS file must have at least SNP, p-value, and block label
    # gwas files used in development came from PLINK output and included:
    # blank,chr, snp, bp, a1, fa, fu, a2 ,x2, P, OR, blank
    gwas = pd.read_csv(gwasPATH,sep='\s+')
    
    top = snpSIZE
    print(f'SNP per block: {top}')
    print('compute subsets')
    # TESTING BLOCKS WERE COMPRISED OF SINGLE CHROMOSOMES. THE FOLLOWING MUST BE EDITED TO INCLUDE OTHER BLOCKS
    sexchr = bim[0].astype(int).eq(blk)
    best=gwas[sexchr].sort_values(by='P',ascending=True)['SNP'][0:top]
    subsetP = bim[1].isin(best)
    subsetP = np.stack(pd.DataFrame(list(range(bim.shape[0])))[subsetP].values,axis=1)[0]

    #load training indeces 
    train = np.loadtxt(trainPATH,dtype=int)
    train_inds = phen['IID'].isin(train.T[0])

    #Loading bed data
    print('load BED data')
    bed_data = Bed(genPATH,count_A1=False)

    #this function is defined in the "pipeline_utilities.py" file
    snpdata = pu.read_bed_file("PATH-TO-BED",
                               phen['IID'],
                               best,
                               snpreader=bed_data,
                               is_sorting_samples=False,
                               is_sorting_snps=True,
                               read_data=True)


    subG=snpdata.val[train_inds]
    target_phen = phen['PHENO'].loc[train_inds].values

                       
    print("Calculate means")
    # calculate column means without including nan values
    center = np.zeros(subG.shape[1])
    spread = np.zeros(subG.shape[1])
    for col in range(0,subG.shape[1]):
        center[col] = np.nanmean(subG[:,col])
        spread[col] = np.nanstd(subG[:,col])

    print("nan replacement")     
    # nan replacement
    missing = np.argwhere(np.isnan(subG))
    for row in range(0,missing.shape[0]):
        ind1 = missing[row,0]
        ind2 = missing[row,1]
        subG[ind1,ind2] = center[ind2]

    print("Standardize")    
    # standardize the columns
    for col in range(0,subG.shape[1]):
        val = spread[col]
        if spread[col] == 0.0:
            val = 1.0
        subG[:,col] = (subG[:,col] - center[col])/val

    y = target_phen
    # standardize the phenotype
    ymu = np.mean(y)
    ysig = np.std(y)
    y = (y-ymu)/ysig
        
    # do the lasso
    print("Begin "+str(MLTYPE),flush=True)    
    t = time.time()
    path = skl.linear_model.lasso_path(subG,y,n_alphas=nstep,eps=lamratio,n_iter=1500)
    elapsed = time.time() - t
    print(str(MLTYPE)+" time:",flush=True)
    print(elapsed)
    
    betas = path[1]
    lamb = path[0]
    gap = path[2]

        
    metadat = bim.iloc[subsetP,:]
    metadat = metadat.reset_index(drop=True)

    #rescale LASSO coefficients to the standardized values
    betas = pd.DataFrame(np.transpose(np.transpose(betas)*np.transpose(ysig/spread)))
    lamb = pd.DataFrame(lamb)
    gap = pd.DataFrame(gap)

    #output results
    gap.to_csv(r''+gapPATH,sep=' ',index=False,header=False)
    out = pd.concat([metadat,pd.DataFrame(center),betas],ignore_index=True,axis=1)
    out.to_csv(r''+betaPATH,sep = ' ',index=False,header=False)
    lamb.to_csv(r''+lamPATH,sep=' ',index=False,header=False)


    return 0





def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     prog='ML',
                                     usage='%(ML)s',
                                     description='''Runs the ML path algo for lasso, enet, and l1-logistic.''')

    # essential arguments
    required_named = parser.add_argument_group('Required named arguments')

    # path to genotype bed matrix
    required_named.add_argument('--geno-path',
                                type=str,
                                required=True,
                                help='path to genotypes')

    # trait name
    required_named.add_argument('--trait',
                                type=str,
                                required=True,
                                help='name of trait')

    # cv fold
    required_named.add_argument('--cv-fold',
                                type=str,
                                required=True,
                                help='index variable, 1-5')

    # ML algorithm
    required_named.add_argument('--ml-type',
                                type=str,
                                required=True,
                                help='start of pheno file name to regress on')

    # working and output path
    required_named.add_argument('--working-path',
                                type=str,
                                required=True,
                                help='Where all the output goes')

    # snp set size per block
    required_named.add_argument('--snp-size',
                                type=int,
                                required=True,
                                help='snp set size')
    # chromosome
    required_named.add_argument('--blk',
                                type=int,
                                required=True,
                                help='block')

    # optional arguments
    # optional_named = parser.add_argument_group('Optional named arguments')
    
    args = parser.parse_args()
    
    

    runML(args.geno_path,args.trait,args.cv_fold,args.gwas_type,args.cov_type,args.ml_type,args.working_path,args.train_size,args.snp_size,args.chrm)


exit(main())
