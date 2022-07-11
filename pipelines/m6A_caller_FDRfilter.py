
"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used for FDR correction"""
"""Input: [.txt]"""



import pandas as pd
import argparse
from statsmodels.sandbox.stats.multicomp import multipletests
import os
import sys
import time
from time import strftime


def get_data(file1,output,adjustP):
    size = os.path.getsize(file1)
    if size != 0:
        df1_t = pd.read_csv(file1, sep='\t', header=None)
        """df1"""
        pvalue1 = df1_t.iloc[:, 14]
        adjustP1 = multipletests(pvalue1, method='fdr_bh')[1]
        x1=pd.DataFrame({'AdjustPvalue':adjustP1})
        df1_t2 = pd.concat([df1_t,x1],axis=1)
        df1 = df1_t2[df1_t2['AdjustPvalue'] < adjustP]
        df1_total_csv = get_cov_ratio(df1, 'm6A.sites')
        df1_total_csv.to_csv(output + '.csv', index=False, sep="\t")
        sys.stderr.write("[%s] FDR filtering finished successfully!\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    else:
        sys.stderr.write("[%s] Empty file!\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())))


def get_cov_ratio(df,sam):
    chr = df.iloc[:,0]
    sites = df.iloc[:,1]
    strand = df.iloc[:,2]
    gene = df.iloc[:,3].fillna('ELSE')
    transcript = df.iloc[:,5].fillna('ELSE')
    nonCR = df.iloc[:,8]
    AGcov = df.iloc[:, 11]
    Acov = df.iloc[:, 12]
    ratio = df.iloc[:,13]
    normedRatio = df.iloc[:, 13]*(1-df.iloc[:, 8])
    Pvalue = df.iloc[:,14]
    Genecov = df.iloc[:, 16]
    p_adjust = df.iloc[:,18]
    pd_t = pd.DataFrame({'Chr':chr,'Sites':sites,'Strand':strand,'Gene':gene,'CR':1-nonCR,'AGcov':AGcov,'Acov':Acov,'Genecov':Genecov,\
     'Ratio':ratio,'Pvalue':Pvalue,'P_adjust':p_adjust})
    return pd_t


if __name__ == "__main__":
    description = """
	"""
    parser = argparse.ArgumentParser(prog="get_commonSites", fromfile_prefix_chars='@', description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)
    # Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i", "--file", dest="file", required=True, help="file1")
    group_required.add_argument("-o", "--output", dest="output", required=True,
                                help="output prefix. file names are [prefix].[A-cutoff].txt")
    group_required.add_argument("-adp", "--adjustpvalue", dest="adjustpvalue", default=0.005, type=float, help="adjustpvalue, default=0.05")

    options = parser.parse_args()

    get_data(options.file,options.output,options.adjustpvalue)




