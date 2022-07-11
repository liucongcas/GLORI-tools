import pandas as pd
import os, sys
import argparse
import time
from statsmodels.sandbox.stats.multicomp import multipletests


def get_data(file1,file2,output,coverage,count,ratio,pvalue,adjustP,signal):
    df1_total = pd.read_csv(file1, sep='\t', header=None)
    df2_total = pd.read_csv(file2, sep='\t', header=None)
    df1_t = df1_total[(df1_total[11] >= coverage) & (df1_total[12] >= count) & (df1_total[13] >= ratio) & (df1_total[14] < pvalue) & (df1_total[15] >= signal)]
    df2_t = df2_total[(df2_total[11] >= coverage) & (df2_total[12] >= count) & (df2_total[13] >= ratio) & (df2_total[14] < pvalue) & (df2_total[15] >= signal)]
    df1_t.reset_index(drop=True, inplace=True)
    df2_t.reset_index(drop=True, inplace=True)
    """df1"""
    pvalue1 = df1_t.iloc[:, 14]
    adjustP1 = multipletests(pvalue1, method='bonferroni')[1]
    x1=pd.DataFrame({'AdjustPvalue':adjustP1})
    df1_t2 = pd.concat([df1_t,x1],axis=1)
    df1 = df1_t2[df1_t2['AdjustPvalue'] < adjustP]
    """df2"""
    pvalue2 = df2_t.iloc[:, 14]
    adjustP2 = multipletests(pvalue2, method='bonferroni')[1]
    x2 = pd.DataFrame({'AdjustPvalue': adjustP2})
    df2_t2 = pd.concat([df2_t, x2], axis=1)
    df2 = df2_t2[df2_t2['AdjustPvalue'] < adjustP]

    """get common"""
    df_common = df1.merge(df2, on=[0, 1, 2],)
    print(df1_t2.columns,df2_t2.columns,df1_t.shape, df1.shape, df2_t.shape, df2.shape,df_common.shape)

    """uniq mutation rate"""
    df1_uniq_total = df1.merge(df2, on=[0, 1, 2], indicator=True, how="left")[lambda x: x._merge == 'left_only'].drop(
        '_merge', 1)
    df2_uniq_total = df2.merge(df1, on=[0, 1, 2], indicator=True, how="left")[lambda x: x._merge == 'left_only'].drop(
        '_merge', 1)
    df1_uniq = df1.merge(df2, on=[0, 1, 2], indicator=True, how="left")[lambda x: x._merge == 'left_only'].drop('_merge', 1).iloc[:, 0:3]
    df2_uniq = df2.merge(df1, on=[0, 1, 2], indicator=True, how="left")[lambda x: x._merge == 'left_only'].drop('_merge', 1).iloc[:, 0:3]
    df2_in1 = df1_total.merge(df2_uniq, on=[0, 1, 2], indicator=True)
    df1_in2 = df2_total.merge(df1_uniq, on=[0, 1, 2], indicator=True)
    return df1,df2,df1_uniq_total,df2_uniq_total,df2_in1,df1_in2,df_common


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
    p_adjust = df.iloc[:,17]
    pd_t = pd.DataFrame({'Chr':chr,'Sites':sites,'Strand':strand,'Gene':gene,'Transcript':transcript,'NonCR':nonCR,'AGcov':AGcov,'Acov':Acov,\
     'Ratio':ratio,'NormeRatio':normedRatio,'Pvalue':Pvalue,'P_adjust':p_adjust,'Sample':[sam for i in range(df.shape[0])]})
    return pd_t


def get_common_cov_mut(output,df1,df2,df1_uniq_total,df2_uniq_total,df2_in1,df1_in2,df_common):
    df_1 = get_cov_ratio(df1, 'm6A_rep1.sites')
    df_1.to_csv(output + '_m6A_rep1.totalsites.csv', index=False,sep="\t")
    df_2 = get_cov_ratio(df2, 'm6A_rep2.sites')
    df_2.to_csv(output + '_m6A_rep2.totalsites.csv', index=False, sep="\t")

    df_1_uniq_total = get_cov_ratio(df1_uniq_total, 'm6A_rep1_uniq.sites')
    df_1_uniq_total.to_csv(output + '_m6A_rep1.uniqSites.csv', index=False, sep="\t")

    df_2_uniq_total = get_cov_ratio(df2_uniq_total, 'm6A_rep2_uniq.sites')
    df_2_uniq_total.to_csv(output + '_m6A_rep2.uniqSites.csv', index=False, sep="\t")

    df_2_in1 = get_cov_ratio(df2_in1, 'm6A_rep2_uniq_in1.sites')
    df_2_in1.to_csv(output + '_m6A_rep2.uniqSites_in1.csv', index=False, sep="\t")

    df_1_in2 = get_cov_ratio(df1_in2, 'm6A_rep1_uniq_in2.sites')
    df_1_in2.to_csv(output + '_m6A_rep1.uniqSites_in2.csv', index=False, sep="\t")

    # print(df_common.columns)
    df_common1 = df_common.iloc[:,:18]
    # print(df_common1.columns)
    df_common2_t1 = df_common.iloc[:, :3]
    df_common2_t2 = df_common.iloc[:, 18:]
    df_common2=pd.concat([df_common2_t1,df_common2_t2],axis=1)
    # print(df_common2.columns)

    df_common11 = get_cov_ratio(df_common1, 'm6A_common_inRep1')
    df_common11.to_csv(output + '_m6A_common_inRep1.csv', index=False, sep="\t")

    df_common22 = get_cov_ratio(df_common2, 'm6A_common_inRep2')
    df_common22.to_csv(output + '_m6A_common_inRep2.csv', index=False, sep="\t")

    pd_data=pd.concat([df_1,df_2,df_1_uniq_total,df_2_uniq_total,df_1_in2,df_1_in2,df_common11,df_common22],axis=0)
    pd_data.to_csv(output + '_m6A_totalCondition.csv', index=False, sep="\t")
    return pd_data


# filedir1='F:/Chengqiyi/m6A-new/improve/'
# file1 = filedir1 + '/m6A-1.total2.sites'
# file2 = filedir1 + '/m6A-2.total2.sites'
# output = filedir1 + '/Data20201218'
# coverage = 10
# count = 3
# ratio = 0.1
# pvalue = 0.05
# adjustp = 0.05
# signal = 0.9

# df1, df2, df1_uniq_total, df2_uniq_total, df2_in1, df1_in2, df_common = get_data(file1,file2,output,coverage,count,\
#                                 ratio,adjustp,pvalue,signal)
#
# get_common_cov_mut(output, df1, df2, df1_uniq_total, df2_uniq_total, df2_in1, df1_in2, df_common)
#


if __name__ == "__main__":
    description = """
	"""
    parser = argparse.ArgumentParser(prog="get_commonSites", fromfile_prefix_chars='@', description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)
    # Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-f1", "--file1", dest="file1", required=True, help="file1")
    group_required.add_argument("-f2", "--file2", dest="file2", required=True, help="file2")
    group_required.add_argument("-o", "--output", dest="output", required=True,
                                help="output prefix. file names are [prefix].[A-cutoff].txt")
    group_required.add_argument("-c", "--coverage", dest="coverage", default=15, type=int, help = "A+G coverage, default=10")
    group_required.add_argument("-C", "--count", dest="count", default=5, type=int,
                            help="A count, below which the site will not count, default=3")
    group_required.add_argument("-r", "--ratio", dest="ratio", default=0.1, type=float, help="m6A level/ratio, default=0.1")
    group_required.add_argument("-p", "--pvalue", dest="pvalue", default=0.005, type=float, help="pvalue, default=0.05")
    group_required.add_argument("-adp", "--adjustpvalue", dest="adjustpvalue", default=0.005, type=float, help="adjustpvalue, default=0.05")
    group_required.add_argument("-s", "--signal", dest="signal", default=0.8, type=float,
                            help="signal ratio, equals coverage(under A-cutoff)/coverage, default=0.9")

    options = parser.parse_args()

    df1, df2, df1_uniq_total, df2_uniq_total, df2_in1, df1_in2, df_common = get_data(options.file1,options.file2,options.output,options.coverage,options.count,\
                                options.ratio,options.pvalue,options.adjustpvalue,options.signal)

    get_common_cov_mut(options.output, df1, df2, df1_uniq_total, df2_uniq_total, df2_in1, df1_in2, df_common)



