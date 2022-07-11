"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: used for annotation files """


import subprocess
import os
import sys
from optparse import OptionParser
from collections import defaultdict
import re
import argparse
import time


def get_change(changefile):
    d2 = [r1.strip().split("\t") for r1 in open(changefile).readlines()]
    dict_change = {}
    for yr in d2:
        if yr[0][0] != "#":
            ncbi = yr[6]
            ucsc = yr[9]
            dict_change[ncbi] =ucsc
    return dict_change


def change_toanno(file,Output_file,dict_change):
    # print(dict_change)
    d1 = [r1.strip().split("\t") for r1 in open(file).readlines()]
    for line in d1:
        if line[0][0] !="#":
            # print(line)
            row2 = [i.split(" ") for i in line if re.search('gene_id', i)][0]
            gene_id = row2[row2.index("gene_id") + 1].replace('\"', '').strip(";")
            if gene_id not in dict_trans:
                dict_trans[gene_id] = {}
            chr = line[0]
            if chr in dict_change:
                if line[2] == "gene":
                    gene_biotype = ' gene_biotype \"'+row2[row2.index("gene_biotype") + 1].replace('\"', '').strip(";")+"\""
                    dict_gene_biotype[gene_id] = gene_biotype
                elif line[2] == "exon":
                    transcript_id = row2[row2.index("transcript_id") + 1].replace('\"', '').strip(";")
                    if transcript_id in dict_trans[gene_id]:
                        dict_trans[gene_id][transcript_id] += [int(line[3]),int(line[4])]
                    else:
                        dict_trans[gene_id][transcript_id] = [int(line[3]),int(line[4])]
    trans_old = 'la'
    Output = open(Output_file, 'w+')
    for x in range(len(d1)):
        xr = d1[x]
        if xr[0][0] != "#":
            row2 = [i.split(" ") for i in xr if re.search('gene_id', i)][0]
            gene_id = row2[row2.index("gene_id") + 1].replace('\"', '').strip(";")
            if xr[0] in dict_change:
                xr[0] = dict_change[xr[0]]
                if xr[2] == "gene":
                    Output.write("\t".join(xr)+"\n")
                else:
                    trans_new = row2[row2.index("transcript_id") + 1].replace('\"', '').strip(";")
                    # print(trans_new,"*******************")
                    if trans_new != trans_old and not re.search('unknown_transcript',trans_new):
                        # print(" ".join(row2[:-2]),dict_gene_biotype[gene_id])
                        trans_newline = xr[:2] + ['transcript', min(dict_trans[gene_id][trans_new]),
                                                  max(dict_trans[gene_id][trans_new])] + xr[5:8] + [" ".join(row2[:-2]) + dict_gene_biotype[gene_id]]
                        Output.write("\t".join(map(str,trans_newline)) + "\n")
                        x_add = xr[-1] + dict_gene_biotype[gene_id]
                        xr[-1] = x_add
                        Output.write("\t".join(xr) + "\n")
                        trans_old = trans_new
                    else:
                        x_add = xr[-1] + dict_gene_biotype[gene_id]
                        xr[-1] = x_add
                        Output.write("\t".join(map(str,xr)) + "\n")
                        trans_old = trans_new
    Output.close()


if __name__ == "__main__":
    # Parser
    usage="Usage: python %prog -i <gtf> > output.anno"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("-i", "--input",nargs="?", type=str, default=sys.stdin, required=True, help="NCBIgtf file")
    parser.add_argument("-j", "--changefile",nargs="?", type=str, default=sys.stdin, required=True, help="chrchange file")
    parser.add_argument("-o", "--output", nargs="?", type=str, default=sys.stdin, required=True)
    options = parser.parse_args()
    global dict_trans, dict_gene_biotype
    dict_trans = defaultdict(dict)
    dict_gene_biotype = defaultdict(dict)
    dict_change = get_change(options.changefile)
    change_toanno(options.input, options.output,dict_change)
