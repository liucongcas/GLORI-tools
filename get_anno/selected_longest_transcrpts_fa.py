"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: used for annotation files """

import argparse
import sys
import time
import pandas as pd
from Bio import SeqIO
import subprocess
import re
from collections import defaultdict,OrderedDict


parser = argparse.ArgumentParser(description="parse the mpileup file")
parser.add_argument("-anno", "--annofile", nargs="?", type=str, default=sys.stdin, help="annofile")
parser.add_argument("-fafile", "--fafile", nargs="?", type=str, default=sys.stdin, help="fafile")
parser.add_argument("-outname_prx", "--outname_prx", nargs="?", type=str, default=sys.stdin, help="outname_prx")
args = parser.parse_args()

annofile = args.annofile
fafile = args.fafile
outname_prx = args.outname_prx

def read_anno(fn):
    output = defaultdict(dict)
    with open(fn, 'r') as input:
        line = input.readline()
        while (line):
            line = line.strip().split("\t")
            trans_id = line[1]
            chr = line[2]
            dir = line[3]
            # All to 0-based
            exonStarts = map(lambda x: int(x), line[9].split(",")[:-1])  # Starts are 0-based
            exonEnds = map(lambda x: int(x) - 1, line[10].split(",")[:-1])  # Ends are 1-based
            # print(line)
            # print(exonStarts)
            gene_id = line[12]
            bins = list(zip(exonStarts, exonEnds))
            # print(bins)
            enst_start = 0
            output[trans_id]['dir'] = dir
            output[trans_id]['ensg'] = gene_id
            output[trans_id]['chr'] = chr
            output[trans_id]['introns'] = OrderedDict()
            output[trans_id]['exons'] = OrderedDict()
            last_end = None
            last_start = None
            enst_start = -1  # 0-based
            if dir == "+":
                for start, end in bins:
                    enst_start += 1
                    enst_end = enst_start + end - start
                    output[trans_id]['exons'][(enst_start, enst_end)] = (start, end)
                    enst_start = enst_end
                    last_end = end
            elif dir == "-":
                bins = bins[::-1]
                for start, end in bins:
                    enst_start += 1
                    enst_end = enst_start + end - start
                    output[trans_id]['exons'][(enst_end, enst_start)] = (start, end)  # Noted that '-' strand is reverse
                    enst_start = enst_end
                    last_start = start
            line = input.readline()
    return output

def get_length(annofile):
    annotation = read_anno(annofile)
    list_output = []
    for i in annotation:
        genome_info = annotation.get(i)
        dir = genome_info['dir']
        gene_name = genome_info['ensg']
        Chr = genome_info['chr']
        genome_info_iter = list(genome_info["exons"].items())
        list_maxend = []
        list_sites = []
        for key, values in genome_info_iter:
            list_maxend += [key[0], key[1]]
            list_sites += [values[0], values[1]]
        len_transcript = max(list_maxend)
        min_sites = min(list_sites)
        max_sites = max(list_sites)
        list_output.append("\t".join(map(str, [Chr, min_sites, max_sites, len_transcript, gene_name, dir, i])))
    genefile = annofile.split("_change2Ens.tbl2")[0] + '.translength'
    OUTPUT = open(genefile, 'w+')
    OUTPUT.writelines("\n".join(list_output) + "\n")
    return genefile


def select_anno(fafile,genefile,outname_prx):
    pd_gene = pd.read_csv(genefile, sep="\t", names=['Chr', 'min_sites', 'max_sites', 'Trans_length', 'Gene', 'Strand', 'Transcript'])
    total_gene_list = list(set(pd_gene['Gene'].values.tolist()))
    pd_gene2 = pd_gene[(~pd_gene['Transcript'].str.contains('XR')) & (~pd_gene['Transcript'].str.contains('XM'))]
    pd_data_partial = pd_gene2.sort_values(by=['Chr', 'Gene', 'Trans_length'], ascending=False)
    xx_p = pd_data_partial.drop_duplicates(['Chr', 'Gene'], keep='first')
    trans_list_p = xx_p['Transcript'].values.tolist()
    gene_list_p = xx_p['Gene'].values.tolist()
    pd_gene3 = pd_gene[~pd_gene['Gene'].isin(gene_list_p)]
    pd_data_other = pd_gene3.sort_values(by=['Chr', 'Gene', 'Trans_length'], ascending=False)
    xx_o = pd_data_other.drop_duplicates(['Chr', 'Gene'], keep='first')
    xx_t = pd.concat([xx_p,xx_o])
    trans_list_o = xx_o['Transcript'].values.tolist()
    trans_list = trans_list_o + trans_list_p
    xx_t.to_csv(genefile + ".longest.trans", sep="\t", index=False)
    index = 0
    subprocess.call('rm ' + outname_prx,shell=True)
    changed_refer2 = open(outname_prx, 'a')
    for record in SeqIO.parse(fafile, "fasta"):
        if record.id in trans_list:
            index += 1
            SeqIO.write(record, changed_refer2, "fasta")
    changed_refer2.close()


if __name__ == "__main__":
    t1 = time.time()
    genefile = get_length(annofile)
    select_anno(fafile,genefile,outname_prx)
    t2 = time.time()
    print(t2-t1)
