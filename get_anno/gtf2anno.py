"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: used for annotation files """


import re
import argparse

def get_geneanno(list_gene):
    lw = open(options.output, 'w+')
    list_tt=[]
    with open(options.input, 'r') as gtf:
        line = gtf.readline()
        while (line):
            row = line.strip().split("\t")
            if line.startswith("#") or row[2] != 'gene':
                line = gtf.readline()
                continue
            chr = row[0]
            dir = row[6]
            # print(row)
            start, end = int(row[3]) - 1, int(row[4])
            # print(list_gene)
            try:
                row2 = [i.split(" ") for i in row if re.search('gene_id', i)][0]
                gene_id = row2[row2.index("gene_id") + 1].replace('\"', '').strip(";")
                if gene_id not in list_gene:
                    trans_id = gene_id
                    gene_biotype = row2[row2.index("gene_biotype") + 1].replace('\"', '').strip(";")
                    # print(gene_id,gene_biotype)
                    if not trans_id in dictInfo:
                        dictInfo[trans_id] = {'start': [],
                                              'end': [],
                                              'dir': "",
                                              'chr': "",
                                              'gene_id': ""}
                        dictInfo[trans_id]['dir'] = dir
                        dictInfo[trans_id]['chr'] = chr
                        dictInfo[trans_id]['gene_id'] = gene_id
                        dictInfo[trans_id]['start'] = dictInfo[trans_id]['start'] + [start]
                        dictInfo[trans_id]['end'] = dictInfo[trans_id]['end'] + [end]
                    line = gtf.readline()
                else:
                    line = gtf.readline()
            except ValueError:
                line = gtf.readline()
        for trans_id in dictInfo.keys():
            exonCount = str(len(dictInfo[trans_id]['start']))
            exonStarts = ','.join([str(i) for i in sorted(dictInfo[trans_id]['start'])]) + ","
            exonEnds = ','.join([str(i) for i in sorted(dictInfo[trans_id]['end'])]) + ","
            list_tt.append("\t".join(
                ['.', trans_id, dictInfo[trans_id]['chr'], dictInfo[trans_id]['dir'], 'txStart', 'txEnd', 'cdsStart',
                 'cdsEnd', exonCount, exonStarts, exonEnds, ".", dictInfo[trans_id]['gene_id'], '.', '.', '.']))
            # print(len(dictInfo),list_tt)
        lw.writelines("\n".join(list_tt)+"\n")


if __name__ == "__main__":
    # Parser
    usage = "Usage: python %prog -i <gtf> > output.anno"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("-i", dest="input", help="Input file")
    parser.add_argument("-o", dest="output", help="Output file")
    options = parser.parse_args()
    # format:
    # transcript_id gene_id gene_name chr dir 5_UTR_len CDS_len 3_UTR_len 5_UTR start_codon CDS stop_codon 3_UTR 5_len CDS_len 3_len Exon
    dictInfo = {}
    list_gene=[]
    with open(options.input, 'r') as gtf:
        line = gtf.readline()
        while (line):
            row = line.strip().split("\t")
            if line.startswith("#") or row[2] == 'gene':
                line = gtf.readline()
                continue
            chr = row[0]
            dir = row[6]
            # print(row)
            start, end = int(row[3]) - 1, int(row[4])
            try:
                row2 = [i.split(" ") for i in row if re.search('transcript_id', i)][0]
                trans_id = row2[row2.index("transcript_id") + 1].replace('\"', '').strip(";")
                gene_id = row2[row2.index("gene_id") + 1].replace('\"', '').strip(";")
                try:
                    gene_name = row2[row2.index("gene") + 1].replace('\"', '').strip(";")
                except ValueError:
                    gene_name = gene_id
                type = row[2]
                if not trans_id in dictInfo:
                    dictInfo[trans_id] = {'start': [],
                                          'end': [],
                                          'dir': "",
                                          'chr': "",
                                          'gene_id': ""}
                    dictInfo[trans_id]['dir'] = dir
                    dictInfo[trans_id]['chr'] = chr
                    dictInfo[trans_id]['gene_id'] = gene_id
                if type == "exon" or type == "tRNAscan":
                    dictInfo[trans_id]['start'] = dictInfo[trans_id]['start'] + [start]
                    dictInfo[trans_id]['end'] = dictInfo[trans_id]['end'] + [end]
                line = gtf.readline()
            except ValueError:
                line = gtf.readline()
    for trans_id in dictInfo.keys():
        gene_name2=dictInfo[trans_id]['gene_id']
        if gene_name2 not in list_gene:
            list_gene.append(gene_name2)
    get_geneanno(list_gene)
