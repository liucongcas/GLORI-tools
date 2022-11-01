"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used to get referbase for mpileup files with ternary genomes"""
"""Input: [.mpileup]"""


import argparse
import sys
import pysam
import time

def FilterAll(referfa,seqment,OUTPUT):
    line = seqment.strip().split("\t")
    if len(line[5].split(",")) > 1:
        site_chr = line[0].split('_AG_converted')[0]
        site_loci = int(line[1])
        sites_base = referfa.fetch(reference=site_chr, start=site_loci - 1, end=site_loci).upper()
        line[3] = sites_base
        OUTPUT.write("\t".join(map(str,line))+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="filered mpilup")
    # Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-input", "--input", nargs="?", type=str, default='default', help="mpileup file")
    group_required.add_argument("-referFa", "--referFa", nargs="?", type=str, default='default', help="reference fa")
    group_required.add_argument("-outname_prx", "--outname_prx", nargs="?", type=str, default=sys.stdin, help="outname_prx")

    options = parser.parse_args()
    referfa = pysam.FastaFile(options.referFa)
    with open(options.input,'r') as INPUT, open(options.outname_prx + ".referbase.mpi",'w+') as OUTPUT:
        for segment in INPUT:
            FilterAll(referfa,segment, OUTPUT)

