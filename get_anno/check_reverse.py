import os, sys, time, argparse, pysam
from Bio.Seq import reverse_complement

def checkTrans(genome,site_loci,bin,chr,seq):
    fp = pysam.FastaFile(genome)
    sites_base = fp.fetch(reference=chr, start=site_loci-1, end=site_loci+bin-1).upper()
    sites_base2=reverse_complement(sites_base).replace('A','G')
    X2=reverse_complement(seq)
    return sites_base2==X2


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get mutation rates for m6A_new")
    parser.add_argument("-genome", "--genome", nargs="?", type=str, default=sys.stdin, help="genome")
    parser.add_argument("-site_loci","--site_loci", nargs="?", type=int, default=sys.stdin, help="site_loci")
    parser.add_argument("-bin", "--bin", nargs="?", type=int, default=sys.stdin, help="bin")
    parser.add_argument("-chr", "--chr", nargs="?", type=str, default=sys.stdin, help="chr")
    parser.add_argument("-seq", "--seq", nargs="?", type=str, default=sys.stdin, help="seq")
    args = parser.parse_args()

    result = checkTrans(args.genome,args.site_loci,args.bin,args.chr,args.seq)
    print(result)