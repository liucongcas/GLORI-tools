import argparse
import sys
import pysam

def FilterAll(CutoffAs,referfa,seqment,OUTPUT,Filtered,depth):
    line = seqment.strip().split("\t")
    mapped_base = line[4].split(",")
    A_counts = map(int,line[5].split(","))
    if len(mapped_base) > depth:
        site_chr = line[0].split("_")[0]
        site_loci = int(line[1])
        sites_base = referfa.fetch(reference=site_chr, start=site_loci - 1, end=site_loci).upper()
        if sites_base == "A":
            N90 = False
            index_As = [i for i in range(len(A_counts)) if A_counts[i] >= CutoffAs]
            ratio = len(index_As)/len(A_counts)
            if ratio <=90:
                N90 = True
            list_base=[]
            for j in range(len(mapped_base)):
                if j not in index_As:
                    list_base.append(mapped_base[j])
            A_ratio = list_base.count("A")/len(list_base)
            T_ratio = list_base.count("T")/len(list_base)
            C_ratio = list_base.count("C")/len(list_base)
            G_ratio = list_base.count("G")/len(list_base)
            line += [A_ratio,G_ratio,T_ratio,C_ratio,N90]
            OUTPUT.write("\t".join(line)+"\n")
        else:
            Filtered.write("\t".join(line)+"\n")
    else:
        Filtered.write("\t".join(line)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="filered mpilup")
    group_required = parser.add_argument_group("Required")
    parser.add_argument("-input", "--input", nargs="?", type=str, default='default', help="mpileup file")
    parser.add_argument("-referFa", "--referFa", nargs="?", type=str, default='default', help="reference fa")
    parser.add_argument("-outname_prx", "--outname_prx", nargs="?", type=str, default=sys.stdin, help="outname_prx")

    group_optional = parser.add_argument_group("Optional")
    parser.add_argument("--FilterAll", "--FilterAll", default=True, action="store_true",
                        help="Filter reads with 3A and 3A/N90")
    parser.add_argument("--FilterAs", "--FilterAs", default=False, action="store_true",
                        help="Filter reads with multiA untransformed")
    parser.add_argument("--FilterPolyA", "--FilterPolyA", default=False, action="store_true",
                        help="Filter reads with untransformed A in polyA")
    parser.add_argument("--FilterContinuousA", "--FilterContinuousA", default=False, action="store_true",
                        help="Filter reads with multiA untransformeed")
    parser.add_argument("-CutoffAs", "--CutoffAs", nargs="?", type=int, default=3, help="CutoffAs")
    parser.add_argument("-CutoffPloyA", "--CutoffPloyA", nargs="?", type=int, default=3, help="CutoffPloyA")
    parser.add_argument("-DepthCut", "--DepthCut", nargs="?", type=int, default=20, help="DepthCut")
    parser.add_argument("--outputfiltered", "--outputfiltered", action="store_true", default='True',
                        help="if output filtered alignment")

    options = parser.parse_args()

    referfa = pysam.FastaFile(options.referFa)
    if options.FilterAll:
        if options.outputfiltered:
            with open(options.input) as INPUT, open(options.outname_prx + ".rmAll.txt") as OUTPUT,\
                    open(options.outname_prx+".rmAll.filtered.txt") as Filtered:
                for segment in INPUT:
                    FilterAll(options.CutoffAs,referfa,segment, OUTPUT, Filtered,options.DepthCut)

