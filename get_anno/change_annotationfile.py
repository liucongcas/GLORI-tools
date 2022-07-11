import argparse
import sys
import time


parser = argparse.ArgumentParser(description="annotation file")
parser.add_argument("-annotationfile", "--annotationfile", nargs="?", type=str, default=sys.stdin, help="Input nucleotidefile")
args = parser.parse_args()

annotationfile = args.annotationfile

mutation_type=["A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G"]


def filter_ann(annotationfile):
    dict_ann = {}
    list_total=[]
    d1 = [r1.strip().split("\t") for r1 in open(annotationfile).readlines()]
    for xr in d1:
        par=xr[1]
        if par in dict_ann:
            dict_ann[par].append(xr)
        else:
            dict_ann[par]=[xr]
    for k in dict_ann:
        list_v=dict_ann[k]
        # print(list_v)
        if len(list_v)==1:
            list_total.append(list_v[0])
        else:
            list_tmp = []
            for yr in list_v:
                # print(yr)
                if len(yr[2].split("_"))==1:
                    list_tmp.append(yr)
            if len(list_tmp)==1:
                list_total.append(list_tmp[0])
            else:
                for yr2 in list_tmp:
                    if yr2[2] !="chrY":
                        # print(yr2)
                        list_total.append(yr2)
    print(len(list_total))
    lw1 = open(annotationfile+".filterRep", 'w+')
    list_total2=["\t".join(i) for i in list_total]
    lw1.writelines("\n".join(map(str,list_total2)) + "\n")
    lw1.close()


if __name__=="__main__":
    filter_ann(annotationfile)