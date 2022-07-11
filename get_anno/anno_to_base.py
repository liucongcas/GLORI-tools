
"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: used for annotation files """


import os,sys
import argparse


description = """
"""
parser = argparse.ArgumentParser(prog="",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
#Require
group_require = parser.add_argument_group("Required")
group_require.add_argument("-i","--input",dest="input",required=True,help="anno file")
group_require.add_argument("-o","--output",dest="output",required=True,help="output bed file") #chr start_0 start
#group_require.add_argument("-","--special",dest="special",required=False,help="")
group_other = parser.add_argument_group("Other")
options = parser.parse_args()

with open(options.input,'r') as input,open(options.output,'w') as output:
	line = input.readline()
	while (line):
		line = line.strip().split("\t")
		trans_id = line[1]
		chr = line[2]
		dir = line[3]
		exonCounts = int(line[8])
		exonStarts = line[9].split(",")
		exonEnds = line[10].split(",")
		gene_id = line[12]
		for i in range(exonCounts):
			start = int(exonStarts[i]) #0-based
			end = int(exonEnds[i]) #1-based
			for pos in range(start,end):
				out = "\t".join([chr,str(pos),str(pos+1),dir,gene_id,trans_id])
				output.write(out)
				output.write("\n")
		line = input.readline()
sorted_out = options.output + ".sorted"
os.system("sort -k 4,4 -k 1,2 -S 10G -T ./ --parallel 10 -T ./ %s > %s" % (options.output,sorted_out))
