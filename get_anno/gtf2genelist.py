"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: used for annotation files """

import sys
import os
from optparse import OptionParser
from Bio import SeqIO
import re
import time


def get_genelist(list,list_gene):
	with open(options.input,'r') as gtf:
		line = gtf.readline()
		while(line):
			row = line.strip().split("\t")
			type = row[2]
			# print(type)
			if line.startswith("#") or type != "gene":
				line = gtf.readline()
				continue
			# print(line)
			chr = row[0]
			dir = row[6]
			start = row[3]
			end = row[4]
			row2 = [i.split(" ") for i in row if re.search('gene_id', i)][0]
			gene_id = row2[row2.index("gene_id")+1].replace('\"','').strip(";")
			if gene_id not in list_gene:
				try:
					gene_name = row2[row2.index("gene")+1].replace('\"','').strip(";")
				except ValueError:
					gene_name = gene_id
				trans_id = gene_id
				trans_name = gene_id
				try:
					gene_biotype = row2[row2.index("gene_biotype")+1].replace('\"','').strip(";")
				except ValueError:
					gene_biotype = "unknown"
				# length = str(dictLength.get(trans_id))
				length = str(abs(int(end)-int(start)))
				list.append([trans_name,gene_name,trans_id,gene_id,gene_biotype,chr,start,end,dir,length])
				print(gene_id, gene_biotype)
				if gene_id not in list_gene:
					list_gene.append(gene_id)
			line = gtf.readline()
	list_sorted = sorted(list, key=lambda x: (x[0], x[1]))
	list_sorted2 = ["\t".join(map(str,i)) for i in list_sorted]
	# print(list_sorted2)
	lw = open(options.output,'w+')
	lw.writelines("\n".join(list_sorted2)+"\n")


if __name__ == "__main__":
	#Parser
	usage = "USage: python %prog -i <gtf> -f <fasta>"
	parser = OptionParser(usage=usage)
	parser.add_option("-i",dest="input",help="GTF file")
	parser.add_option("-f",dest="fasta",help="Fasta file")
	parser.add_option("-o",dest="output",help="output")
	(options,args) = parser.parse_args()
	dictLength = {}

	list_gene = []
	for seq in SeqIO.parse(options.fasta,"fasta"):
		dictLength[seq.id] = len(seq.seq)
	
	with open(options.input,'r') as gtf:
		line = gtf.readline()
		list = []
		while(line):
			if line.startswith("#"):
				line = gtf.readline()
				continue
			row = line.strip().split("\t")
			type = row[2]
			if type == "transcript":
				chr = row[0]
				dir = row[6]
				start = row[3]
				end = row[4]
				row2 = [i.split(" ") for i in row if re.search('transcript_id', i)][0]
				gene_id = row2[row2.index("gene_id")+1].replace('\"','').strip(";")
				try:
					gene_name = row[row.index("gene")+1].replace('\"','').strip(";")
				except ValueError:
					gene_name = gene_id
				trans_id = row2[row2.index("transcript_id")+1].replace('\"','').strip(";")
				try:
					trans_name = row2[row2.index("transcript_name")+1].replace('\"','').strip(";")
				except ValueError:
					trans_name = trans_id
				try:
					gene_biotype = row2[row2.index("gene_biotype")+1].replace('\"','').strip(";")
				except ValueError:
					gene_biotype = "unknown"
				# length = str(dictLength.get(trans_id))
				length = str(abs(int(end) - int(start)))
				list.append([trans_name,gene_name,trans_id,gene_id,gene_biotype,chr,start,end,dir,length])
				if gene_id not in list_gene:
					list_gene.append(gene_id)
			line = gtf.readline()

	get_genelist(list,list_gene)
