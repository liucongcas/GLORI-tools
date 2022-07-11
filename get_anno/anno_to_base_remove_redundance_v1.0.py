"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: used for annotation files """



import os,sys
import argparse
from collections import defaultdict
import time


description = """
"""
parser = argparse.ArgumentParser(prog="",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)

#Require
group_require = parser.add_argument_group("Required")
group_require.add_argument("-i","--input",dest="input",required=True,help="anno file")
group_require.add_argument("-o","--output",dest="output",required=True,help="output bed file") #chr start_0 start
group_require.add_argument("-g","--genelist",dest="genelist",required=True,help="genelist, mark the longest transcript")
group_require.add_argument("--score",dest="biotype_score",required=False,help="Score table for biotypes: biotype[tab]score. If not given, use preset")
group_other = parser.add_argument_group("Other")
options = parser.parse_args()


biotype = {}
transLen = {}
geneName = {}
isoformName = {}
GENES = {}

def changeDict(key,gene_id,trans_id):
	GENES[key] = {'gene': gene_id,
			   'trans': trans_id,
			   'biotype': biotype.get(gene_id),
			   'length': transLen.get(trans_id),
			   'geneName': geneName.get(gene_id),
			   'isoformName': isoformName.get(trans_id),
			   'order':order[biotype.get(gene_id)]
			 }


order = {'3prime_overlapping_ncrna' : 1,
		'3prime_overlapping_ncRNA': 1,
		'antisense' : 0,
		'IG_C_gene' : 9,
		'IG_C_pseudogene' : -1,
		'IG_D_gene' : 2,
		'IG_J_gene' : 2,
		'IG_J_pseudogene' : -1,
		'IG_V_gene' : 2,
		'IG_V_pseudogene' : -1,
		'lincRNA' : 1,
		 'lncRNA' :1,
		 'miRNA' : 1,
		'misc_RNA' : -2,
		'Mt_rRNA' : 9,
		'Mt_tRNA' : 9,
		'polymorphic_pseudogene' : -1,
		'processed_pseudogene' : -1,
		'processed_transcript' : 5,
		'protein_coding' : 10,
		'pseudogene' : -1,
		'rRNA' : 9,
		'sense_intronic' : 5,
		'sense_overlapping' : -2,
		'snoRNA' : 9,
		'snRNA' : 9,
		'Spike_in' : 10,
		'TR_C_gene' : 2,
		'TR_D_gene' : 2,
		'IG_D_pseudogene': -1,
		'TR_J_gene' : 2,
		'TR_J_pseudogene' : -1,
		'TR_V_gene' : 2,
		'TR_V_pseudogene' : -1,
		'IG_LV_gene': 2,
		'IG_pseudogene': -1,
		'TEC':-2,
		'unprocessed_pseudogene':-1,
		'transcribed_processed_pseudogene':-1,
		'transcribed_unprocessed_pseudogene':-1,
		'transcribed_pseudogene':-1,
		'transcribed_unitary_pseudogene':-1,
		'scaRNA':9,
		'ribozyme':0,
		'scRNA':9,
		'bidirectional_promoter_lncRNA':0,
		'unitary_pseudogene':-1,
		'macro_lncRNA':0,
		'sRNA':0,
		'guide_RNA':0,
		'pre_miRNA':9,
		'tRNA':9,
		'SRP_RNA':8,
		'ncRNA':8,
		'nontranslating_CDS':8,
		'RNase_MRP_RNA':8,
		'antisense_RNA':7,
		'C_region':2,
		'C_region_pseudogene':-1,
		'D_segment':2,
		'D_segment_pseudogene':-1,
		'J_segment':2,
		'J_segment_pseudogene':-1,
		'ncRNA_pseudogene':-1,
		'other':10,
		'RNase_P_RNA':8,
		'telomerase_RNA':8,
		'vault_RNA':2,
		'V_segment':2,
		'V_segment_pseudogene':-1,
		'Y_RNA':2
}


if options.biotype_score:
	with open(options.biotype_score,'r') as input:
		for line in input.readlines():
			line = line.strip().split("\t")
			order[line[0]] = float(line[1])

with open(options.genelist,'r') as input:
	line = input.readline()
	while (line):
		line = line.strip().split("\t")
		if line[-1] != "None":
			length = int(line[-1])
		else:
			length = 0
		trans = line[2]
		gene = line[3]
		transLen[trans] = length
		biotype[gene] = line[4]
		geneName[gene] = line[1]
		isoformName[trans] = line[0]
		
		line = input.readline()

with open(options.input,'r') as input,open(options.output,'w') as output:
	line = input.readline()
	line = line.strip()
	row = line.split("\t")
	chr = row[0]
	pos_0 = row[1]
	pos_1 = row[2]
	gene_id = row[4]
	type = biotype.get(gene_id) # type =="unknown"
	dir = row[3]
	trans_id = row[5]
	key = (chr,pos_0,dir)
	print(line)
	changeDict(key,gene_id,trans_id)
	LINE = "\t".join([chr,pos_0,pos_1,dir,gene_id,trans_id,geneName.get(gene_id),isoformName.get(trans_id),type])
	tmp = LINE
	line = input.readline()
	while (line):
		line = line.strip()
		row = line.split("\t")
		chr = row[0]
		pos_0 = row[1]
		pos_1 = row[2]
		gene_id = row[4]
		type = biotype.get(gene_id)
		dir = row[3]
		trans_id = row[5]
		key = (chr,pos_0,dir)
		LINE = "\t".join([chr,pos_0,pos_1,dir,gene_id,trans_id,geneName.get(gene_id),isoformName.get(trans_id),type])
		if key not in GENES: #Finished
			output.write(tmp)
			output.write("\n")
			tmp = LINE
			GENES = {}
			changeDict(key,gene_id,trans_id)
		else:
			if GENES[key]['gene'] == gene_id:
				if GENES[key]['length'] > transLen.get(trans_id):  #different transcipt with same gene,select the shortest
					tmp = LINE
					changeDict(key,gene_id,trans_id)
			else:
				if type == "protein_coding":
					if GENES[key]['biotype'] != "protein_coding":
						tmp = LINE
						changeDict(key,gene_id,trans_id)
					else:
						print(line,GENES[key]['length'],transLen.get(trans_id))
						if GENES[key]['length'] < transLen.get(trans_id):
							tmp = LINE
							changeDict(key,gene_id,trans_id)
				else:
					#if order[type] > GENES[key]['order']:
					if order.get(type) > GENES[key]['order']:
						tmp = LINE
						changeDict(key,gene_id,trans_id)
					#elif order[type] == GENES[key]['order']:
					elif order.get(type) == GENES[key]['order']:
						if GENES[key]['length'] < transLen.get(trans_id):
							tmp = LINE
							changeDict(key,gene_id,trans_id)
		line = input.readline()
	output.write(tmp)
	output.write("\n")
