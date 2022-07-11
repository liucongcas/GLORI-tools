

"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used to merge multiple BAM filtes to one, then sort and index it"""
"""Input: [.bam]"""


import sys,os
import argparse
import pysam
import time
from time import gmtime, strftime


def merge_bam(fin,fout):
	with pysam.AlignmentFile(fn, 'rb') as INPUT: 
		for read in INPUT:
			# name = read.reference_name
			# new_id = hid_dict.get(name)
			# read.reference_id = new_id
			# print(fin,read.reference_id,lift_over[fin][read.reference_id],read)
			# time.sleep(1)
			read.reference_id = lift_over[fin][read.reference_id]
			read.next_reference_id = -1
			read.next_reference_start = 0
			# read.set_tag("TS", "-")
			# if read.next_reference_id is not None:
				# read.next_reference_id = lift_over[fin][read.next_reference_id]
			fout.write(read)

def read_headers(fn,hid,new_header,hid_dict,lift_over):
	lift_over[fn] = {}
	with pysam.AlignmentFile(fn, 'rb') as INPUT:
		n = 0
		for header in INPUT.header["SQ"]:
			if header['SN'] not in hid_dict:
				# print(header['SN'])
				hid_dict[header['SN']] = hid
				new_header['SQ'].append(header)
				hid += 1
			lift_over[fn][n] = hid_dict[header['SN']]
			n += 1
	return hid,new_header,hid_dict,lift_over


if __name__ == "__main__":
	description = """
	"""
	parser = argparse.ArgumentParser(prog="concat_bam",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	# Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("-i","--input",dest="input", nargs='*',required=True,help="Input bam files")
	group_required.add_argument("-o","--output",dest="output",required=True,help="Output bam")
	group_optional = parser.add_argument_group("Optional")
	group_optional.add_argument("--sort",dest="sort",default=False,action="store_true",help="Sort bam (and delete unsort)")
	group_optional.add_argument("--no-del-bam",dest="no_del_bam",default=False,action="store_true",help="Do not del bam file after sorting")
	group_optional.add_argument("--index",dest="index",default=False,action="store_true",help="Index sorted bam")
	group_optional.add_argument("-t","--threads",dest="threads",default=1,type=int,help="Threads for samtools sort, default=1")
	group_optional.add_argument("-m","--memory",dest="memory",default="1G",help="Memory for samtools sort, default=4G")
	options = parser.parse_args()

	hid = 0
	hid_dict = {}
	lift_over = {}
	new_header = {}
	new_header['HD'] = {'SO': 'unsorted', 'VN': '1.0'}
	new_header['SQ'] = []
	
	for fn in options.input:
		hid,new_header,hid_dict,lift_over = read_headers(fn,hid,new_header,hid_dict,lift_over)
	
	with pysam.AlignmentFile(options.output, 'wb', header = new_header) as OUTPUT: 
		for fn in options.input:
			merge_bam(fn,OUTPUT)
	
	if options.sort == True:
		sys.stderr.write("[%s]Sorting bam...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		if options.threads > 1:
			pysam.sort("-@",str(options.threads),"-m",options.memory,"-o", options.output.replace(".bam",".sorted.bam"),options.output)
		else:
			pysam.sort("-m",options.memory,"-o", options.output.replace(".bam",".sorted.bam"), options.output)
		if options.index == True:
			sys.stderr.write("[%s]Indexing bam...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
			pysam.index(options.output.replace(".bam",".sorted.bam"))
		if options.no_del_bam == False:
			os.remove(options.output)
