
"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used to get mpileup files"""
"""Input: [merged .bam files]"""


import sys
import os
import argparse
import pysam
from Bio.Seq import reverse_complement
from Bio import SeqIO
import multiprocessing
from multiprocessing import Process,Pool
import signal
import time
from time import strftime


def signal_handler(sig,frame):
	pool.terminate()
	sys.exit()


def worker(conting,start,stop):
	global parent_pid,dictRefSeq,options
	pid = os.getpid()
	tmp_file_name = 'tmp_' + str(parent_pid) + '_' + str(pid)
	max_depth = options.max_depth
	if "GL" in conting:
		max_depth = options.rRNA_max_depth
	else:
		max_depth = options.max_depth
	with pysam.AlignmentFile(options.input, 'rb') as samfile, open(tmp_file_name,'a') as tmp_file:
		try:
			for pileupcolumn in samfile.pileup(conting,start=start,stop=stop,max_depth=max_depth,ignore_orphans=False,ignore_overlaps=False,min_base_quality=0,truncate=True):
				chr = pileupcolumn.reference_name
				ref_seq = dictRefSeq.get(pileupcolumn.reference_name)
				ref_base = ref_seq[pileupcolumn.pos].upper()
				pos = pileupcolumn.pos + 1
				PF_positive_base = {}
				PF_negative_base = {}
				PF_positive_qual = {}
				PF_negative_qual = {}
				PF_positive_A = {}
				PF_negative_A = {}
				for pileupread in pileupcolumn.pileups:
					if pileupread.query_position is not None and pileupread.alignment.query_qualities[pileupread.query_position] >= options.qual:
						not_at_end = True
						if pileupread.alignment.is_reverse:
							if options.trim_tail and pileupread.query_position < options.trim_tail:
								not_at_end = False
							if not_at_end == True and options.trim_head and pileupread.alignment.query_length - options.trim_head <= pileupread.query_position:
								not_at_end = False
						else:
							if options.trim_tail and pileupread.alignment.query_length - options.trim_tail < pileupread.query_position:
								not_at_end = False
							if not_at_end == True and pileupread.query_position < options.trim_head:
								not_at_end = False
						if not_at_end:
							"""pileupread.alignment.query_name: reads name;
							   pileupread.query_position: positions in mapped reads(1-len(reads))
							   pileupread.alignment.query_length: mapped length(includ softclip)"""
							query_name,query_base,query_qual = pileupread.alignment.query_name,\
								pileupread.alignment.query_sequence[pileupread.query_position],\
								pileupread.alignment.query_qualities[pileupread.query_position]
							if pileupread.alignment.is_reverse:
								A_counts = pileupread.alignment.query_sequence.count('T')
								query_base = reverse_complement(query_base)
								if PF_negative_base.get(query_name) is None:
									PF_negative_base[query_name] = query_base
									PF_negative_qual[query_name] = query_qual
									PF_negative_A[query_name] = A_counts
								else:
									lastRead = PF_negative_base.get(query_name)
									lastAcount = PF_negative_A.get(query_name)
									if lastRead != query_base:
										if options.omit:
											PF_negative_base.pop(query_name)
											PF_negative_qual.pop(query_name)
											PF_negative_A.pop(query_name)
										else:
											lastQual = PF_negative_qual.get(query_name)
											if lastQual < query_qual:
												PF_negative_base[query_name] = query_base
												PF_negative_qual[query_name] = query_qual
												PF_negative_A[query_name] = A_counts
									elif lastRead == query_base:
										if lastAcount < PF_negative_A.get(query_name):
											PF_negative_base[query_name] = query_base
											PF_negative_qual[query_name] = query_qual
											PF_negative_A[query_name] = A_counts
							else:
								A_counts = pileupread.alignment.query_sequence.count('A')
								if PF_positive_base.get(query_name) is None:
									PF_positive_base[query_name] = query_base
									PF_positive_qual[query_name] = query_qual
									PF_positive_A[query_name] = A_counts  # A counts in reads
								else:
									lastRead = PF_positive_base.get(query_name)
									lastAcount = PF_positive_A.get(query_name)
									if lastRead != query_base:
										if options.omit:
											PF_positive_base.pop(query_name)
											PF_positive_qual.pop(query_name)
											PF_positive_A.pop(query_name)
										else:
											lastQual = PF_positive_qual.get(query_name)
											if lastQual < query_qual:
												PF_positive_base[query_name] = query_base
												PF_positive_qual[query_name] = query_qual
												PF_positive_A[query_name] = A_counts
									elif lastRead == query_base:
										if lastAcount <= PF_positive_A.get(query_name):
											PF_positive_base[query_name] = query_base
											PF_positive_qual[query_name] = query_qual
											PF_positive_A[query_name] = A_counts

				list_tmp_positive = PF_positive_base.values()
				list_tmp_negative = PF_negative_base.values()
				if len(list_tmp_positive) > 0:
					positive_seq = PF_positive_base.values()
					positive_A_seg = [str(i) for i in PF_positive_A.values()]
					tmp_file.write("\t".join([chr,str(pos),'+',ref_base,','.join(positive_seq),','.join(positive_A_seg)]))
					tmp_file.write("\n")
				if len(list_tmp_negative) > 0:
					negative_seq = PF_negative_base.values()
					negative_C_seg = [str(i) for i in PF_negative_A.values()]
					ref_base = reverse_complement(ref_base)
					tmp_file.write("\t".join([chr,str(pos),'-',ref_base,','.join(negative_seq),','.join(negative_C_seg)]))
					tmp_file.write("\n")
		except ValueError:
			print("Conting [%s] does not exist in @SQ header, pass" % conting)


if __name__ == "__main__":
	description = """
	"""
	parser = argparse.ArgumentParser(prog="pileup_genome_multiprocessing",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Required
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("-i",dest="input",required=True,help="Sorted bam and indexed bam input")
	group_required.add_argument("-o",dest="output",default="genome_merge.pileup.txt",help="Output text")
	group_required.add_argument("-f",dest="fasta",nargs="*",required=True,help="fasta, can be multiple (-f a.fa b.fa ...)")
	#Filter
	group_optional = parser.add_argument_group("Optional")
	group_optional.add_argument("-P",dest="process",type=int,default=1,help="Process number, default=1")
	group_optional.add_argument("-q",dest="qual",type=int,default=10,help="Quality filter, default=30")
	group_optional.add_argument("-s",dest="step",type=int,default=100000,help="Steps for pileup, default=100000")
	group_optional.add_argument("-m","--max-depth",dest="max_depth",type=int,default=10000000,help="Max depth for pileup, default=10000")
	group_optional.add_argument("-M","--max-depth-rRNA",dest="rRNA_max_depth",type=int,default=10000,help="Max depth for pileup, default=10000")
	group_optional.add_argument("--omit-confilct",dest="omit",action="store_true",default=False,help="Omit confilct base, default is using higher sequencing quality")
	group_optional.add_argument("--trim-head",dest="trim_head",type=int,default=0,help="Trim the head of reads, default=0")
	group_optional.add_argument("--trim-tail",dest="trim_tail",type=int,default=0,help="Trim the tail of reads, default=0")
	group_other = parser.add_argument_group("Other")
	group_other.add_argument("--version",action="version",version="%(prog)s 1.3")
	options = parser.parse_args()
	
	parent_pid = os.getpid()
	sys.stderr.write("[%s] Pileup genome, processers: %d, pid: %d\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),options.process ,parent_pid))
	t1=time.time()
	dictRefSeq = {}
	RefBins = []
	for genome in options.fasta:
		for seq in SeqIO.parse(genome,"fasta"):
			if seq.id not in dictRefSeq:
				dictRefSeq[str(seq.id)] = seq.seq
				if "GL" in seq.id:
					for start in range(0,len(seq.seq),100):
						RefBins.append((seq.id,start,start+100))
				else:
					for start in range(0,len(seq.seq),options.step):   # split genome into bins
						RefBins.append((seq.id,start,start+options.step))
	index = 1
	signal.signal(signal.SIGINT,signal_handler)
	multiprocessing.freeze_support()
	pool = multiprocessing.Pool(options.process)
	try:
		for item in RefBins:
			chr,start,stop = item
			pool.apply_async(func=worker,args=(chr,start,stop,))
		pool.close()
		print("Merging tmp files222...")
		pool.join()
		print("Merging tmp files333...")
		tmp_names = "tmp_" + str(parent_pid) + "*" 
		tmp_merge_name = options.output
		print("Merging tmp files...")
		sys.stderr.write("[%s] Merging TEMPs\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		os.system("cat %s > %s" % (tmp_names,tmp_merge_name))
		print("Removing pileup tmp files...")
		os.system('rm %s' % tmp_names)
	finally:
		pool.terminate()
	t2 = time.time()
	print("time spended:",t2-t1)
	sys.stderr.write("[%s] Genome pileup finished.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))