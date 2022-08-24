
"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used for aligning reads to Genome and trancriptom"""
"""Input: [fastq files]"""

# import mxnet
import sys
import os
import re
import argparse
import subprocess
import itertools
import time
from heapq import merge
import glob
from time import gmtime, strftime
from Bio.Seq import reverse_complement


parser = argparse.ArgumentParser(description = "building index for reference with three bases")
parser.add_argument("-q", "--fastq", nargs="?", type=str, default=sys.stdin, help = "fastqfiles with surfix as _1.fq;_1.fastq;_2.fq;_2.fastq")
parser.add_argument("-p", "--Threads", nargs="?", type=str, default='1', help = "number of alignment threads to launch")
parser.add_argument("-f", "--reference", nargs="?", type=str, default=sys.stdin, help = "indexfile")
parser.add_argument("-rvs", "--rvsref", nargs="?", type=str, default=sys.stdin, help = "transcriptom reference indexfile")
parser.add_argument("-Tf", "--transref", nargs="?", type=str, default=sys.stdin, help = "transcriptom reference indexfile")
parser.add_argument("-t", "--tools", nargs="?", type=str, default=sys.stdin, help="bowtie,STAR")
parser.add_argument("-m", "--mismatch", nargs="?", type=int, default=2, help="mapping mismatch")
parser.add_argument("-F", "--FilterN", nargs="?", type=str, default=0.5, help="MinOverLread")
parser.add_argument("-mulMax", "--mulMax", nargs="?", type=int, default=1, help="suppress all alignments if > <int> exist")
parser.add_argument("--combine", "--combine", help="whether mapping with changed reads",action="store_true")
parser.add_argument("--untreated", "--untreated", help="if the input is untreated",action="store_true")
parser.add_argument("--local", "--local", help="whether use local algrithm",action="store_true")
parser.add_argument("--readsname", "--readsname", help="whether use readsname as changedfastqname",action="store_true")
parser.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default',help = "--outname_prefix")
parser.add_argument("-o", "--outputdir", nargs="?", type=str, default=sys.stdin, help="outputdir")
parser.add_argument("--rvs_fac", "--rvs_fac",help="whether use rvs_fac algrithm", action="store_true")
args = parser.parse_args()


fastq = args.fastq
Threads = args.Threads
reference = args.reference
rvsref = args.rvsref
transref = args.transref
tools = args.tools
global FilterN
FilterN = args.FilterN

mismatch = args.mismatch
mulMax = args.mulMax
outputdir = args.outputdir
outname_prx = args.outname_prefix


re_digits = re.compile(r'(\d+)')

def embedded_numbers(s):
    s2=s.strip().split("\t")
    pieces = re_digits.split(s2[0])
    pieces[1::2] = map(int, pieces[1::2])
    return pieces


def sort_bedfiles(bedfiles,outputfiles):
    prx2 = bedfiles[:-4]
    path = prx2+"_chunk_*.bed"
    chunksize = 5000000
    fid = 1
    lines = []
    with open(bedfiles, 'r') as f_in:
        f_out = open(prx2+'_chunk_{}.bed'.format(fid), 'w')
        for line_num, line in enumerate(f_in, 1):
            lines.append(line)
            if not line_num % chunksize:
                lines = sorted(lines, key=embedded_numbers)
                f_out.writelines(lines)
                f_out.close()
                lines = []
                fid += 1
                f_out = open(prx2+'_chunk_{}.bed'.format(fid), 'w')
        # last chunk
        if lines:
            lines = sorted(lines, key=embedded_numbers)
            f_out.writelines(lines)
            f_out.close()
            lines = []

    chunks = []
    for filename in glob.glob(path):
        chunks += [open(filename, 'r')]

    with open(outputfiles, 'w') as f_out:
        print('merging bedfiles')
        f_out.writelines(merge(*chunks, key=embedded_numbers))
    subprocess.call("rm -f " + prx2 + "_chunk_*.bed", shell=True)


def change_reads(fastq,changename,output_bed,outputdir,change_fac):
    file = open(fastq,'r')
    fac_t=change_fac[1]
    fac_q=change_fac[0]
    if os.path.exists(outputdir):
        pass
    else:
        os.makedirs(outputdir)
    subprocess.call("rm -f " + changename + " 2>/dev/null",shell=True)
    subprocess.call("rm -f " + output_bed + " 2>/dev/null",shell=True)
    file_change = open(changename,'a+')
    file_bed = open(output_bed,'a+')
    fac = True
    while fac:
        fr = list(itertools.islice(file, step))
        list_change = []
        list_bed = []
        if len(fr) != 0:
            list_fr = [lines.strip().split("\t") for lines in fr]
            for x in range(0,len(list_fr),4):
                yr = list_fr[x+1][0].upper()
                reads_name=list_fr[x][0].split(" ")[0]
                A_sites = [m.start() for m in re.finditer('A', yr)]
                if len(A_sites)>=1:
                    A_sites2 = "_".join(map(str,A_sites))
                    list_bed.append([reads_name[1:],A_sites2])
                    list_change += [reads_name,yr.replace(fac_q,fac_t), \
                            list_fr[x+2][0],list_fr[x+3][0]]
                else:
                    list_bed.append([reads_name[1:],'NA'])
                    list_change += [reads_name, yr, list_fr[x + 2][0], list_fr[x + 3][0]]
            file_change.writelines("\n".join(list_change) + "\n")
            list_bed1 = ['\t'.join(map(str, it)) for it in list_bed]
            file_bed.writelines("\n".join(list_bed1) + "\n")
        else:
            file_change.close()
            file_bed.close()
            fac = False
    sort_bedfiles(output_bed, output_bed + "_sorted")
    subprocess.call("rm -f " + output_bed, shell=True)


def mapping_files(tool,fastq,reference,Threads,muta_N,fqname,outputdir,mulMax,flag):
    outputfile = outputdir +fqname+".sam"
    unmapfastq = outputdir +fqname+"_un_2.fq"
    if tool == "bowtie":
        if args.local:
            print("Erro: bowtie do not support --local parameter")
            os._exit(0)
        else:
            para_0 = 'bowtie -k 1 -m '+ str(mulMax)
            para_A = ' -v '+ str(muta_N)
            para_B = ' --best --strata -p ' + Threads
            para_C = ' -x '+ reference +" "+ fastq +' -S ' + outputfile
            para_unmap = ' --un ' + unmapfastq
            para_end = ' 2>' + outputfile +'.output'
            command = para_0+para_A+para_B+para_C+para_unmap+para_end
            print(command)
            subprocess.call(command,shell=True)
    elif tool == "STAR":
        para_0 = "STAR --runThreadN "+ Threads
        para_g = " --genomeDir "+ reference[:-3]
        para_A = " --limitOutSJcollapsed 5000000 "
        para_B = " --outFilterMismatchNmax " + str(muta_N)
        # para_B_2 = " --outFilterMismatchNoverLmax 0.3"
        # para_B_3 = " --outFilterMismatchNoverReadLmax 1"
        para_B_2=''
        para_B_3=' --outFilterScoreMinOverLread '+FilterN+' --outFilterMatchNminOverLread '+FilterN+' --seedSearchStartLmax 30 '# increase overall mapping sensitivity
        para_C = " --outSAMattributes All --outSAMprimaryFlag AllBestScore --outMultimapperOrder Random --outSAMmultNmax 1 "
        para_D = " --outFilterMultimapNmax " + str(mulMax)
        para_E = " --outFileNamePrefix " + outputfile[:-3] + " --readFilesIn " + fastq
        para_unmap = " --outSAMunmapped Within --outReadsUnmapped Fastx"
        line_command = para_0+para_g+para_A+para_B+para_B_2+para_B_3+para_C+para_D+para_E + para_unmap
        print(line_command)
        subprocess.call(line_command, shell=True)
        print("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile[:-3] + 'Aligned.out.sam | samtools sort -n -O SAM > ' + outputfile)
        subprocess.call("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile[:-3] + 'Aligned.out.sam | samtools sort -n -O SAM > ' + outputfile, shell=True)
        subprocess.call("mv " + outputfile[:-3] + 'Unmapped.out.mate1 ' + unmapfastq, shell=True)
        # subprocess.call("rm -f " + outputfile[:-3] + 'Aligned.out.sam', shell=True)
    return outputfile,unmapfastq


def getbamfiles(outputfile,fac,Threads,flag):
    output_bam = outputfile[:-4] + fac
    print("samtools view -F " + flag + " -bS -@ " + Threads+" -h "+outputfile + \
                    " | samtools sort > " + output_bam)
    subprocess.call("samtools view -F " + flag + " -bS -@ " + Threads+" -h "+outputfile + \
                    " | samtools sort > " + output_bam,shell=True)
    subprocess.call("samtools index " + output_bam, shell=True)
    subprocess.call("rm -f " + outputfile, shell=True)
    return output_bam



def multi_sub(string,sitesA,repl):
    string_change = ''
    if sitesA!="NA":
        A_list = map(int, sitesA.split("_"))
        new = []
        for s in string:
            new.append(s)
        for index in A_list:
            new[index] = repl
        string_change=''.join(new)
    else:
        string_change=string
    return string_change


def reverseReads2(outputfile_change,output_bed,reverse_fac,Threads,flag):
    sorted_sam=outputfile_change[:-4] +"_sorted.sam"
    print("samtools view -F " + flag + " -q 255 -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam)
    subprocess.call("samtools view -F " + flag + " -q 255 -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam,shell=True)
    # print("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > "+sorted_sam )
    # subprocess.call("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam,shell=True)
    subprocess.call("rm -f "+outputfile_change,shell=True)
    reverse_sam = outputfile_change[:-4] + "_r.sam"
    f1 = open(sorted_sam,'r')
    f2 = open(output_bed+"_sorted",'r')
    print(sorted_sam,output_bed+"_sorted")
    subprocess.call("rm -f " + reverse_sam + " 2>/dev/null",shell=True)
    file_reverse = open(reverse_sam,'a+')
    fac=True
    index=0
    index2=0
    old_items='la'
    while fac:
        fr1 = list(itertools.islice(f1, step))
        list_reverse = []
        if len(fr1) != 0:
            list_fr = [lines.strip().split("\t") for lines in fr1]
            for items in list_fr:
                index2 += 1
                reads_A = items[0]
                if reads_A[0] != '@' and reads_A != old_items[0]:
                    for row in f2:
                        its = row.strip().split("\t")
                        reads_S = its[0]
                        if reads_A == reads_S:
                            index += 1
                            reverse_1 = reverse_complement(items[9])
                            reverse_reads = multi_sub(reverse_1, its[1], reverse_fac)
                            items[9] = reverse_complement(reverse_reads)
                            list_reverse.append(items)
                            old_items = items
                            break
                elif reads_A[0] == '@':
                    list_reverse.append(items)
                    index += 1
                elif reads_A == old_items[0]:
                    index += 1
                    list_reverse.append(old_items)
            list_reverse1 = ['\t'.join(map(str, it)) for it in list_reverse]
            file_reverse.writelines("\n".join(list_reverse1) + "\n")
        else:
            f1.close()
            f2.close()
            file_reverse.close()
            fac = False
    print("************reversed_reads == mapped reads***************",index,index2)
    subprocess.call("rm -f " + sorted_sam, shell=True)
    return reverse_sam


def reverseReads(outputfile_change,output_bed,reverse_fac,Threads,flag):
    sorted_sam=outputfile_change[:-4] +"_sorted.sam"
    # print("samtools view -F " + flag + " -q 255 -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam)
    # subprocess.call("samtools view -F " + flag + " -q 255 -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam,shell=True)
    print("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam)
    subprocess.call("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam,shell=True)

    # subprocess.call("rm -f " + outputfile_change, shell=True)
    reverse_sam = outputfile_change[:-4] + "_r.sam"
    f1 = open(sorted_sam,'r')
    f2 = open(output_bed+"_sorted",'r')
    print(sorted_sam,output_bed+"_sorted")
    subprocess.call("rm -f " + reverse_sam + " 2>/dev/null",shell=True)
    file_reverse = open(reverse_sam,'a+')
    fac=True
    index=0
    index2=0
    old_items='la'
    while fac:
        fr1 = list(itertools.islice(f1, step))
        list_reverse = []
        if len(fr1) != 0:
            list_fr = [lines.strip().split("\t") for lines in fr1]
            for items in list_fr:
                index2 += 1
                reads_A = items[0]
                if reads_A[0] != '@' and reads_A != old_items[0]:
                    for row in f2:
                        its = row.strip().split("\t")
                        reads_S = its[0]
                        if reads_A == reads_S:
                            index += 1
                            reverse_reads = multi_sub(items[9], its[1], reverse_fac)
                            items[9] = reverse_reads
                            list_reverse.append(items)
                            old_items = items
                            break
                elif reads_A[0] == '@':
                    list_reverse.append(items)
                    index += 1
                elif reads_A == old_items[0]:
                    index += 1
                    list_reverse.append(old_items)
            list_reverse1 = ['\t'.join(map(str, it)) for it in list_reverse]
            file_reverse.writelines("\n".join(list_reverse1) + "\n")
        else:
            f1.close()
            f2.close()
            file_reverse.close()
            fac = False
    print("************reversed_reads == mapped reads***************",index,index2)
    subprocess.call("rm -f " + sorted_sam, shell=True)
    return reverse_sam


if __name__ == "__main__":
    global step
    step = 10000
    global change_fac,fqname2
    change_fac = 'AG'
    if outname_prx != 'default':
        fqname = outname_prx
    else:
        fqname = "_".join(os.path.basename(fastq).split(".")[:-1])
    outputdir2 = outputdir+"/"
    if os.path.exists(outputdir2):
        pass
    else:
        os.makedirs(outputdir2)
    if args.readsname:
        fqname2 = "_".join(os.path.basename(fastq).split(".")[:-1])
    else:
        fqname2= outname_prx

    if args.untreated:
        sys.stderr.write("[%s]untreated...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        outputfile_untreated, unmapfastq = mapping_files(tools, fastq, reference, Threads, mismatch,
                                                        fqname2, outputdir2, mulMax,'4')
        untreated_bam = getbamfiles(outputfile_untreated,"_s.bam",Threads,'4')
        if args.combine:
            sys.stderr.write("[%s]untreated map to transcriptome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            outputfile_untreated_unmap, _, = mapping_files('bowtie', unmapfastq,  transref, Threads,
                                                mismatch,fqname2 + "_un", outputdir2, mulMax,'4')
            bamAG_unmap = getbamfiles(outputfile_untreated_unmap,"_s.bam", Threads,'4')
    else:
        changefastq = outputdir2 + "/" + fqname2 + "_" + change_fac + "changed_2.fq"
        output_bed = outputdir2 + "/" + fqname2 + "_A.bed"
        if os.path.exists(changefastq):
            sys.stderr.write("[%s] Warning：the changed files already exists, please make sure the input file is correct\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            pass
        else:
            sys.stderr.write("[%s] change to A>G...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            change_reads(fastq, changefastq,output_bed,outputdir2, change_fac)
        sys.stderr.write("[%s] map to genome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        outputfile_changeAG, unmapfastq = mapping_files(tools, changefastq, reference, Threads, mismatch,
                                                      fqname2, outputdir2,mulMax,'20')
        outputfile_changeAG = outputdir2 + fqname2 + ".sam"
        reverse_samAG = reverseReads(outputfile_changeAG,output_bed, 'A',Threads,'20')
        reversed_bamAG = getbamfiles(reverse_samAG,'s.bam', Threads,'20')
        if args.rvs_fac:
            sys.stderr.write("[%s]map to reversegenome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            rvsmapfastq = outputdir2 + fqname2 + "_un_2.fq"
            outputfile_changeAG_rvsmap, _, = mapping_files(tools, rvsmapfastq, rvsref, Threads, mismatch,
                                                          fqname2 + "_rvs", outputdir2, mulMax, '4')
            reverse_samAG_rvsmap = reverseReads2(outputfile_changeAG_rvsmap, output_bed, 'A', Threads, '4')
            reversed_bamAG_rvsmap = getbamfiles(reverse_samAG_rvsmap, 's.bam', Threads, '4')
        else:
            pass

        if args.combine:
            unmapfastq = outputdir2 + fqname2 + "_rvs_un_2.fq"
            sys.stderr.write("[%s]map to transcriptome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            outputfile_changeAG_unmap2, _, = mapping_files('bowtie', unmapfastq, transref, Threads, mismatch,
                                                          fqname2 + "_tf", outputdir2, mulMax,'20')
            reverse_samAG_unmap2 = reverseReads(outputfile_changeAG_unmap2, output_bed, 'A', Threads,'20')
            reversed_bamAG_unmap2 = getbamfiles(reverse_samAG_unmap2, 's.bam', Threads,'20')
