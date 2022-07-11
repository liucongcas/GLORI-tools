Glori-tools currently works, but is still being optimized for a better user experience

一： required
STAR > v2.7.5c
bowtie (bowtie1) version>v1.3.0
samtools > v1.10
python > v3.8.3
python package: 
pysam,pandas,argparse,time,collections,os,sys,re,subprocess,multiprocessing,numpy,scipy,math,sqlite3,Bio,statsmodels


二： processing:

1 download fastafiles including genome file and transcriptom file and the gtf annotation file


wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001405.39_GRCh38.p13/GCF_000001635.26_GRCm38.p6_assembly_report.txt
wget https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/annotation_releases/108/GCF_000001405.39_GRCh38.p13/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/annotation_releases/108/GCF_000001405.39_GRCh38.p13/GCF_000001635.26_GRCm38.p6_rna.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/annotation_releases/108/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz

2 change transcriptom fastafile: select longest transform for each gene 


1)python ../change_UCSCgtf.py -i $gtf -j GCF_000001405.39_GRCh38.p13_assembly_report.txt -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens
2)python ../gtf2anno.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens_change2Ens.tbl
3)awk '$3!~/_/&&$3!="na"' GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens_change2Ens.tbl | sed '/unknown_transcript_1/d'  > GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens_change2Ens.tbl2
4)python selected_longest_transcrpts_fa.py -anno GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens_change2Ens.tbl2 -fafile GCF_000001405.39_GRCh38.p13_rna.fa --outname_prx GCF_000001405.39_GRCh38.p13_rna2.fa

3 building index (required)

1) building genome index using STAR
python build_genome_index.py -f ../genome_humanhg38/GRh38_only.fa -pre GRh38_only.fa
you will get:

GRh38_only.fa.rvsCom.fa
GRh38_only.fa.AG_conversion.fa
and the corresbonding index from STAR

2) building transcriptom index using bowtie

python build_transcriptome_index.py -f ../genome_humanhg38/GCF_000001405.39_GRCh38.p13_rna2.fa -pre GCF_000001405.39_GRCh38.p13_rna2.fa

4 get_base annotation (optional)


1) python ../anno_to_base.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens_change2Ens.tbl2 -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens_change2Ens.tbl2.baseanno

2) python ../gtf2genelist.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens -f GCF_000001405.39_GRCh38.p13_rna.fa -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens_change2Ens.genelist > output2

3) awk '$6!~/_/&&$6!="na"' GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens_change2Ens.genelist > GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens_change2Ens.genelist2

4)python ../anno_to_base_remove_redundance_v1.0.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens_change2Ens.tbl2.baseanno -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens_change2Ens.tbl2.noredundance.base -g GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens_change2Ens.genelist2

5 mapping and callsites (required)


Thread=1
genomdir=your_dir
genome=${genomdir}/GRh38_only.fa.AG_conversion.fa
genome2=${genomdir}/GRh38_only.fa.AG_conversion.fa
rvsgenome=${genomdir}/GRh38_only_revCom_2.fa
TfGenome=${genomdir}/GCF_000001405.39_GRCh38.p13_rna2.fa.AG_conversion.fa

annodir=your_dir
baseanno=${annodir}/GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.noredundance.base
anno=${annodir}/GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2
outputdir=your_dir
tooldir=/yourdir/NS-seq-tools
filedir=your_dir

prx=your_prefix
file=your_trimmed reads

1) call m6A sites annotated the with genes,

python ${tooldir}/run_GLORI.py -i $tooldir -q ${file} -T $Thread -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a $anno -b $baseanno -pre ${prx3} -o $outputdir --combine --rvs_fac

2) call m6A sites without annotated genes, in this situation,the background for each m6A sites are the overall conversion rate:

python ${tooldir}/run_GLORI.py -i $tooldir -q ${file} -T $Thread -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -pre ${prx3} -o $outputdir --combine --rvs_fac

The site list obtained by the above two methods is basically the similiar, and there may be a few differential sites in the list.

3) mapping with samples without GLORI treatment

python ${tooldir}/run_GLORI.py -i $tooldir -q ${file} -T $Thread -f ${genome} -rvs ${rvsgenome} -Tf ${TfGenome} -a $anno -b     $baseanno -pre ${prx3} -o $outputdir --combine --untreated


6 resultes:

1 ${your_prefix}_merged.sorted.bam
Overall mapping reads
2 ${your_prefix}_referbase.mpi
miplup files
3 ${your_prefix}.totalCR.txt
txt file for the overall conversion rate
4 finally sites files
${your_prefix}.totalm6A.FDR.csv
Chr: chromosome
Sites: genomic loci
Strand: strand
Gene: annotated gene
CR: conversion rate for genes
AGcov: reads coverage with A and G
Acov: reads coverage with A
Genecov: mean coverage for the whole gene
Ratio: A rate for the sites/ or methylation level for the sites
Pvalue: test for A rate based on the background
P_adjust: FDR ajusted P value
