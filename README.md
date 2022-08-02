# GLORI-tools v1.0
**GLORI-tools currently works, but is still being optimized for a better user experience.**

## Table of content
* Background
* Installation and Requirement
* Example and Usage
* Maintainers and Contributing
* License

## Background
### GLORI
We developed an absolute m6A quantification method (“GLORI”) that is conceptually similar bisulfite sequencing-based quantification of DNA 5-methylcytosine.
GLORI relies on glyoxal and nitrite-mediated deamination of unmethylated adenosines while keeping m6A intact, thereby achieving specific and efficient m6A detection.

### GLORI-tools

GLORI-tools is a bioinformatics pipeline tailored for the analysis of high-throughput sequencing data generated by GLORI.
![GLORI pipeline]( https://github.com/liucongcas/GLORI-tools/blob/main/GLORI-pipeline.jpg " GLORI pipeline ")

## Installation and Requirement
### Installation
GLORI-tools is written in Python3 and is executed from the command line. To install GLORI-tools simply download the code and extract the files into a GLORI-tools installation folder.

### GLORI-tools needs the following tools to be installed and ideally available in the PATH environment:
* STAR ≥ v2.7.5c
* bowtie (bowtie1) version ≥ v1.3.0
* samtools ≥ v1.10
* python ≥ v3.8.3

### GLORI-tools needs the following python package to be installed:
pysam,pandas,argparse,time,collections,os,sys,re,subprocess,multiprocessing,numpy,scipy,math,sqlite3,Bio,statsmodels

## Example and Usage:

### 1. Generated annotation files 
1.1 download files for annotation (required): 
* ``` wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt ```
* ``` wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz ```

1.2 Unify chromosome naming in GTF file and genome file:
 
* ```python ./get_anno/change_UCSCgtf.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf -j GCF_000001405.39_GRCh38.p13_assembly_report.txt -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens ```

### 2. get reference for reads alignment (required)

2.1 download reference genome and transcriptome
* ```wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz```
* ```wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.fna.gz```

2.2 get the longest transcript for genes 

2.2.1 Get the required annotation table files:

* ```python ./get_anno/gtf2anno.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl``` 

* ```awk '$3!~/_/&&$3!="na"' GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl | sed '/unknown_transcript_1/d'  > GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2```

2.2.2 Get the longest transcript:

* ```python ./get_anno/selected_longest_transcrpts_fa.py -anno GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2 -fafile GCF_000001405.39_GRCh38.p13_rna.fa --outname_prx GCF_000001405.39_GRCh38.p13_rna2.fa```

2.3 build genome index using STAR

* ```python ./pipelines/build_genome_index.py -f $genome_fafile -pre $output```

you will get:
* $output.rvsCom.fa
* $output.AG_conversion.fa
* the corresponding index from STAR

2.4 build transcriptome index using bowtie

* ```python ./pipelines/build_transcriptome_index.py -f $ -pre GCF_000001405.39_GRCh38.p13_rna2.fa```

### 3. get_base annotation (optional)

3.1 get annotation at single-nucleotide resolution

* ```python ./get_anno/anno_to_base.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2 -o GCF_000001405.39_GRCh38.p13_genomic.gtf_ _change2Ens.tbl2.baseanno```

3.2 get required annotation file for further removal of duplicated loci

* ``` python ./get_anno/gtf2genelist.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens -f GCF_000001405.39_GRCh38.p13_rna.fa -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist > output2```

* ```awk '$6!~/_/&&$6!="na"' GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist > GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist2```

3.3 Removal of duplicated loci in the annotation file

* ```python ./get_anno/anno_to_base_remove_redundance_v1.0.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.baseanno -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.noredundance.base -g GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist2```


### 4. alignment and call sites (required)
### 4.1 Example shell scripts

|Used files |
| :--- |
| Thread=1 |
| genomdir=your_dir |
| genome=${genomdir}/GRh38_only.fa.AG_conversion.fa |
| genome2=${genomdir}/GRh38_only.fa.AG_conversion.fa |
| rvsgenome=${genomdir}/GRh38_only_revCom_2.fa |
| TfGenome=${genomdir}/GCF_000001405.39_GRCh38.p13_rna2.fa.AG_conversion.fa |
| annodir=your_dir |
| baseanno=${annodir}/GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.noredundance.base |
| anno=${annodir}/GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2 |
| outputdir=your_dir |
| tooldir=/yourdir/NS-seq-tools |
| filedir=your_dir |
| prx=your_prefix |
| file=your_trimmed reads | 

### 4.2 Call m6A sites annotated the with genes
* ``` python ${tooldir}/run_GLORI.py -i $tooldir -q ${file} -T $Thread -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a $anno -b $baseanno -pre ${prx3} -o $outputdir --combine --rvs_fac```

### 4.3 Call m6A sites without annotated genes.

* ``` python ${tooldir}/run_GLORI.py -i $tooldir -q ${file} -T $Thread -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -pre ${prx3} -o $outputdir --combine --rvs_fac ```

* In this situation, the background for each m6A sites are the overall conversion rate.

* The site list obtained by the above two methods is basically the similar, and there may be a few differential sites in the list.

### 4.4 mapping with samples without GLORI treatment

* ```python ${tooldir}/run_GLORI.py -i $tooldir -q ${file} -T $Thread -f ${genome} -rvs ${rvsgenome} -Tf ${TfGenome} -a $anno -b     $baseanno -pre ${prx3} -o $outputdir --combine –untreated ```


### 5 Resultes:
### 5.1 Output files of GLORI

| Output files | Interpretation |
| :---: | :---: |
| ${your_prefix}_merged.sorted.bam | Overall mapping reads |
|  ${your_prefix}_referbase.mpi |  miplup files |
|  ${your_prefix}.totalCR.txt | txt file for the overall conversion rate |
|  finally sites files  | ${your_prefix}.totalm6A.FDR.csv |

### 5.2 GLORI sites files

| Columns | Interpretation |
| :---: | :---: |
| Chr | chromosome |
| Sites | genomic loci |
| Strand | strand |
| Gene | annotated gene |
| CR | conversion rate for genes |
| AGcov | reads coverage with A and G |
| Acov | reads coverage with A |
| Genecov | mean coverage for the whole gene |
| Ratio | A rate for the sites/ or methylation level for the sites |
| Pvalue | test for A rate based on the background |
| P_adjust | FDR ajusted P value |

## Maintainers and Contributing
* GLORI-tools is developed and maintained by Cong Liu (liucong-1112@pku.edu.cn).
* The development of GLORI-tools is inseparable from the open source software RNA-m5C (https://github.com/SYSU-zhanglab/RNA-m5C).

## Licences
* Released under MIT license







