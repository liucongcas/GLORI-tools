# GLORI-tools v1.0

* Background
* Installation and Requirement
* Example and Usage
* Maintainers and Contributing
* License

# Background
## GLORI
We developed an absolute m6A quantification method (“GLORI”) that is conceptually similar bisulfite sequencing-based quantification of DNA 5-methylcytosine.
GLORI relies on glyoxal and nitrite-mediated deamination of unmethylated adenosines while keeping m6A intact, thereby achieving specific and efficient m6A detection.

## GLORI-tools

GLORI-tools is a bioinformatics pipeline tailored for the analysis of high-throughput sequencing data generated by GLORI。
GLORI-tools currently works, but is still being optimized for a better user experience.

# Installation and Requirement
## Installation
GLORI-tools is written in Python3 and is executed from the command line. To install GLORI-tools simply download the code and extract the files into a GLORI-tools installation folder.

## GLORI-tools needs the following tools to be installed and ideally available in the PATH environment:
* STAR ≥ v2.7.5c
* bowtie (bowtie1) version ≥ v1.3.0
* samtools ≥ v1.10
* python ≥ v3.8.3

## GLORI-tools needs the following python package to be installed:
pysam,pandas,argparse,time,collections,os,sys,re,subprocess,multiprocessing,numpy,scipy,math,sqlite3,Bio,statsmodels

# Example and Usage:

## 2.1 Generated annotation files 
1) download files for annotation (required): 
* '<wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt>'
