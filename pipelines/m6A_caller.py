
"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used to call m6A sites based on background"""
"""Input: [.txt]"""

import os, sys
import argparse
from collections import defaultdict,Counter
import time
from time import strftime
import numpy as np
import math
import scipy.stats


def read_CR(input):
    overall = []
    n = 0
    line = input.readline()
    while (line):
        if line.startswith("#"):
            line = line.strip().split("\t")
            if len(line) == 5:
                n += 1
                gene, total, Cs, nonCR, Alength = line
                gene = gene.replace("#", "")
                if int(total)/int(Alength) < options.AG_number or float(nonCR) >= options.gene_CR:
                    nonCRs[gene] = 1.0
                    averageDepths[gene] = int(total) / int(Alength)
                else:
                    nonCRs[gene] = float(nonCR)
                    averageDepths[gene] = int(total)/int(Alength)
                    overall.append(float(nonCR))
                if gene in control:
                    if int(total) > options.AG_number:
                        if gene != "ALL":
                            control[gene] = float(nonCR)
                if gene == "ELSE":
                    sys.stderr.write("#No annotation\t%.6f\n" % (float(nonCR)))
                elif gene == "ALL":
                    sys.stderr.write("#Overall\t%s\n" % (float(nonCR)))
            else:
                stat, value = line
            # sys.stderr.write("%s\t%s\n" % (stat,value) )
            # if stat == "#Median":
            # nonCRs["Median"] = float(value)
            # elif stat == "#Mean":
            # nonCRs["Mean"] = float(value)
        else:
            if len(overall) > 0:
                sys.stderr.write("#%d of %d genes have AG coverage >= %d\n" % (len(overall), n, options.AG_number))
                nonCRs["Median"] = np.median(overall)
                nonCRs["Mean"] = np.mean(overall)
                nonCRs["discard"] = 1.0
                sys.stderr.write("#Median\t%.8f\n" % nonCRs["Median"])
                sys.stderr.write("#Mean\t%.8f\n" % nonCRs["Mean"])
                sys.stderr.write("#90%%\t%.8f\n" % np.percentile(overall, 90))
                sys.stderr.write("#75%%\t%.8f\n" % np.percentile(overall, 75))
                sys.stderr.write("#50%%\t%.8f\n" % np.percentile(overall, 50))
                sys.stderr.write("#25%%\t%.8f\n" % np.percentile(overall, 25))
                sys.stderr.write("#10%%\t%.8f\n" % np.percentile(overall, 10))
                break
            else:
                sys.stderr.write("#%d of %d genes have AG coverage >= %d\n" % (len(overall), n, options.AG_number))
                break
        line = input.readline()  # Note that the first line of data was read
    control_nonCR = []
    if control:
        for key, value in control.items():
            if value != "None":
                control_nonCR.append(value)
        nonCRs["control"] = np.median(control_nonCR)
        sys.stderr.write("[%s] Control number: %d of %d; median: %.6f; mean: %.6f\n" \
                         % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()), len(control_nonCR), len(control),
                            np.median(control_nonCR), np.mean(control_nonCR)))
    return line, input


def call_m6A(line, col):
    # print("************1",line)
    chr, pos_1, dir, gene, name, trans, isoform, biotype, total_cov, AG, A_count, T_count, C_count, G_count = line[0:14]
    total_cov = int(total_cov)
    AG = int(AG)
    if col == -1:
        total_cov_col = total_cov
        AG_col = AG
        A_count_col = int(A_count)
    else:
        total_cov_col, AG_col, A_count_col = map(int,line[col].split(";")[1].split(","))
    if AG_col >= options.coverage and A_count_col >= options.count:
        # print("************2",line)
        # time.sleep(2)
        ratio = A_count_col / (AG_col + 0.0)
        if ratio < options.ratio or AG / (total_cov + 0.0) < options.var_ratio or AG_col / (total_cov_col + 0.0) < options.var_ratio:
            return None
        else:
            nonCR_gene = nonCRs.get(gene)
            averageDepth_gene = averageDepths.get(gene)
            if nonCR_gene is None:
                nonCR_gene = nonCRs[options.non_anno]
                averageDepth_gene = averageDepths[options.non_anno]
            if options.conversion_rate == "gene":
                nonCR = nonCR_gene
                averageDepth = averageDepth_gene
            elif options.conversion_rate == "overall":
                if nonCR_gene == 1.0:
                    nonCR = 1.0
                    averageDepth = averageDepth_gene
                else:
                    nonCR = nonCRs.get("ALL")
                    averageDepth = averageDepths.get("ALL")
            elif options.conversion_rate == "control":
                if nonCR_gene == 1.0:
                    nonCR = 1.0
                else:
                    nonCR = nonCRs.get("control")
            if options.method == "binomial":
		try:
                    pvalue = scipy.stats.binom_test(A_count_col, n=AG_col, alternative='greater', p=nonCR)
                except AttributeError:
                    pvalue = scipy.stats.binomtest(A_count_col, n=AG_col, alternative='greater', p=nonCR)
            elif options.method == "fisher":
                pass
            elif options.method == "poisson":
                pvalue = scipy.stats.poisson.sf(A_count_col, int(math.ceil(nonCR * AG_col)))
            # print("************3")
            if pvalue < options.pvalue:
                # print("************4")
                Signal = total_cov_col / (total_cov + 0.0)
                if Signal >= options.signal:
                    # print("************5")
                    if nonCR_gene < 0.00001:
                        nonCR_format = format(nonCR_gene, '.1e')
                    else:
                        nonCR_format = str(round(nonCR_gene, 5))
                    if pvalue < 0.001:
                        P_format = format(pvalue, '.1e')
                    else:
                        P_format = str(round(pvalue, 3))
                    # print("**********2")
                    LINE = "\t".join([chr, pos_1, dir, gene, name, trans, isoform, biotype, nonCR_format])
                    LINE = LINE + "\t" + ",".join(line[8:14])
                    LINE = LINE + "\t" + "\t".join(
                        [str(total_cov_col), str(AG_col), str(A_count_col), str(round(ratio, 5)), P_format,
                         str(round(Signal, 6)),str(averageDepth)])
                    LINE = LINE + "\t" + "|".join([i for i in line[14:]]) + "\n"
                    # print("**********6")

                    return LINE
    else:
        return None


if __name__ == "__main__":
    description = """
	"""
    parser = argparse.ArgumentParser(prog="m6A_caller", fromfile_prefix_chars='@', description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)
    # Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i", "--input", dest="input", required=True, help="input pileup, sorted")
    group_required.add_argument("-o", "--output", dest="output", required=True,
                                help="output prefix. file names are [prefix].[A-cutoff].txt")
    # Filter
    group_site = parser.add_argument_group("m6A filter")
    group_site.add_argument("-c", "--coverage", dest="coverage", default=15, type=int, help="A+G coverage, default=10")
    group_site.add_argument("-C", "--count", dest="count", default=5, type=int,
                            help="A count, below which the site will not count, default=3")
    group_site.add_argument("-r", "--ratio", dest="ratio", default=0.1, type=float, help="m6A level/ratio, default=0.1")
    group_site.add_argument("-p", "--pvalue", dest="pvalue", default=0.005, type=float, help="pvalue, default=0.05")
    group_site.add_argument("-s", "--signal", dest="signal", default=0.8, type=float,
                            help="signal ratio, equals coverage(under A-cutoff)/coverage, default=0.8")
    group_site.add_argument("-R", "--var_ratio", dest="var_ratio", default=0.8, type=float,
                            help="the ratio cutoff of AG/Total to filter sequencing/mapping errors, default=0.8")
    group_site.add_argument("-g", "--gene", dest="gene_CR", default=0.2, type=float,
                            help="conversion rate, over which a gene will be discarded, default=0.1")
    group_site.add_argument("-N", "--AG", dest="AG_number", default=0, type=int,
                            help="AG count, below which a gene will be discarded, default=0")
    # Statistics
    group_stat = parser.add_argument_group("Statistic method")
    group_stat.add_argument("--method", dest="method", default="binomial", choices=['binomial', 'fisher', 'poisson'],
                            help="statistical method: binomial, fisher exact test, or poisson, default=binomial")
    group_stat.add_argument("--CR", dest="conversion_rate", default="gene", choices=['gene', 'overall', 'control'],
                            help="conversion rate used: gene or overall, default=gene")
    group_stat.add_argument("--NA", dest="non_anno", default="ELSE",
                            choices=['ELSE', 'Median', 'Mean', 'ALL', 'discard'],
                            help="which CR to use if no aene annotation, default=ELSE")
    group_stat.add_argument("--control", dest="control", help="control list, median non-conversion rate will be used")
    group_stat.add_argument("--cutoff", dest="A_cutoffs", default="3",
                            help="A-cutoffs, 1-10,15,20 or None, seperated by comma, default=3,None")
    # Version
    group_other = parser.add_argument_group("Other")
    group_other.add_argument("--version", action="version", version="%(prog)s 1.0")
    options = parser.parse_args()

    sys.stderr.write("[%s] Running, pid [%s]\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()), str(os.getpid())))
    sys.stderr.write("[%s] Reading header...\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    control = {}
    if options.conversion_rate == "control":
        if not options.control:
            raise IOError("control list not given, exit.\n")
        else:
            with open(options.control, 'r') as input:
                for line in input.readlines():
                    control[line.strip()] = "None"

    A_cutoffs = options.A_cutoffs.split(",")
    for i in range(len(A_cutoffs)):
        if A_cutoffs[i] not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '15', '20', 'None']:
            raise ValueError("%s not in A-cutoffs, exit.\n" % A_cutoffs[i])

    sys.stderr.write(
        "[%s] A-cutoffs = [%s]\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()), ", ".join(A_cutoffs)))

    A_cutoff_cols = {}
    col_start = 14
    n = 0
    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '15', '20', 'None'][:-1]:
        A_cutoff_cols[i] = 14 + n
        n += 1
    A_cutoff_cols["None"] = -1

    sys.stderr.write("[%s] Analyzing, using [%s] conversion rates...\n" % (
    strftime("%Y-%m-%d %H:%M:%S", time.localtime()), options.conversion_rate))
    size = os.path.getsize(options.input)
    if size > 0:
        with open(options.input, 'r') as input:
            nonCRs = {}
            averageDepths = {}
            results = defaultdict(list)
            line, input = read_CR(input)
            while (line):
                line = line.strip().split("\t")
                for A_cutoff in A_cutoffs:
                    col = A_cutoff_cols.get(A_cutoff)
                    # print(line)
                    result = call_m6A(line, col)
                    # print("*********1",result)
                    if result is not None:
                        results[A_cutoff].append(result)
                line = input.readline()
        sys.stderr.write("[%s] Writing to disk...\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        for key, records in results.items():
            with open(options.output + "." + key + ".txt", 'w') as output:
                for record in records:
                    output.write(record)
            sys.stderr.write("[%s] %s written.\n" % (
            strftime("%Y-%m-%d %H:%M:%S", time.localtime()), options.output + "." + key + ".txt"))
        sys.stderr.write("[%s] Finished successfully!\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
