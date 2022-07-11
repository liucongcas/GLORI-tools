#!/usr/bin/env python

"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used for formating mpileup files"""
"""Input: [mpileup files]"""


import os, sys
import argparse
from collections import defaultdict, Counter
import copy
import sqlite3
import time
from time import strftime
import numpy as np



def create_table(cursor):  # id text PRIMARY KEY
    sql = """ CREATE TABLE IF NOT EXISTS 
			locations (
				chr text NOT NULL,
				pos integer NOT NULL,
				dir text NOT NULL,
				gene text NOT NULL,
				trans text NOT NULL,
				name text NOT NULL,
				isoform text NOT NULL,
				biotype text NOT NULL
			); 
		   """
    cursor.execute(sql)


def insert_var(cursor, chr, pos, dir, gene, trans, name, isoform, biotype):
    cursor.execute(
        "insert into locations (chr, pos, dir, gene, trans, name, isoform, biotype) values ( ?, ?, ?, ?, ?, ?, ?, ?)",
        (chr, pos, dir, gene, trans, name, isoform, biotype))  # id


def insert_many_var(cursor, rows):
    cursor.executemany("""
	insert into locations (chr, pos, dir, gene, trans, name, isoform, biotype) values ( ?, ?, ?, ?, ?, ?, ?, ?)
	""", rows)


def create_index(cursor):
    cursor.execute("CREATE INDEX locations_index ON locations (chr,dir,pos)")


def build_database():
    global database, cursor, conn
    with open(options.database, 'r') as input:
        line = input.readline()
        if options.method == "sql":
            conn = sqlite3.connect(":memory:")
            cursor = conn.cursor()
            create_table(cursor)
            cursor.execute("PRAGMA cache_size=65536")
            cursor.execute("PRAGMA page_size=65536")
            cursor.execute('PRAGMA temp_store=MEMORY')
            cursor.execute("PRAGMA synchronous=OFF")
            cursor.execute('PRAGMA journal_mode=MEMORY')
        elif options.method == "dict":
            pass
        n = 0
        rows = []
        while (line):
            line = line.strip().split("\t")
            chr = line[0]
            pos_1 = line[2]
            dir = line[3]
            gene = line[4]
            trans = line[5]
            name = line[6]
            isoform = line[7]
            biotype = line[8]
            if options.method == "sql":
                rows.append([chr, int(pos_1), dir, gene, trans, name, isoform, biotype])
            elif options.method == "dict":
                database[(chr, int(pos_1), dir)] = [gene, name, trans, isoform,
                                                    biotype]  # SELECT gene,name,trans,isoform,biotype
            n += 1
            if n % 1000000 == 0:
                if options.method == "sql":
                    insert_many_var(cursor, rows)
                    rows = []
                sys.stderr.write("[%s] %d items processed...\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()), n))
            line = input.readline()
        if options.method == "sql":
            if rows:
                insert_many_var(cursor, rows)
                rows = []
    if options.method == "sql":
        sys.stderr.write("[%s] All loaded. Creating index...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        create_index(cursor)
        sys.stderr.write("[%s] All finished!\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    elif options.method == "dict":
        pass


def write_CR():
    overall = []
    lines = []
    with open(options.conversion, 'w') as CR_output:
        if all_Aov != 0:
            ratio_nonCR = str(all_A / (all_Aov + 0.0))
            CR_output.write("\t".join(["#ALL", str(all_Aov), str(all_A), ratio_nonCR,str(all_mappinglength)]))
            CR_output.write("\n")
            for key in sorted(result_cov.keys()):
                Acounts = result_A.get(key)
                cov = result_cov.get(key)
                gene_length = result_genelength.get(key)
                if cov is not None and cov > 0:
                    if Acounts is None:
                        Acounts = 0
                    nonCR = Acounts / (cov + 0.0)
                    line = "\t".join([key, str(cov), str(Acounts), str(nonCR),str(gene_length)]) + "\n"
                    lines.append(line)
                    if cov > 0:
                        overall.append(nonCR)
            if len(overall)>0:
                CR_output.write("#Median\t%.8f\n" % np.median(overall))
                CR_output.write("#Mean\t%.8f\n" % np.mean(overall))
                CR_output.write("#90%%\t%.8f\n" % np.percentile(overall, 90))
                CR_output.write("#75%%\t%.8f\n" % np.percentile(overall, 75))
                CR_output.write("#50%%\t%.8f\n" % np.percentile(overall, 50))
                CR_output.write("#25%%\t%.8f\n" % np.percentile(overall, 25))
                CR_output.write("#10%%\t%.8f\n" % np.percentile(overall, 10))
            else:
                pass
            for line in lines:
                CR_output.write("#")
                CR_output.write(line)
        else:
            pass


def search(chr, dir, pos, database=None, cursor=None):
    if options.method == "sql":
        cursor.execute("SELECT gene,name,trans,isoform,biotype FROM locations WHERE chr=? AND dir=? AND pos=?",
                       (chr,dir,pos))
        rows = cursor.fetchall()
        if rows:
            return [str(i) for i in rows[0]]
        else:
            return None
    else:
        return database.get((chr,dir,pos))


def Tideup(chr, pos_1, dir, PIPLEUPS, SURROUNDINGS, database=None, cursor=None):
    global output, all_A, all_Aov,all_mappinglength
    if database != 'Notusingdatabase':
        if options.method == "sql":
            query = search(chr, dir, pos_1, cursor=cursor)
        elif options.method == "dict":
            query = search(chr, dir, pos_1, database=database)
        if query is not None:
            gene, name, trans, isoform, biotype = query
            GENE = gene
        else:
            gene, name, trans, isoform, biotype = ["NA"] * 5
            GENE = "ELSE"
    else:
        gene, name, trans, isoform, biotype = ["NA"] * 5
        GENE = "ELSE"
    # if dir=="-":
    #     PIPLEUPS2 = ['G' if i == "C" else i for i in PIPLEUPS]
    #     PIPLEUPS3 = ['A' if i == "T" else i for i in PIPLEUPS2]
    #     PIPLEUPS4 = ['T' if i == "A" else i for i in PIPLEUPS3]
    #     PIPLEUPS5 = ['C' if i == "G" else i for i in PIPLEUPS4]
    #     PIPLEUPS = PIPLEUPS5
    bases = sorted(Counter(zip(PIPLEUPS, SURROUNDINGS)).items(), key=lambda x: x[0][1])
    counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'total': 0, 'AG': 0}
    bins = {i: {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'total': 0, 'AG': 0} for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20]}
    limits = [99999, 20, 15, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    """bases: [(('A', '0'), 2), (('A', '1'), 2), (('A', '2'), 1), (('G', '2'), 1), (('G', '3'), 1), (('G', '4'), 1)]"""
    for (base, surr), count in bases:
        if surr > limits[-1]:
            while surr > limits[-1]:
                index = limits.pop()
                bins[index] = copy.copy(counts)
        if base != "N":
            counts[base] += count
            counts["total"] += count
        if base == "A" or base == "G":
            counts["AG"] += count
    for item in limits:
        bins[item] = copy.copy(counts)
    # write to disk
    output.write("\t".join([chr, str(pos_1), dir, gene, name, trans, isoform, biotype]))  # information
    output.write("\t")
    output.write("\t".join([str(counts['total']), str(counts['AG']), str(counts['A']), str(counts['T']), str(counts['C']),
         str(counts['G'])]))  # Overall count
    for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 99999][:-1]:
        output.write("\t")
        output.write(str(i))
        output.write(";")
        output.write(",".join(
            [str(bins[i]['total']), str(bins[i]['AG']), str(bins[i]['A'])]))  # SP=bins[i]['total']/counts['total']
    output.write("\n")

    # Add values to result_cov and result_A
    if counts['AG'] >= 15:
        result_A[GENE] += counts['A']
        result_cov[GENE] += counts['AG']
        result_genelength[GENE] += 1
        all_A += counts['A']
        all_Aov += counts['AG']
        all_mappinglength += 1


def check(chr, pos_1, dir, pileup, surrounding, database, cursor):
    PIPLEUPS = pileup
    SURROUNDINGS = surrounding
    Tideup(chr, pos_1, dir, PIPLEUPS, SURROUNDINGS, database=database, cursor=cursor)


if __name__ == "__main__":

    description = """
	Input should be no redundance
	"""
    parser = argparse.ArgumentParser(prog="", fromfile_prefix_chars='@', description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)
    # Require
    group_require = parser.add_argument_group("Required")
    group_require.add_argument("-i", "--input", dest="input", required=True, help="input pileup, sorted")
    group_require.add_argument("-o", "--output", dest="output", required=True, help="output")  # chr start_0 start
    group_require.add_argument("--CR", dest="conversion", required=True, help="gene conversion rate")
    group_require.add_argument("--method", dest="method", default="sql", choices=["sql", "dict"],
                               help="database type, sql or dictionary")
    # anno
    group_anno = parser.add_argument_group("m6A anno")
    group_anno.add_argument("--db", dest="database",default='Notusingdatabase', help="database, base-gene annotation")
    # group_optional = parser.add_argument_group("Optional")
    group_other = parser.add_argument_group("Other")
    options = parser.parse_args()

    sys.stderr.write("[%s] Reading input...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    if options.database != 'Notusingdatabase':
        if options.method == "sql":
            database = None
            cursor = None
            conn = None
        elif options.method == "dict":
            database = {}
            cursor = None
            conn = None
        build_database()
        sys.stderr.write("[%s] In-memory database connection setup\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    else:
        database=options.database
        cursor = None
        pass
    # init
    result_cov = defaultdict(int)
    result_A = defaultdict(int)
    result_genelength = defaultdict(int)
    all_Aov = 0
    all_A = 0
    all_mappinglength = 0
    with open(options.input, 'r') as input, open(options.output + "_tmp", 'w') as output:
        line = input.readline()
        while (line):
            row = line.strip().split("\t")
            chr = row[0].split("_AG_converted")[0]
            pos_1 = row[1]
            dir = row[2]
            base = row[3]  # as transcript
            pileup = row[4].split(",")
            surrounding = [int(i) for i in row[5].split(",")]
            if dir == "+" and base == "A":
                check(chr, pos_1, dir, pileup, surrounding, database, cursor)
            elif dir == "-" and base == "T":
                check(chr, pos_1, dir, pileup, surrounding, database, cursor)
            line = input.readline()
            # if options.method == "sql":
            # check(chr,pos_1,dir,pileup,surrounding,database=database,cursor=cursor)
            # elif options.method == "dict":
            # check(chr,pos_1,dir,pileup,surrounding,database=database,cursor=cursor)


    sys.stderr.write("[%s] Finished reading input.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    # get conversion rates
    write_CR()
    sys.stderr.write("[%s] Calculating conversion rates...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    os.system("cat %s %s > %s" % (options.conversion, options.output + "_tmp", options.output))
    os.remove(options.output + "_tmp")
    # os.system("mv %s %s " % (options.output + "_tmp", options.output))

    if options.database == "sql":
        database.close()
    else:
        del database
    # sys.stderr.write("[%s] All finished!\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
