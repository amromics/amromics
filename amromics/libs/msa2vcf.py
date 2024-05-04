#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys
import traceback
from operator import itemgetter
from itertools import groupby

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def iupac_to_base(base):
    lookup_iupac = { 'R' : ['A', 'G'],
                     'Y' : ['C', 'T'],
                     'S' : ['G', 'C'],
                     'W' : ['A', 'T'],
                     'K' : ['G', 'T'],
                     'M' : ['A', 'C'],
                     'B' : ['C', 'G', 'T'],
                     'D' : ['A', 'G', 'T'],
                     'H' : ['A', 'C', 'T'],
                     'V' : ['A', 'C', 'G'] }

    if lookup_iupac.get(base):
        return lookup_iupac.get(base)
    else:
        return base


def get_ref_seq(msa, refname):
    with open(msa) as f:
        for name, seq, qual in readfq(f):
            if name == refname:
                return seq.upper()
                break


def group_dels(dels):
    grouped_dels = []
    for k,g in groupby(enumerate(dels),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        grouped_dels.append((group[0] - 1 ,len(group) + 1))
    return grouped_dels


def group_ins(ins):
    res =  [(el - 1, ins.count(el) + 1) for el in ins]
    return set(res)


def fix_complex_vars(all_vars, rseq, qseq):
    positions = [i[2] for i in all_vars]
    duplicate_positions = set([x for x in positions if positions.count(x) > 1])

    for pos in duplicate_positions:
        duplicate_variants = [x for x in all_vars if x[2] == pos]
        for variant in duplicate_variants:
            if variant[0] == "snp":
                all_vars.remove(variant)

            elif variant[0] == "ins":
                variant[5] += 1
                variant[2] -= 1
                variant[4] -= 1
                variant[1] = rseq[variant[4]]
                variant[3] = qseq[variant[4]:variant[4]+variant[5]]

            elif variant[0] == "del":
                variant[5] += 1
                variant[2] -= 1
                variant[4] -= 1
                variant[3] = rseq[variant[4]]
                variant[1] = rseq[variant[4]:variant[4]+variant[5]]
    return all_vars


def update_snps(qseq, rseq,isProt=False):
    dels = []
    ins = []
    snps = []

    for qpos, qbase in enumerate(qseq):
        ins_aware_pos = qpos - len(ins)

        if qbase != rseq[qpos]:
            if qbase == "-":
                # deletion
                dels.append(ins_aware_pos)

            elif rseq[qpos] == "-":
                # insertion
                ins.append(ins_aware_pos)

            else:
                # snp
                snps.append((ins_aware_pos, qbase))

        elif qbase == "-" and rseq[qpos] == "-":
            # insertion but not in this sequence
            ins.append(ins_aware_pos)

    all_vars = []

    if dels:
        grouped_dels = group_dels(dels)
        for start,length in group_dels(dels):
            ins_correction = len([i for i in ins if i < start])
            var = ["del", rseq[start+ins_correction:start+ins_correction+length], start, rseq[start+ins_correction], start+ins_correction, length ]

            all_vars.append(var)

    if snps:
        for position,base in snps:
            ins_correction = len([i for i in ins if i < position])
            var = ["snp", rseq[position+ins_correction], position, base, position+ins_correction, 1]
            all_vars.append(var)

    if ins:
        for start,length in group_ins(ins):
            ins_correction = len([i for i in ins if i < start])
            if rseq[start+ins_correction:start+ins_correction+length] ==  qseq[start:start+length]:
                # insertion but not in this sequence
                continue
            var = ["ins", rseq[start+ins_correction], start, qseq[start:start+length], start+ins_correction, length]
            all_vars.append(var)

    if len(all_vars) >= 1:
        fixed_vars = fix_complex_vars(all_vars, rseq, qseq)
        if not isProt:
            for var in fixed_vars:
                deconvolute_IUPAC(var)
                #print(var)
        fixed_vars_s = sorted(fixed_vars, key=lambda x: x[2])
        return fixed_vars_s
    else:
        return None


def deconvolute_IUPAC(var):
    num_alts = 1

    for base in var[3]:
        no_iupac = iupac_to_base(base)
        if isinstance(no_iupac, list):
            num_alts = len(no_iupac)
            if var[1] not in no_iupac:
                print(var[1])
            no_iupac.remove(var[1])
            var[3] = ','.join(no_iupac)

    var.append(round(1 / num_alts, 2))
    return var


def make_vcf(snps, qname, rname, keep_n):

    vcflines = []

    # header
    vcflines.append("##fileformat=VCFv4.2")
    vcflines.append("##source="+os.path.basename(sys.argv[0]))
    vcflines.append("##contig=<ID=" + rname + ">")
    vcflines.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">")
    vcflines.append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">")
    vcflines.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    vcflines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+qname)

    # variants
    for line in snps:
         vcf_line = None

         alt_no_n = re.sub(r'(N|-)+', '', line[3])

         alt_no_gap = re.sub(r'-+', '', line[3])

         if alt_no_n:
             if not line[1] == alt_no_n:
                 vcf_line = '\t'.join([rname, str(line[2] + 1 ), ".", line[1], alt_no_gap, ".", "PASS", "DP=1", "GT", "1"])

         else:
             if keep_n:
                 vcf_line = '\t'.join([rname, str(line[2] + 1 ), ".", line[1], alt_no_gap, ".", "PASS", "DP=1", "GT", "1"])

         if vcf_line:
             vcflines.append(vcf_line)
    return vcflines


def write_vcf(vcflines, qname,output_dir):
    filename=output_dir+"/"+qname +".vcf"
    with open(filename, 'w') as f:
        for line in vcflines:
            f.write(line+"\n")
    return filename


def remove_terminal_gapns(seq):
    return re.sub(r'(N|-)*$', '', seq)


def go(msa,refname,output_dir,isProt=False,keep_n=False):
    map_gene_vcf={}
    try:
        refseq = get_ref_seq(msa, refname)

        if not refseq:
            print(refname + " not found\n")
            sys.exit(1)

        with open(msa) as f:
            for name, seq, qual in readfq(f):
                if name == refname:
                    continue
                else:
                    #print(name)
                    snps = update_snps(remove_terminal_gapns(seq.upper()), refseq,isProt)

                    if snps:
                        vcflines = make_vcf(snps, name, refname, keep_n)
                        if isProt:
                            gene_vcf_file=write_vcf(vcflines, name+".prot",output_dir)
                        else:
                            gene_vcf_file=write_vcf(vcflines, name,output_dir)
                        map_gene_vcf[name]=gene_vcf_file

    except Exception as ex:
        traceback.print_exc()
        print('Error call vcf '+msa+":"+str(ex))
    return map_gene_vcf
