# -*- coding: utf-8 -*-
"""
    Wrapper for tools

----------------
Revision history:
----------------
2019-08-17: Amromics created

"""
from __future__ import division, print_function, absolute_import

import os
import shutil
import csv
import logging
import multiprocessing
import amromics.libs.bioseq as bioseq
import amromics.libs.element_finder as element_finder
import amromics.libs.mlst as mlst
import amromics.libs.amr as amr
import amromics.libs.annotation as annotation
import amromics.libs.assembler as assembler
import amromics.libs.preprocess as preprocess
import amromics.libs.qc as qc
import amromics.libs.taxonomy as taxonomy
from Bio import SeqIO
import pandas as pd
import json
import datetime
import gzip

import amromics.libs.bioseq as bioseq
import amromics.libs.pangenome as pangenome
import amromics.libs.alignment as alignment
import amromics.libs.phylogeny as phylogeny
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def run_single_sample(sample,extraStep=False, sample_dir='.', threads=0, memory=50, trim=False, overwrite=None, timing_log=None):
    #handle assembly input, ignore spades and bwa:
    sample['execution_start'] =  datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    reads=None
    if sample['input_type'] in ['asm', 'assembly']:
        sample['assembly'] = assembler.get_assembly(sample['id'], sample['files'],base_dir=sample_dir)
    elif sample['input_type'] in ['pacbio-raw', 'pacbio-hifi', 'pacbio-corr','nano-raw', 'nano-hq', 'nano-corr']:
        read_file = sample['files'].split(';')
        if len(read_file) > 1:
            raise Exception('All reads should be put in one file only')
        reads={}
        reads['long-read'] = read_file[0]
        sample['assembly'] = assembler.assemble_flye(sample['id'],reads, input_type=sample['input_type'], base_dir=sample_dir, threads=threads,timing_log=timing_log,gsize=sample['gsize'])
    elif sample['input_type'] in ['Illumina paired-end', 'Illumina single-end', 'Illumina']:
        pe_files = sample['files'].split(';')
        reads={}
        if len(pe_files) == 1:
            reads['se']=pe_files[0]
        elif len(pe_files) == 2:
            reads['pe1']=pe_files[0]
            reads['pe2']=pe_files[1]
        else:
            raise Exception('There should be one or two input files')
        #preprocesing
        preprocess.rename_reads(reads)
        if trim:
            if 'pe1' in reads and 'pe2' in reads:
                reads['pe1'],reads['pe2'] = preprocess.trim_fastp(sample['id'],reads, overwrite=overwrite, base_dir=sample_dir, timing_log=timing_log,threads=threads)
            elif 'se' in reads:
                reads['se'] = preprocess.trim_fastp(sample['id'],reads,
                        base_dir=sample_dir, overwrite=overwrite, timing_log=timing_log,threads=threads)
        #estimate gsize if not provided
        if sample['gsize']==None:
            sample['gsize']=preprocess.estimate_gsize_mash(sample['id'], reads,
                    base_dir=sample_dir, overwrite=overwrite, threads=threads, memory=memory,timing_log=timing_log)
        #subsampling to 100X if needed
        reads=preprocess.subsample_seqtk(sample['id'], reads, base_dir=sample_dir, overwrite=overwrite, threads=threads, memory=memory, gsize=sample['gsize'], timing_log=timing_log)
        #run assembly by spades
        sample['assembly'] = assembler.assemble_spades(sample['id'], reads, base_dir=sample_dir, threads=threads, memory=memory,timing_log=timing_log)
        #sample['assembly'] = assembler.assemble_shovill(sample['id'],reads, base_dir=sample_dir, threads=0, memory=memory,trim=trim,timing_log=timing_log,gsize=sample['gsize'])
    elif sample['input_type'] in ['gff']:
        sample['annotation_gff'],sample['annotation_faa'],sample['annotation_ffn'],sample['annotation_fna']=annotation.parseGFF(sample['id'],sample['files'],base_dir=sample_dir, overwrite=overwrite)
        sample['assembly'] = assembler.get_assembly_from_gff(sample['id'], sample['files'],base_dir=sample_dir)

    # if extraStep and not reads==None :
    #     sample['se_bam']=qc.map_reads_to_assembly_bwamem(sample['id'],sample['assembly'],reads, base_dir=sample_dir, threads=0, memory=memory,timing_log=timing_log)
    if extraStep and not reads==None:
        sample['qc'] =qc.qc_reads(sample['id'],reads, base_dir=sample_dir, threads=0, timing_log=timing_log)
    #QUAST to check quality
    # if extraStep:
    #     sample['quast']=qc.assembly_eval(sample['id'],sample['assembly'], base_dir=sample_dir, threads=0, timing_log=timing_log)
    if extraStep:
        sample['taxonomy']=taxonomy.species_identification_kraken(sample['id'],sample['assembly'], base_dir=sample_dir, timing_log=timing_log,threads=threads)
    #QUAST to check quantity
    #FastQC, + MultiQC
    if not 'gram' in sample.keys():
        sample['gram']=None
    if not 'annotation_gff' in sample.keys() or sample['annotation_gff']==None:
        sample['annotation_gff'],sample['annotation_faa'],sample['annotation_ffn'],sample['annotation_fna'],sample['annotation_gbk'] = annotation.annotate_prokka(sample['id'],sample['assembly'],sample['genus'],sample['species'],sample['strain'],sample['gram'], base_dir=sample_dir,timing_log=timing_log, threads=threads)
    sample['mlst']  = taxonomy.detect_mlst(sample['id'],sample['assembly'], base_dir=sample_dir, threads=threads)
    if extraStep:
        sample['resistome'],sample['point'],sample['virulome'] = amr.detect_amr_amrfinder(sample['id'],sample['annotation_faa'],sample['annotation_fna'],sample['annotation_gff'],sample['genus'],sample['species'], base_dir=sample_dir,timing_log=timing_log, threads=threads)
    else:
        sample['resistome'] = amr.detect_amr_abricate(
        sample['id'],sample['assembly'], base_dir=sample_dir, threads=threads, timing_log=timing_log)
    sample['virulome'] = amr.detect_virulome(sample['id'],sample['assembly'], base_dir=sample_dir, threads=threads)
    sample['plasmid']  = amr.detect_plasmid(sample['id'],sample['assembly'], base_dir=sample_dir, threads=threads)
    if extraStep:
        sample['pmlst']=amr.detect_pmlst(sample['id'],sample['assembly'], base_dir=sample_dir, threads=threads)
    # if extraStep:
    #     sample['is']=amr.detect_insertion_sequence(sample['id'],sample['assembly'], base_dir=sample_dir, threads=threads)
    # if extraStep:
    #     sample['integron']=amr.detect_integron(sample['id'],sample['assembly'], base_dir=sample_dir,timing_log=timing_log, threads=threads)
    #sample=detect_prophage(sample, base_dir=base_dir, threads=threads)
    sample['execution_end'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return sample
def run_collection(report,gff_dir,ffn_dir, base_dir='.',threads=8,progressive=False, overwrite=None,memory=50, timing_log=None,method='roary',genetree=True):
    try:
        if method=='roary':
            report['roary'] = pangenome.run_roary(gff_dir, threads=threads, base_dir=base_dir,overwrite=overwrite,timing_log=timing_log)
            report['alignments'] = alignment.runGeneAlignment(report['roary'],14, ffn_dir,overwrite=overwrite,collection_dir=base_dir, threads=threads,timing_log=timing_log)

        if method=='panta':
            report['roary'] = pangenome.run_panta_cmd(gff_dir, threads=threads, base_dir=base_dir,progressive=progressive,overwrite=overwrite,timing_log=timing_log)
            report['alignments'] = alignment.runGeneAlignment(report['roary'],8, ffn_dir,overwrite=overwrite,collection_dir=base_dir, threads=threads,timing_log=timing_log)
            report['alignments']=runVCFCallingFromGeneAlignment(report['roary'],overwrite=overwrite,collection_dir=base_dir, threads=threads,timing_log=timing_log)
        if genetree:
            report['alignments']  = phylogeny.run_gene_phylogeny_iqtree(report['roary'], collection_dir=base_dir,overwrite=overwrite, threads=threads,timing_log=timing_log)
        report['coregene'] = alignment.create_core_gene_alignment(report['roary'], collection_dir=base_dir,overwrite=overwrite, threads=threads,timing_log=timing_log)

        #report['phylogeny']  = phylogeny.run_species_phylogeny_iqtree(report['roary'] ,collection_dir=base_dir,overwrite=False, threads=threads,timing_log=timing_log)
        if 'phylogeny' not in report.keys() or report['phylogeny']==None:

            report['phylogeny']  = phylogeny.run_species_phylogeny_fastree(report['roary'] ,collection_dir=base_dir,overwrite=False, threads=threads,timing_log=timing_log)
    except Exception as error:
        logger.error(error)
    return report
