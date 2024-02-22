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
import traceback
from Bio import SeqIO
import pandas as pd
import json
import gzip
from datetime import datetime
import amromics.libs.bioseq as bioseq
import amromics.libs.pangenome as pangenome
import amromics.libs.alignment as alignment
import amromics.libs.phylogeny as phylogeny
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def run_single_sample(sample,extraStep=False, sample_dir='.', assembly_method='spades', threads=0, memory=50, trim=False, overwrite=None, timing_log=None):
    #handle assembly input, ignore spades and bwa:
    sample['execution_start'] =  datetime.now().strftime("%Y-%m-%d %H:%M:%S")
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
        #TODO: get the number of bases from fastq json file. If not trim, iterate through the read file to count the number of bases
        if trim:
            if 'pe1' in reads and 'pe2' in reads:
                reads['pe1'],reads['pe2'] = preprocess.trim_fastp(sample['id'],reads, overwrite=overwrite, base_dir=sample_dir, timing_log=timing_log,threads=threads)
            elif 'se' in reads:
                reads['se'] = preprocess.trim_fastp(sample['id'],reads,
                        base_dir=sample_dir, overwrite=overwrite, timing_log=timing_log,threads=threads)
        #estimate gsize if not provided
        if sample['gsize']==None:            
            sample['gsize'] = preprocess.estimate_gsize_mash(sample['id'], reads,
                    base_dir=sample_dir, overwrite=overwrite, threads=threads, memory=memory,timing_log=timing_log)
            logger.info('Genome size of {sample_id} is estimated to be {gsize}'.format(sample_id=sample['id'],gsize=sample['gsize']))
        #subsampling to 100X if needed
        #TODO: by now, we should know the genome size and the coverage. If genome size is not known use the upper bound= 8M            
        reads=preprocess.subsample_seqtk(sample['id'], reads, base_dir=sample_dir, overwrite=overwrite, threads=threads, memory=memory, gsize=sample['gsize'], timing_log=timing_log)
        #run assembly by spades
        if assembly_method=='spades':
            sample['assembly'] = assembler.assemble_spades(sample['id'], reads, base_dir=sample_dir, threads=threads, memory=memory,overwrite=overwrite, timing_log=timing_log)
        elif assembly_method=='skesa':
            sample['assembly'] = assembler.assemble_skesa(sample['id'], reads, base_dir=sample_dir, threads=threads, memory=memory,overwrite=overwrite, timing_log=timing_log)
        else:
            raise Exception(f'Unknown assembly method {assembly_method}')
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
    sample['execution_end'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return sample
def run_collection(report,gff_dir,ffn_dir,faa_dir, base_dir='.',threads=8,progressive=False, overwrite=None,memory=50, timing_log=None,method='roary',genetree=True, tree='fasttree'):
    try:
        starttime = datetime.now()
        if method=='roary':
            stime = datetime.now()
            report['pan'] = pangenome.run_roary(gff_dir, threads=threads, base_dir=base_dir,overwrite=overwrite,timing_log=timing_log)
            elapsed = datetime.now() - stime
            logger.info(f'Roary -- time taken {str(elapsed)}')
            stime = datetime.now()
            report['alignments'] = alignment.runGeneAlignment(report['pan'],14, ffn_dir,faa_dir,overwrite=overwrite,collection_dir=base_dir, threads=threads,timing_log=timing_log)
            elapsed = datetime.now() - stime
            logger.info(f'Alignment from roary -- time taken {str(elapsed)}')
        if method=='panta':
            stime = datetime.now()
            report['pan'] = pangenome.run_panta_cmd(gff_dir, threads=threads, base_dir=base_dir,progressive=progressive,overwrite=overwrite,timing_log=timing_log)
            elapsed = datetime.now() - stime
            logger.info(f'Panta -- time taken {str(elapsed)}')
            stime = datetime.now()
            #report['alignments'] = alignment.runGeneAlignment(report['pan'],8, ffn_dir,faa_dir,overwrite=overwrite,collection_dir=base_dir, threads=threads,timing_log=timing_log)
            report['alignments'] = alignment.copyClustersPanta(base_dir,report['pan'])

            elapsed = datetime.now() - stime
            logger.info(f'Gene alignment -- time taken {str(elapsed)}')
            stime = datetime.now()
            report['alignments']=alignment.runVCFCallingFromGeneAlignment(report['pan'],overwrite=overwrite,collection_dir=base_dir, threads=threads,timing_log=timing_log)
            elapsed = datetime.now() - stime
            logger.info(f'Call VCF -- time taken {str(elapsed)}')
        if genetree:
            stime = datetime.now()
            report['alignments']  = phylogeny.run_gene_phylogeny_iqtree(report['pan'], collection_dir=base_dir,overwrite=overwrite, threads=threads,timing_log=timing_log)
            elapsed = datetime.now() - stime
            logger.info(f'IQTree for genes -- time taken {str(elapsed)}')
        stime = datetime.now()
        report['coregene'] = alignment.create_core_gene_alignment(report['pan'], collection_dir=base_dir,overwrite=overwrite, threads=threads,timing_log=timing_log)
        elapsed = datetime.now() - stime
        logger.info(f'Create core gene -- time taken {str(elapsed)}')
        if tree=='iqtree':
            stime = datetime.now()
            report['phylogeny']  = phylogeny.run_species_phylogeny_iqtree(report['pan'] ,collection_dir=base_dir,overwrite=False, threads=threads,timing_log=timing_log)
            elapsed = datetime.now() - stime
            logger.info(f'IQTree for species -- time taken {str(elapsed)}')
        if 'phylogeny' not in report.keys() or report['phylogeny']==None:
            stime = datetime.now()
            report['phylogeny']  = phylogeny.run_species_phylogeny_fastree(report['pan'] ,collection_dir=base_dir,overwrite=False, threads=threads,timing_log=timing_log)
            elapsed = datetime.now() - stime
            logger.info(f'Fasttree for species -- time taken {str(elapsed)}')
        elapsed = datetime.now() - starttime
        logger.info(f'Run Collection -- time taken {str(elapsed)}')
    except Exception as error:
        traceback.print_exc()
        logger.error(error)
        return None
    return report
