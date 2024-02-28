#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    The entry point
"""
from __future__ import division, print_function, absolute_import

import argparse
import logging
import multiprocessing
import os
import sys
import pandas as pd

from amromics.pipeline.analysis import single_genome_analysis, pan_genome_analysis
from amromics.utils.utils import valid_id, software_version
from amromics.db.utils import setup_db,setup_minidb
from amromics import __version__

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s [%(name)s] %(levelname)s : %(message)s')
logger = logging.getLogger(__name__)


def version_func(args):
    software_version([
        'python',
        'blast',
        'trimmomatic',
        'spades',
        'shovill',
        'prokka',
        'mlst',
        'abricate',
#        'roary',
#        'parsnp',
        'iqtree',
        'fastqc',
        'kraken2',
        'multiqc',
#        'quast',
        'isescan',
        'amrfinder'
    ])

def setup_db_func(args):
    if args.initdb:
        setup_db()
    else:
        setup_minidb()

def is_true(t):
    return t.lower() in ['y', 'yes', '1', 'true', 't', 'ok']


def input_file_to_samples(input_file, sep='\t'):
    sample_report = []
    sample_df = pd.read_csv(input_file, sep=sep, dtype='str')
    sample_df.fillna('', inplace=True)

    if 'trim' not in sample_df.columns:
        sample_df['trim'] = False
    else:
        sample_df['trim'] = sample_df['trim'].apply(lambda x: is_true(x))

    if 'strain' not in sample_df.columns:
        sample_df['strain'] = None

    if 'metadata' not in sample_df.columns:
        sample_df['metadata'] = ''
    if 'gsize' not in sample_df.columns:
        sample_df['gsize'] = None
    # 1. validate input
    for i, row in sample_df.iterrows():
        sample_id = row['sample_id']
        if not valid_id(sample_id):
            raise Exception(
                '{} invalid: sample ID can only contain alpha-numerical or underscore charactors'.format(
                    sample_id))
        input_files = row['files'].split(';')
        input_files = [input_file.strip() for input_file in input_files]
        if len(input_files) <= 0:
            raise Exception(
                'No input file for sample {}'.format(sample_id))

        for input_file in input_files:
            if not os.path.isfile(input_file):
                raise Exception(
                    'Input file {} (sample {}) not found!'.format(input_file, sample_id))

        metadata = row['metadata'].split(';')
        mt = {}
        if len(metadata) > 0:
            for kv in metadata:
                if len(kv.split(':')) == 2:
                    k, v = kv.split(':')
                    mt[k] = v
        sample = {
            'id': sample_id,
            'name': row['sample_desc'].strip(),
            'input_type': row['input_type'].strip(),
            'files': ';'.join(input_files),  # Re-join to make sure no white characters slipped in
            'genus': row['genus'].strip(),
            'species': row['species'].strip(),
            'strain': row['strain'].strip(),
            'trim': row['trim'],
            'metadata': mt,
            'updated': False,
            'gsize':int(row['gsize']) if row['gsize'] else None,
        }
        sample_report.append(sample)
    return sample_report


def single_genome_clean_func(args):
    """
    Remove a single genome from the file structure, report any collections that might be affected
    Parameters
    ----------
    args

    Returns
    -------

    """
    pass


def collection_clean_func(args):
    pass


def single_genome_analysis_func(args):
    """
    Parse command line and call single genome analysis
    Parameters
    ----------
    args

    Returns
    -------

    """
    work_dir = args.work_dir
    threads = args.threads
    memory = args.memory
    timing_log = args.time_log
    overwrite=args.overwrite
    if threads <= 0:
        threads = multiprocessing.cpu_count()
    #auto setup db if db not exists
    if not os.path.exists('db') :
        if args.initdb:
            setup_db()
        else:
            setup_minidb()
    # run single sample pipeline
    samples = input_file_to_samples(args.input)
    single_genome_analysis(samples, work_dir, overwrite, threads, memory, timing_log)
    return samples


def pan_genome_analysis_func(args):
    """
    Parse command line and call single genome analysis and pan-genome analysis
    """
    collection_id = args.collection_id
    collection_name = args.collection_name
    if not collection_name:
        collection_name = collection_id

    work_dir = args.work_dir
    threads = args.threads
    memory = args.memory
    timing_log = args.time_log
    pangenome_method = args.pangenome_method
    if not pangenome_method:
        pangenome_method='panta'
    # if pangenome_method == 'roary':
    #     raise Exception('roary is no longer supported due to incompatibility of software installation. Please use panta instead!')
    #overwrite = False
    overwrite=args.overwrite

    if not valid_id(collection_id):
        raise Exception('{} invalid: collection ID can only contain alpha-numerical or underscore charactors'.format(collection_id))

    if threads <= 0:
        threads = multiprocessing.cpu_count()
    #auto setup db if db not exists
    if not os.path.exists('db') :
        if args.initdb:
            setup_db()
        else:
            setup_minidb()
    samples = input_file_to_samples(args.input)

    # First run single analysis
    samples = single_genome_analysis(samples, work_dir, assembly_method=args.assembly_method, overwrite=overwrite, threads=threads, memory=memory, timing_log=timing_log)
    report = pan_genome_analysis(
        samples, work_dir,
        collection_id, collection_name, overwrite=overwrite,
        threads=threads, memory=memory, timing_log=timing_log,method=pangenome_method,rate_coverage=args.ratio_coverage, genetree=args.genetree, progressive=args.progressive,tree=args.tree_method)
    logger.info('Congratulations, collection {} is done!'.format(collection_id))
def main(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        prog='amromics',
        description='Tool for managing and analyzing antibiotic resistant bacterial datasets')
    parser.add_argument('-V', '--version', action='version', version=__version__)

    subparsers = parser.add_subparsers(title='sub command', help='sub command help')
    version_cmd = subparsers.add_parser(
        'dep', description='Check dependency versions',
        help='Print versions dependencies if exist',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    version_cmd.set_defaults(func=version_func)

    pg_cmd = subparsers.add_parser(
        'pg',
        description='Pan-genome analysis of a collection',
        help='Pan-genome analysis of a collection',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pg_cmd.set_defaults(func=pan_genome_analysis_func)
    pg_cmd.add_argument('-t', '--threads', help='Number of threads to use, 0 for all', default=0, type=int)
    pg_cmd.add_argument('-m', '--memory', help='Amount of memory in Gb to use', default=30, type=float)
    pg_cmd.add_argument('-c', '--collection-id', help='Collection ID', required=True, type=str)
    pg_cmd.add_argument('-n', '--collection-name', help='Collection name', type=str, default='')
    pg_cmd.add_argument('-i', '--input', help='Input file', required=True, type=str)
    pg_cmd.add_argument('--work-dir', help='Working directory', default='data/work')
    pg_cmd.add_argument('--time-log', help='Time log file', default=None, type=str)

    pg_cmd.add_argument('--pangenome-method', choices=['panta', 'roary'], default='panta', help='Pangenome method')
    pg_cmd.add_argument('--assembly-method', choices=['spades', 'skesa'], default='skesa', help='Short read assembly methods')
    pg_cmd.add_argument('--tree-method', choices=['fasttree', 'iqtree'], default='fasttree', help='Tree building method')
    pg_cmd.add_argument('-r', '--ratio-coverage', help='Ratio of coverage to align', default=0.25, type=float)
    pg_cmd.add_argument('--genetree', help='Run phylogenty for each gene cluster or not', default=False)
    pg_cmd.add_argument('--progressive', help='Run pangenome in progressive mode', default=False)
    pg_cmd.add_argument('--overwrite', help='Force overwrite exist results', default=False)
    pg_cmd.add_argument('--initdb', help='Init full database', required=False,type=eval,choices=[True, False],default='False')

    sg_cmd = subparsers.add_parser(
        'sg',
        description='Single-genome analysis',
        help='Single-genome analysis of a collection',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    sg_cmd.set_defaults(func=single_genome_analysis_func)
    sg_cmd.add_argument('-t', '--threads', help='Number of threads to use, 0 for all', default=0, type=int)
    sg_cmd.add_argument('-m', '--memory', help='Amount of memory in Gb to use', default=30, type=float)
    sg_cmd.add_argument('-i', '--input', help='Input file', required=True, type=str)
    sg_cmd.add_argument('--work-dir', help='Working directory', default='data/work')
    sg_cmd.add_argument('--overwrite', help='Force overwrite exist results', default=False)
    sg_cmd.add_argument('--time-log', help='Time log file', default=None, type=str)
    sg_cmd.add_argument('--initdb', help='Init full database', required=False,type=eval,choices=[True, False],default='False')

    db_cmd = subparsers.add_parser(
        'download_db',
        description='Force download and setup databases',
        help='Force download and setup databases',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    db_cmd.add_argument('--initdb', help='Init full database', required=False,type=eval,choices=[True, False],default='False')
    db_cmd.set_defaults(func=setup_db_func)


    args = parser.parse_args(arguments)
    return args.func(args)


if __name__ == "__main__":
    main()
