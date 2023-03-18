import os,shutil
import logging
import multiprocessing
from Bio import SeqIO
import csv
import pandas as pd
import json
import gzip
from panta import *
from amromics.utils.command import run_command
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def run_roary(gff_folder,overwrite=False,threads=0, base_dir='.', timing_log=None):
    """
        Run roay make pangeome analysis (using prokka results in previous step)
        :param read_data: result holder
        :param base_dir: working directory
        :param threads: number of core CPU
        :return:
    """

    roary_folder=os.path.join(base_dir,'pangenome/roary')
    temp_folder = os.path.join(base_dir, 'pangenome/temp_roary')
    roary_output = os.path.join(roary_folder, 'summary_statistics.txt')
    if os.path.isfile(roary_output) and (not overwrite):
        logger.info('roary has run and the input has not changed, skip roarying')
        return roary_folder
    if not os.path.isdir(temp_folder):
        os.makedirs(temp_folder)
    gff_list = []
    for filename in os.listdir(gff_folder):
        if filename.endswith('.gz'):
            sample_id = filename.replace('.gff.gz','')
            #gffgz_file = os.path.join(sample['annotation'], sample_id + '.gff.gz')
            gff_file = os.path.join(temp_folder, sample_id + '.gff')
            if run_command('gunzip -c {} > {}'.format(os.path.join(gff_folder, filename), gff_file)) != 0:
                raise Exception('Cannot get {}'.format(os.path.join(gff_folder, filename)))
            gff_list.append(gff_file)
        else:
            gff_list.append(os.path.join(gff_folder, filename))

    # Make sure the directory is not there or roary will add timestamp
    if os.path.isfile(roary_folder):
        os.remove(roary_folder)
    if os.path.exists(roary_folder):
        shutil.rmtree(roary_folder)
    cmd = 'roary -p {} -f {} -v '.format(threads, roary_folder) + ' '.join(gff_list)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('roary fail to run!')

    cmd = 'gzip ' + os.path.join(roary_folder, 'gene_presence_absence.csv')
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running {}'.format(cmd))

    shutil.rmtree(temp_folder)

    return roary_folder
def run_panta(gff_folder,overwrite=False,threads=0, base_dir='.', timing_log=None):
    starttime = datetime.now()
    out_folder=os.path.join(base_dir,'pangenome/panta')
    temp_dir = os.path.join(base_dir, 'pangenome/temp_panta')

    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv.gz')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv.gz')
    gff_list = []
    samples=[]
    for filename in os.listdir(gff_folder):
        if filename.endswith('.gz'):
            sample_id = filename.replace('.gff.gz','')
            #gffgz_file = os.path.join(sample['annotation'], sample_id + '.gff.gz')
            gff_file = os.path.join(temp_folder, sample_id + '.gff')
            if run_command('gunzip -c {} > {}'.format(os.path.join(gff_folder, filename), gff_file)) != 0:
                raise Exception('Cannot get {}'.format(os.path.join(gff_folder, filename)))
            gff_list.append(gff_file)
        else:
            gff_list.append(os.path.join(gff_folder, filename))

    for gff in gff_list:
        base_name = os.path.basename(gff)
        if gff.endswith('.gff'):
            sample_id = base_name[:-4]
        elif gff.endswith('.gff.gz'):
            sample_id = base_name[:-7]
        else:
            raise Exception(f'{gff} should be a gff3 file')

        sample_id = sample_id.replace('-','_')#Make sure that - is not part of sample_id
        sample_id_list.append(sample_id)
        samples.append({'id':sample_id, 'gff_file':gff, 'assembly':None})
    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_folder,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    combined_faa = data_preparation.combine_proteins(
        out_dir=out_folder,
        samples=samples)

    # main_pipeline
    cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
        faa_file=combined_faa,
        out_dir=temp_dir,
        threads=threads)
    # logger.info(f'len cd_hit_clusters = {len(cd_hit_clusters)}')


    #print(f'Diamond = {args.diamond}')
    blast_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = cd_hit_represent_fasta,
        query_fasta = cd_hit_represent_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)
    #blast_result = os.path.join(os.path.join(temp_dir, 'blast'), 'blast_results')
    #TODO check here

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result)

    inflated_clusters, clusters = main_pipeline.reinflate_clusters(
        cd_hit_clusters=cd_hit_clusters,
        mcl_file=mcl_file)
    logger.info(f'len inflated_clusters = {len(inflated_clusters)} len clusters = {len(clusters)}')


    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len split_clusters = {len(split_clusters)}')

    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)


    output.create_outputs(annotated_clusters,samples,out_folder)
    # if args.alignment != None:
    #     post_analysis.run_gene_alignment(annotated_clusters, samples, out_folder, args.alignment, threads)

    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_folder, 'gene_annotation.csv.gz')
    main_gene_position_fn = os.path.join(out_folder, 'gene_position.csv.gz')

    shutil.copy(gene_annotation_fn, main_gene_annotation_fn)
    shutil.copy(gene_position_fn, main_gene_position_fn)

    json.dump(samples, open(os.path.join(out_folder, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.copy(cd_hit_represent_fasta, os.path.join(out_folder, 'representative.fasta'))
    json.dump(clusters, open(os.path.join(out_folder, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    cmd = f'gzip -c {blast_result} > ' + os.path.join(out_folder, 'blast.tsv.gz')
    os.system(cmd)


    # shutil.rmtree(temp_dir)

    elapsed = datetime.now() - starttime
    logger.info(f'Done -- time taken {str(elapsed)}')
