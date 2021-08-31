import argparse
import os
import re
import shutil
import logging
import json
from glob import glob
from datetime import datetime
from pan_genome import *

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)s : %(message)s',
    datefmt='%I:%M:%S')
logger = logging.getLogger(__name__)


def run_post_analysis(gene_annotation,gene_position,inflated_clusters,dontsplit,samples, out_dir):
    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_annotation=gene_annotation,
        gene_position=gene_position,
        unsplit_clusters= inflated_clusters,
        dontsplit=dontsplit
        )
    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters, 
        gene_annotation=gene_annotation)

    # output
    spreadsheet_file = output.create_spreadsheet(
        annotated_clusters=annotated_clusters, 
        gene_annotation=gene_annotation,
        samples=samples,
        out_dir=out_dir
    )
    rtab_file = output.create_rtab(
        annotated_clusters=annotated_clusters, 
        gene_annotation=gene_annotation,
        samples=samples,
        out_dir=out_dir
    )
    summary_file = output.create_summary(
        rtab_file=rtab_file, 
        out_dir=out_dir
    )


def collect_sample(sample_id_list, args):
    samples = []
    fasta_ext = ('.fasta', '.fna', 'ffn')
    for path in args.gff_files:
        if not path.endswith('gff'):
            raise Exception(f'{path} should be a gff3 file')
        sample = {'gff_file':path}
        dir_name = os.path.dirname(path)
        base_name = os.path.basename(path)
        sample_id = base_name.rsplit('.', 1)[0]
        if args.fasta == True:
            try:
                fasta_file = [f for f in glob(dir_name+'/'+sample_id+'*') if f.endswith(fasta_ext)][0]
                sample['fasta_file'] = fasta_file
            except:
                raise Exception(f'The corresponding fasta file of {sample_id} does not exist')
        sample_id = re.sub(r'\W', '_', sample_id)
        if sample_id in sample_id_list:
            logging.info(f'{sample_id} already exists -- skip')
            continue
        else:
            sample_id_list.append(sample_id)
        sample['id'] = sample_id
        samples.append(sample)
    samples.sort(key= lambda x:x['id'])
    return samples

def run_main_pipeline(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    dontsplit = args.dont_split
    diamond = args.diamond
    identity = args.identity
    
    temp_dir = os.path.join(out_dir, 'temp')
    timing_log = os.path.join(out_dir, 'time.log')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)    

    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')
    
    # data preparation
    gene_annotation = {}
    gene_position = {}
    data_preparation.extract_proteins(
        samples=samples,
        out_dir=temp_dir,
        gene_annotation = gene_annotation,
        gene_position = gene_position,
        fasta=args.fasta,
        table=args.table,
        threads=threads)

    combined_faa = data_preparation.combine_proteins(
        out_dir=temp_dir,
        samples=samples,
        timing_log=timing_log)

    # main_pipeline
    cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
        faa_file=combined_faa,
        out_dir=temp_dir,
        threads=threads, 
        timing_log=timing_log)

    blast_result = main_pipeline.pairwise_alignment(
        diamond=diamond,
        database_fasta = cd_hit_represent_fasta,
        query_fasta = cd_hit_represent_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        identity=identity,
        threads=threads, 
        timing_log=timing_log)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = blast_result,
        threads=threads,
        timing_log=timing_log)

    inflated_clusters, clusters = main_pipeline.reinflate_clusters(
        cd_hit_clusters=cd_hit_clusters,
        mcl_file=mcl_file)
    
    # post analysis
    run_post_analysis(gene_annotation, gene_position, inflated_clusters, dontsplit, samples, out_dir)

    # output for next run
    output.export_gene_annotation(gene_annotation, out_dir)
    json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.copy(cd_hit_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))

    # shutil.rmtree(temp_dir)
    shutil.rmtree(os.path.join(temp_dir, 'samples'))
    
    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')

    

def run_add_sample_pipeline(args):
    starttime = datetime.now()

    collection_dir = args.collection_dir
    out_dir = args.outdir
    if out_dir == None:
        out_dir = collection_dir
    threads = args.threads
    dontsplit = args.dont_split
    diamond = args.diamond
    identity = args.identity

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    temp_dir = os.path.join(out_dir, 'temp')
    timing_log = os.path.join(out_dir, 'time.log')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
    else:
        os.makedirs(temp_dir)
    
    # Check required files
    gene_annotation_file = os.path.join(collection_dir, 'gene_annotation.tsv')
    if not os.path.isfile(gene_annotation_file):
        raise Exception(f'{gene_annotation_file} does not exist')
    gene_annotation = output.import_gene_annotation(gene_annotation_file)
    
    gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))
    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    old_clusters = json.load(open(os.path.join(collection_dir, 'clusters.json'), 'r'))

    old_represent_faa = os.path.join(collection_dir, 'representative.fasta')
    if not os.path.isfile(old_represent_faa):
        raise Exception(f'{old_represent_faa} does not exist')
    
    old_blast_result = os.path.join(collection_dir, 'blast.tsv')
    if not os.path.isfile(old_blast_result):
        raise Exception(f'{old_blast_result} does not exist')

    # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')
    
    # data preparation
    data_preparation.extract_proteins(
        samples=new_samples,
        out_dir=temp_dir,
        gene_annotation = gene_annotation,
        gene_position = gene_position,
        fasta=args.fasta,
        table=args.table,
        threads=threads
        )
    new_combined_faa = data_preparation.combine_proteins(
        out_dir=temp_dir,
        samples=new_samples,
        timing_log=timing_log)

    not_match_faa, cd_hit_2d_clusters = add_sample_pipeline.run_cd_hit_2d(
        database_1 = old_represent_faa,
        database_2 = new_combined_faa,
        out_dir = temp_dir,
        threads=threads, 
        timing_log=timing_log)

    not_match_represent_faa, not_match_clusters = main_pipeline.run_cd_hit(
        faa_file=not_match_faa,
        out_dir=temp_dir,
        threads=threads, 
        timing_log=timing_log)

    blast_1_result = main_pipeline.pairwise_alignment(
        diamond=diamond,
        database_fasta = old_represent_faa,
        query_fasta = not_match_represent_faa,
        out_dir = os.path.join(temp_dir, 'blast1'),
        identity=identity,
        threads=threads, 
        timing_log=timing_log
        )

    blast_2_result = main_pipeline.pairwise_alignment(
        diamond=diamond,
        database_fasta = not_match_represent_faa,
        query_fasta = not_match_represent_faa,
        out_dir = os.path.join(temp_dir, 'blast2'),
        identity=identity,
        threads=threads, 
        timing_log=timing_log
        )

    combined_blast_result = add_sample_pipeline.combine_blast_results(
        blast_1=old_blast_result, 
        blast_2=blast_1_result, 
        blast_3=blast_2_result, 
        outdir=temp_dir)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = combined_blast_result,
        threads=threads, 
        timing_log=timing_log)

    inflated_clusters, new_clusters = add_sample_pipeline.reinflate_clusters(
        old_clusters=old_clusters, 
        cd_hit_2d_clusters=cd_hit_2d_clusters,
        not_match_clusters=not_match_clusters,
        mcl_file=mcl_file
        )

    # post analysis
    new_samples.extend(old_samples)
    new_samples.sort(key= lambda x:x['id'])
    
    run_post_analysis(gene_annotation, gene_position, inflated_clusters, dontsplit, new_samples, out_dir)

    # output for next run
    output.export_gene_annotation(gene_annotation, out_dir)
    json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    json.dump(new_samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    add_sample_pipeline.combine_representative(not_match_represent_faa, old_represent_faa, out_dir)
    json.dump(new_clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    shutil.copy(combined_blast_result, os.path.join(out_dir, 'blast.tsv'))

    # shutil.rmtree(temp_dir)
    shutil.rmtree(os.path.join(temp_dir, 'samples'))

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    
    main_cmd = subparsers.add_parser(
        'main',
        description='Main pipeline: run pan-genome analysis for the first time',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    main_cmd.set_defaults(func=run_main_pipeline)
    main_cmd.add_argument('gff_files', help='a.gff b.gff ... (*.gff)', type=str, nargs='+')
    main_cmd.add_argument('-o', '--outdir', help='output directory', required=True, type=str)
    main_cmd.add_argument('-s', '--dont-split', help='dont split paralog clusters', default=False, action='store_true')
    main_cmd.add_argument('-d', '--diamond', help='use Diamond for all-agaist-all alignment instead of Blastp', default=False, action='store_true')
    main_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=95, type=float)
    main_cmd.add_argument('-f', '--fasta', help='fasta files are seperated from gff files. (fasta file must have the same name, be in the same folder of coresponding gff file, and have one of following extension: .fasta .fna .fnn)', default=False, action='store_true')
    main_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    main_cmd.add_argument('--table', help='codon table', default=11, type=int)

    add_cmd = subparsers.add_parser(
        'add',
        description='Add pipeline: add sample into previous collection',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_cmd.set_defaults(func=run_add_sample_pipeline)
    add_cmd.add_argument('gff_files', help='a.gff b.gff ... (*.gff)', type=str, nargs='+')
    add_cmd.add_argument('-c', '--collection-dir', help='previous collection directory', required=True, type=str)
    add_cmd.add_argument('-o', '--outdir', help='output directory',default=None, type=str)
    add_cmd.add_argument('-s', '--dont-split', help='dont split paralog clusters', default=False, action='store_true')
    add_cmd.add_argument('-d', '--diamond', help='use Diamond for all-agaist-all alignment instead of Blastp', default=False, action='store_true')
    add_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=95, type=float)
    add_cmd.add_argument('-f', '--fasta', help='fasta files are seperated from gff files. (fasta file must have the same name, be in the same folder of coresponding gff file, and have one of following extension: .fasta .fna .fnn)', default=False, action='store_true')
    add_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    add_cmd.add_argument('--table', help='codon table', default=11, type=int)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()