# -*- coding: utf-8 -*-
import json
import os
import shutil
import re
import glob
import gzip
import datetime
import logging
import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from amromics.utils import get_open_func, get_compress_type, translate_dna, run_command

logger = logging.getLogger(__name__)

# TODOs:
# - Can make it faster with using fasttree (parsnps need to have this option specifically set
# - By default, parsnp use bootstrap of 1000. See if we can change the value and get the boottrap values


def assemble_shovill(sample, sample_dir, threads=4, memory=50, overwrite=False, timing_log=None):
    """
    Run assembly using shovill

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    memory: float
        maximum memory used
    timing_log: str
        log file
    Returns str
    -------
        Path to the assembly
    """
    sample_id = sample['id']
    path_out = os.path.join(sample_dir, 'assembly')
    assembly_file = os.path.join(path_out, sample_id + '_contigs.fasta.gz')

    if not os.path.exists(path_out):
        os.makedirs(path_out)
    elif os.path.isfile(assembly_file) and (not overwrite):
        return assembly_file

    cmd = 'shovill  --ram {memory} --cpus {threads} --outdir {path_out}'.format(
        memory=int(memory), threads=threads, path_out=path_out)

    if 'trim' in sample and sample['trim']:
        cmd += ' --trim --depth 200'
    else:
        cmd += ' --depth 120'

    pe_files = sample['files'].split(';')
    if len(pe_files) > 1:
        pe1 = pe_files[0]
        pe2 = pe_files[1]
        cmd += ' --force  --R1 {pe1} --R2 {pe2}'.format(pe1=pe1, pe2=pe2)
    else:
        raise Exception('Only support pair-end reads!')

    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running shovill!')

    # Read in list of contigs
    contigs = list(SeqIO.parse(os.path.join(path_out, 'contigs.fa'), "fasta"))
    contigs = sorted(contigs, key=len, reverse=True)
    with gzip.open(assembly_file, 'wt') as f:
        for i, contig in enumerate(contigs):
            contig.id = sample['id']+'_C'+str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")

    run_command('rm -f ' + os.path.join(path_out, 'spades.fasta'))
    run_command('rm -f ' + os.path.join(path_out, 'contigs.fa'))
    run_command('rm -f ' + os.path.join(path_out, 'contigs.gfa'))
    run_command('gzip ' + os.path.join(path_out, 'shovill.log'))

    sample['updated'] = True
    return assembly_file


def get_assembly(sample, sample_dir, overwrite=False):
    """
    Get the assembly from user input

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    Returns
    -------
        path to assembly file
    """
    path_out = os.path.join(sample_dir, 'assembly')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    open_func = get_open_func(sample['files'])
    with open_func(sample['files'], 'rt') as fn:
        contigs = list(SeqIO.parse(fn, 'fasta'))

    assembly_file = os.path.join(path_out, sample['id'] + '_contigs.fasta.gz')
    if os.path.isfile(assembly_file) and (not overwrite):
        logger.info(f'Not overwrite {assembly_file}')
        return assembly_file

    contigs = sorted(contigs, key=len, reverse=True)
    with gzip.open(assembly_file, 'wt') as f:
        for i, contig in enumerate(contigs):
            contig.id = sample['id'] + '_C' + str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")

    sample['updated'] = True
    return assembly_file


def annotate_prokka(sample, sample_dir,  threads=8, overwrite=False, timing_log=None):
    """
    Annotate the sample using prokka

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing annotation
    timing_log: str
        log file
    Returns
    -------
        path to prokka folder
    """

    path_out = os.path.join(sample_dir, 'prokka')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    gff_file_out = os.path.join(path_out, sample['id'] + '.gff.gz')
    gbk_file_out = os.path.join(path_out, sample['id'] + '.gbk.gz')

    if os.path.isfile(gff_file_out) and os.path.isfile(gbk_file_out) and (not overwrite):
        # Dont run again if gff/gbk file exists
        logger.info('GFF and GBK files found, skip annotating')
        return path_out

    gunzip_fasta = os.path.join(path_out, sample['id'] + '.fin')
    cmd = 'gunzip -c {} > {}'.format(sample['assembly'], gunzip_fasta)
    run_command(cmd)
    cmd = 'prokka --force --cpus {threads} --addgenes --mincontiglen 200'.format(threads=threads)
    cmd += ' --prefix {sample_id} --locus {sample_id} --outdir {path} '.format(sample_id=sample['id'], path=path_out)
    if sample['genus']:
        cmd += ' --usegenus --genus ' + sample['genus']
    if sample['species']:
        cmd += ' --species ' + sample['species']
    if sample['strain']:
        cmd += ' --strain ' + sample['strain']

    # Disable this for now so that we dont have to install signalp
    # if sample['gram']:
    #    cmd += ' --gram ' + sample['gram']
    cmd += ' ' + gunzip_fasta
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Command {} returns non-zero ()!'.format(cmd, ret))

    for file_name in glob.glob(os.path.join(path_out, '*')):
        ext = file_name[-3:]
        if ext in ['gff', 'gbk', 'ffn']: # fna?
            run_command('gzip {}'.format(file_name))
        else:
            os.remove(file_name)

    # for ext in ['err', 'faa', 'fsa', 'log', 'sqn', 'tbl', 'tsv', 'txt']:
    #     # What about ffn and fna?
    #     file_name = os.path.join(path_out, sample['id'] + '.' + ext)
    #     if os.path.isfile(file_name):
    #         os.remove(file_name)

    sample['updated'] = True
    return path_out


def mlst(sample, sample_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run mlst

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing result
    timing_log: str
        log file
    Returns
    -------
        path to mlst result file
    """
    path_out = os.path.join(sample_dir, 'mlst')
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    mlst_out = os.path.join(path_out, sample['id'] + '_mlst.tsv')
    if os.path.isfile(mlst_out) and (not overwrite):
        logger.info('MLST for {} exists, skip mlsting'.format(sample['id']))
        return mlst_out

    cmd = 'mlst --quiet --threads {threads} --nopath {infile} > {outfile}'.format(
        threads=threads,
        infile=sample['assembly'],
        outfile=mlst_out)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running mlst')
    sample['updated'] = True
    return mlst_out


def detect_amr(sample, sample_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run abricate to identify resistant genes

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing result
    timing_log: str
        log file
    Returns
    -------
        path to resistant gene file
    """
    path_out = os.path.join(sample_dir, 'abricate')
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    # TODO: replace by consensus db later
    amr_out = os.path.join(path_out, sample['id'] + '_resistome.tsv')
    if os.path.isfile(amr_out) and (not overwrite):
        logger.info('Resistome for {} exists, skip analysis'.format(sample['id']))
        return amr_out

    cmd = 'abricate --quiet --threads {threads} --nopath --db card {infile} > {outfile}'.format(
        threads=threads,
        infile=sample['assembly'],
        outfile=amr_out)
    if run_command(cmd, timing_log) != 0:
        raise Exception('Error running amr')
    sample['updated'] = True
    return amr_out


def detect_virulome(sample, sample_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run abricate to identify virulent genes using VFDB

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing result
    timing_log: str
        log file
    Returns
    -------
        path to virulent gene file
    """
    path_out = os.path.join(sample_dir, 'abricate')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    vir_out = os.path.join(path_out, sample['id'] + '_virulome.tsv')
    if os.path.isfile(vir_out) and (not overwrite):
        logger.info('Virulome for {} exists, skip analysis'.format(sample['id']))
        return vir_out

    cmd = 'abricate --quiet --threads {threads} --nopath --db vfdb {infile} > {outfile}'.format(
        threads=threads,
        infile=sample['assembly'],
        outfile=vir_out)
    if run_command(cmd, timing_log) != 0:
        raise Exception('Error running virulent detection')
    sample['updated'] = True
    return vir_out


def detect_plasmid(sample, sample_dir,  threads=8, overwrite=False, timing_log=None):
    """
    Detect plasmids

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing result
    timing_log: str
        log file
    Returns
    -------
        path to origin of replication file
    """

    path_out = os.path.join(sample_dir, 'abricate')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    # Plasmid finder
    oriREP_out = os.path.join(path_out, sample['id'] + '_plasmid.tsv')
    if os.path.isfile(oriREP_out) and (not overwrite):
        logger.info('ORI for {} exists, skip analysis'.format(sample['id']))
        return oriREP_out

    cmd = 'abricate --quiet --threads {threads} --nopath --db plasmidfinder {infile} > {outfile}'.format(
        threads=threads, infile=sample['assembly'], outfile=oriREP_out)
    if run_command(cmd, timing_log) != 0:
        return None

    sample['updated'] = True
    return oriREP_out


def run_single_sample(sample, sample_dir, threads=8, memory=50, timing_log=None):
    """
    Run the pipeline to analyse one single sample. Steps:
        0. Assembly if the input are read data
        1. Annotate the assembly using prokka
        2. Identify resistant, virulent genes, origins of replication, and plasmids

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    memory: float
        maximum memory
    timing_log: str
        log file
    Returns
    -------
        the sample object
    """

    sample['execution_start'] = str(datetime.datetime.now())

    if sample['input_type'] not in ['asm', 'assembly']:
        sample['assembly'] = assemble_shovill(
            sample, sample_dir=sample_dir, threads=threads, memory=memory, timing_log=timing_log)
    else:
        sample['assembly'] = get_assembly(sample, sample_dir=sample_dir)

    sample['annotation'] = annotate_prokka(
        sample, sample_dir=sample_dir,  threads=threads, timing_log=timing_log)
    sample['mlst'] = mlst(sample, sample_dir=sample_dir, threads=threads)
    sample['resistome'] = detect_amr(
        sample, sample_dir=sample_dir, threads=threads, timing_log=timing_log)
    sample['virulome'] = detect_virulome(
        sample, sample_dir=sample_dir, threads=threads, timing_log=timing_log)
    sample['plasmid'] = detect_plasmid(
        sample, sample_dir=sample_dir, threads=threads, timing_log=timing_log)
    sample['execution_end'] = str(datetime.datetime.now())
    return sample


def run_roary(report, collection_dir='.', threads=8, overwrite=False, timing_log=None):
    """
    Run roary for pangenome analysis. If the list of samples has not changed, and
    none of the samples has changed, the existing tree will be kept unless overwrite is
    set to True

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """

    for sample in report['samples']:
        # Check if any sample has been updated
        overwrite = overwrite or sample['updated']
    # any([overwrite] + [sample['updated'] for sample in report['samples']])

    roary_folder = os.path.join(collection_dir, 'roary')
    temp_folder = os.path.join(collection_dir, 'temp_roary')
    roary_output = os.path.join(roary_folder, 'summary_statistics.txt')

    # Check if roary has run for the same dataset ID and the same set of samples
    report['pan_genome'] = roary_folder
    if os.path.isfile(roary_output) and (not overwrite):
        logger.info('roary has run and the input has not changed, skip roarying')
        return report

    if not os.path.isdir(temp_folder):
        os.makedirs(temp_folder)

    gff_list = []
    for sample in report['samples']:
        sample_id = sample['id']
        gffgz_file = os.path.join(sample['annotation'], sample_id + '.gff.gz')
        gff_file = os.path.join(temp_folder, sample_id + '.gff')
        if run_command('gunzip -c {} > {}'.format(gffgz_file, gff_file)) != 0:
            raise Exception('Cannot get {}'.format(gffgz_file))
        gff_list.append(gff_file)

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
    return report


def run_phylogeny(report, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run parsnp to create phylogeny tree. If the list of samples has not changed, and
    none of the samples has changed, the existing tree will be kept unless overwrite is
    set to True

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    phylogeny_folder = os.path.join(collection_dir, 'phylogeny')
    if not os.path.exists(phylogeny_folder):
        os.makedirs(phylogeny_folder)
    genome_dir = os.path.join(collection_dir, 'temp/fasta')
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)

    report['phylogeny'] = phylogeny_folder
    phylogeny_file = os.path.join(phylogeny_folder, 'parsnp.tree')
    if os.path.isfile(phylogeny_file) and (not overwrite):
        logger.info('phylogeny tree exists and input has not changed, skip phylogeny analysis')
        return report

    temp_folder = os.path.join(collection_dir, 'temp_phylo')
    if not os.path.isdir(temp_folder):
        os.makedirs(temp_folder)
    reference_genome = None
    sample_list = []
    for i, sample in enumerate(report['samples']):
        fasta_file = os.path.join(temp_folder, sample['id'] + '.fasta')
        cmd = 'gunzip -c {} > {}'.format(sample['assembly'], fasta_file)
        run_command(cmd)
        if i == 0:
            reference_genome = fasta_file
        else:
            sample_list.append(fasta_file)
    cmd = 'parsnp -r {} -d {} -o {} -p {}'.format(
        reference_genome,
        ' '.join(sample_list),
        phylogeny_folder, threads)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running parsnp')
    run_command('gzip {}'.format(os.path.join(phylogeny_folder, 'parsnp.xmfa')))
    run_command('gzip {}'.format(os.path.join(phylogeny_folder, 'parsnp.ggr')))
    shutil.rmtree(temp_folder)
    return report


def run_alignment(report, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run phylogenetic analysis of gene clusters. If the list of samples has not changed, and
    none of the samples has changed, the existing tree will be kept unless overwrite is
    set to True

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    gene_cluster_file = report['pan_genome'] + '/gene_presence_absence.csv.gz'
    dict_cds = {}
    for sample in report['samples']:
        with gzip.open(os.path.join(sample['annotation'], sample['id'] + '.ffn.gz'), 'rt') as fn:
            for seq in SeqIO.parse(fn, 'fasta'):
                dict_cds[seq.id] = seq

    # make folder contains sequences for each gene
    alignment_dir = os.path.join(collection_dir, 'alignments')
    gene_df = pd.read_csv(gene_cluster_file, dtype=str, compression='gzip')
    gene_df.fillna('', inplace=True)

    sample_columns = list(gene_df.columns)[14:]
    for _, row in gene_df.iterrows():
        gene_id = row['Gene']
        gene_list = []
        for sample_column in sample_columns:
            if row[sample_column]:
                # roary can pool together genes from the same sample and tab-separate them
                for sample_gene in row[sample_column].split('\t'):
                    gene_list.append(sample_gene)
                    # TODO: make sure all samples in this gene have not updated

        gene_list = sorted(gene_list)
        # Only analyse if there are at least 3 genes
        if len(gene_list) < 3:
            logger.info('There are too few genes for {} skipping'.format(gene_id))
            continue

        gene_dir = os.path.join(alignment_dir, gene_id)
        # Check if done before
        gene_list_json = os.path.join(gene_dir, 'gene_list.json')
        # if os.path.isfile(os.path.join(gene_dir, 'parsnp.tree')) and (not overwrite):
        if not overwrite:
            if os.path.isfile(gene_list_json):
                with open(gene_list_json) as fn:
                    existing_gene_list = json.load(fn)
                    if gene_list == existing_gene_list:
                        logger.info('Phylogeny for gene {} done, skipping'.format(gene_id))
                        continue  # for _, row

        gene_file_dir = os.path.join(gene_dir, 'files')
        if not os.path.exists(gene_file_dir):
            os.makedirs(gene_file_dir)

        gene_files = []
        for sample_gene in gene_list:
            gene_file = os.path.join(gene_file_dir, sample_gene + '.fasta')
            SeqIO.write(dict_cds[sample_gene], gene_file, 'fasta')
            gene_files.append(gene_file)

        # Use the first gene as the reference
        cmd = 'parsnp -d {} -r {} -o {} -p {}'.format(
            ' '.join(gene_files[1:]), gene_files[0], gene_dir, threads)
        ret = run_command(cmd, timing_log)
        # if ret != 0:
        #     raise Exception('error')

        with open(gene_list_json, 'w') as fn:
            json.dump(gene_list, fn)
        run_command('gzip {}'.format(os.path.join(gene_dir, 'parsnp.xmfa')))
        run_command('gzip {}'.format(os.path.join(gene_dir, 'parsnp.ggr')))

        if os.path.exists(gene_file_dir):
            shutil.rmtree(gene_file_dir)
        #clean up
        run_command('rm -f ' + os.path.join(gene_dir, '*.ini ') + os.path.join(gene_dir, '*block* '))
        shutil.rmtree(os.path.join(gene_dir, 'blocks'), True)
        shutil.rmtree(os.path.join(gene_dir, 'tmp'), True)

    report['alignments'] = alignment_dir
    return report


def run_species_phylogeny(report, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run iqtree to create phylogeny tree from core gene alignment. If the list of samples has 
    not changed, and none of the samples has changed, the existing tree will be kept unless 
    overwrite is set to True

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    phylogeny_folder = os.path.join(collection_dir, 'phylogeny')
    if not os.path.exists(phylogeny_folder):
        os.makedirs(phylogeny_folder)
    report['phylogeny'] = phylogeny_folder

    phylogeny_file = os.path.join(phylogeny_folder, 'core_gene_alignment.treefile')
    if os.path.isfile(phylogeny_file) and (not overwrite):
        logger.info('phylogeny tree exists and input has not changed, skip phylogeny analysis')
        return report

    aln_file = os.path.join(phylogeny_folder, 'core_gene_alignment.aln.gz')
    if not os.path.isfile(aln_file):
        aln_file = os.path.join(report['roary'], 'core_gene_alignment.aln.gz')

    cmd = 'iqtree -s {alignment} --prefix {prefix} -B 1000 -T {threads}'.format(
        alignment=aln_file, prefix=phylogeny_folder+'/core_gene_alignment', threads=threads)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('iqtree fail to create phylogeny tree from core gene alignment!')

    return report


def run_gene_phylogeny(report, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run phylogenetic analysis of gene clusters. If the list of samples has not changed, and
    none of the samples has changed, the existing tree will be kept unless overwrite is
    set to True

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    alignment_dir = os.path.join(collection_dir, 'alignments')
    if not os.path.exists(alignment_dir):
        os.makedirs(alignment_dir)
    gene_cluster_file = report['rtab']
    gene_df = pd.read_csv(gene_cluster_file, sep='\t', index_col='Gene')
    gene_df.fillna('', inplace=True)
    
    cmds_file = os.path.join(alignment_dir,"phylo_cmds")
    with open(cmds_file,'w') as cmds:
        for gene_id, row in gene_df.iterrows():
            # Only analyse if there are at least 3 genes
            if row.sum() < 3:
                continue
                
            gene_id = re.sub(r'\W+', '', gene_id)
            gene_dir = os.path.join(alignment_dir, gene_id)
            if not os.path.exists(gene_dir):
                os.makedirs(gene_dir)
            # check if done before
            iqtree_output = os.path.join(gene_dir, gene_id + '.treefile')
            if (not overwrite) and os.path.isfile(iqtree_output):
                continue

            gene_aln_file = os.path.join(gene_dir, gene_id + '.fna.aln')
            if not os.path.isfile(gene_aln_file):
                logger.info('{} does not exist'.format(gene_aln_file))
                continue

            cmd = f"iqtree -s {gene_aln_file} --prefix {gene_dir+'/'+gene_id} -m GTR -quiet -T 1 -B 1000 2> /dev/null"
            cmd += f" || iqtree -s {gene_aln_file} --prefix {gene_dir+'/'+gene_id} -m GTR -quiet -T 1"
            # translate to protein alignment
            #protein_aln_file = os.path.join(gene_dir, gene_id + '.faa.aln')
            #with open(protein_aln_file, 'w') as fh:
            #    for record in SeqIO.parse(gene_aln_file, 'fasta'):
            #        trans = translate_dna(str(record.seq))
            #        new_record = SeqRecord(Seq(trans), id=record.id,)
            #        SeqIO.write(new_record, fh, 'fasta')
            #cmd = f"iqtree -s {protein_aln_file} --prefix {gene_dir+'/'+gene_id} -m LG -quiet -T 1"
            #cmd = f"fasttree -nt -gtr -quiet {gene_aln_file} > {gene_dir+'/'+gene_id+'.treefile'} && echo '{gen_list_string}' > {gene_list_json}"
            cmds.write(cmd + '\n')
        
    cmd = f"parallel --bar -j {threads} -a {cmds_file}"
    ret = run_command(cmd, timing_log)
    report['alignments'] = alignment_dir
    return report


def get_gene_sequences(report, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Create protein sequences and nucleotide sequences for each gene cluster

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    logger.info('Getting sequences of gene clusters')
    gene_cluster_file = report['spreadsheet']
    dict_nucleotide = {}
    for sample in report['samples']:
        with gzip.open(os.path.join(sample['annotation'], sample['id'] + '.ffn.gz'), 'rt') as fn:
            for seq_record in SeqIO.parse(fn, 'fasta'):
                seq_record.seq = seq_record.seq[:-3]
                seq_record = SeqRecord(seq_record.seq, id=seq_record.id, description = '')
                dict_nucleotide[seq_record.id] = seq_record

    # make folder contains sequences for each gene
    alignment_dir = os.path.join(collection_dir, 'alignments')
    if not os.path.exists(alignment_dir):
        os.makedirs(alignment_dir)
    report['alignments'] = alignment_dir

    gene_df = pd.read_csv(gene_cluster_file, dtype=str, compression='gzip')
    gene_df.fillna('', inplace=True)
    sample_columns = list(gene_df.columns)[14:]
    for _, row in gene_df.iterrows():
        gene_id = row['Gene']
        gene_id = re.sub(r'\W+', '', gene_id)
        gene_dir = os.path.join(alignment_dir, gene_id)
        if not os.path.exists(gene_dir):
            os.makedirs(gene_dir)

        # check if done before
        protein_seq_file = os.path.join(gene_dir, gene_id + '.faa')
        nucleotide_seq_file = os.path.join(gene_dir, gene_id + '.fna')
        if (not overwrite) and os.path.isfile(protein_seq_file) and os.path.isfile(nucleotide_seq_file):
            continue

        gene_list = []
        for sample_column in sample_columns:
            if row[sample_column]:
                # roary can pool together genes from the same sample and tab-separate them
                for sample_gene in row[sample_column].split('\t'):
                    gene_list.append(sample_gene)
        gene_list = sorted(gene_list)

        with open(protein_seq_file, 'w') as prot_fh, open(nucleotide_seq_file, 'w') as nucl_fh:
            for sample_gene in gene_list:
                nu_seq_record = dict_nucleotide[sample_gene]
                SeqIO.write(nu_seq_record, nucl_fh, 'fasta')
                pro_seq = nu_seq_record.seq.translate(table=11)
                pro_seq_record = SeqRecord(pro_seq, id = nu_seq_record.id, description = '')
                SeqIO.write(pro_seq_record, prot_fh, 'fasta')

    return report


def run_protein_alignment(report, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Align protein sequence by mafft

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    alignment_dir = os.path.join(collection_dir, 'alignments')

    gene_cluster_file = report['rtab']   
    gene_df = pd.read_csv(gene_cluster_file, sep='\t', index_col='Gene')
    gene_df.fillna('', inplace=True)

    cmds_file = os.path.join(alignment_dir,"align_cmds")
    with open(cmds_file,'w') as cmds:
        for gene_id, row in gene_df.iterrows():
            # Only align if there are at least 2 sequences
            if row.sum() < 2:
                continue

            gene_id = re.sub(r'\W+', '', gene_id)
            gene_dir = os.path.join(alignment_dir, gene_id)

            # check if done before
            gene_aln_file = os.path.join(gene_dir, gene_id + '.faa.aln')
            if (not overwrite) and os.path.isfile(gene_aln_file):
                continue
            
            gene_seq_file = os.path.join(gene_dir, gene_id + '.faa')
            if not os.path.isfile(gene_seq_file):
                logger.info('{} does not exist'.format(gene_aln_file))
                continue
            
            cmd = f"mafft --auto --quiet --thread 1 {gene_seq_file} > {gene_aln_file}"
            cmds.write(cmd + '\n')
        
    cmd = f"parallel --bar -j {threads} -a {cmds_file}"
    ret = run_command(cmd, timing_log)
    report['alignments'] = alignment_dir
    return report


def create_nucleotide_alignment(report, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Create nucleotide alignment according to protein alignment

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    logger.info('Creating nucleotide alignment')
    alignment_dir = os.path.join(collection_dir, 'alignments')
    
    gene_cluster_file = report['rtab']  
    gene_df = pd.read_csv(gene_cluster_file, sep='\t', index_col='Gene')
    gene_df.fillna('', inplace=True)
    
    for gene_id, row in gene_df.iterrows():
        # Only run if there are at least 2 sequences
        if row.sum() < 2:
            continue
        gene_id = re.sub(r'\W+', '', gene_id)
        gene_dir = os.path.join(alignment_dir, gene_id)

        # check if done before
        nucleotide_aln_file = os.path.join(gene_dir, gene_id + '.fna.aln')
        if (not overwrite) and os.path.isfile(nucleotide_aln_file):
            continue

        protein_aln_file = os.path.join(gene_dir, gene_id + '.faa.aln')
        if not os.path.isfile(protein_aln_file):
            logger.info('{} does not exist'.format(protein_aln_file))
            continue
        protein_dict = {}
        for seq_record in SeqIO.parse(protein_aln_file, 'fasta'):
            protein_dict[seq_record.id] = str(seq_record.seq)

        nucleotide_seq_file = os.path.join(gene_dir, gene_id + '.fna')
        nucleotide_dict = {}
        for seq_record in SeqIO.parse(nucleotide_seq_file, 'fasta'):
            nucleotide_dict[seq_record.id] = str(seq_record.seq)

        with open(nucleotide_aln_file, 'w') as fh:
            for seq_id in protein_dict.keys():
                protein = protein_dict[seq_id]
                nucleotide = nucleotide_dict[seq_id]
                result = ''
                codon_pos = 0
                for c in protein:
                    if c == '-':
                        result += '---'
                    else:
                        result += nucleotide[codon_pos * 3: codon_pos * 3 + 3]
                        codon_pos += 1
                new_record = SeqRecord(Seq(result), id = seq_id, description = '')
                SeqIO.write(new_record, fh, 'fasta')
    return report


def create_core_gene_alignment(report, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Concatenate all the nucleotide alignment of core genes to create core gene alignment

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    logger.info('Creating core gene alignment')
    alignment_dir = os.path.join(collection_dir, 'alignments')
    phylogeny_folder = os.path.join(collection_dir, 'phylogeny')
    report['phylogeny'] = phylogeny_folder
    if not os.path.exists(phylogeny_folder):
        os.makedirs(phylogeny_folder)
    # check if done before
    core_gene_aln_file = os.path.join(phylogeny_folder, 'core_gene_alignment.aln.gz')
    if os.path.isfile(core_gene_aln_file) and (not overwrite):
        logger.info('Core gene alignment exists and input has not changed, skipping')
        return report

    gene_cluster_file = report['rtab']  
    gene_df = pd.read_csv(gene_cluster_file, sep='\t', index_col='Gene')
    gene_df.fillna('', inplace=True)

    seq_dict = {}
    for sample in report['samples']:
        seq_dict[sample['id']]= ''
    sample_list = seq_dict.keys()
    for gene_id, row in gene_df.iterrows():
        # Only run if it is core gene
        if len(row[row == 0]) != 0:
            continue
        gene_id = re.sub(r'\W+', '', gene_id)
        gene_dir = os.path.join(alignment_dir, gene_id)

        nucleotide_aln_file = os.path.join(gene_dir, gene_id + '.fna.aln')
        if not os.path.isfile(nucleotide_aln_file):
            logger.info('{} does not exist'.format(nucleotide_aln_file))
            continue

        for seq_record in SeqIO.parse(nucleotide_aln_file, 'fasta'):
            sample_name = re.findall(r'^(.+)_', seq_record.id)
            sample_name = sample_name[0]
            if sample_name not in sample_list:
                raise Exception(f'Error concatenating gene alignment: {sample_name} is not a sample id')
            seq_dict[sample_name] += str(seq_record.seq)

    with gzip.open(core_gene_aln_file, 'wt') as fh:
        for sample in sample_list:
            new_record = SeqRecord(Seq(seq_dict[sample]), id = sample, description = '')
            SeqIO.write(new_record, fh, 'fasta')

    return report