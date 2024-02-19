import os
from glob import glob
import shutil
import csv
import logging
import gzip
import multiprocessing
from Bio import SeqIO
import amromics.libs.bioseq as bioseq
from amromics.utils.command import run_command
from amromics.utils.utils import get_open_func
from amromics.utils.utils import get_compress_type
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
logger = logging.getLogger(__name__)

###NGS assembly using SPAdes
def assemble_spades(prefix_name,reads, base_dir = '.', threads=0, memory=50,overwrite=False, timing_log=None,trim=False,gsize=None, **kargs):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, prefix_name + '_spades')
    assembly_file = os.path.join(path_out, prefix_name + '_contigs.fasta.gz')
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    elif os.path.isfile(assembly_file) and (not overwrite):
        logger.info(f'Assembly file {assembly_file} found, skip assembling')
        return assembly_file
   
    #rename reads zip/unzip now move to preprocess.py
    #cmd = 'spades.py -m {memory} -t {threads} -k 77,99,127 --isolate --disable-gzip-output -o {path_out}'.format(
    cmd = 'spades.py -m {memory} -t {threads} --only-assembler --isolate --cov-cutoff auto -o {path_out}'.format(
        memory=int(memory), threads=threads, path_out=path_out)
    if 'pe1' in reads and 'pe2' in reads:
        cmd += ' -1 {pe1} -2 {pe2}'.format(pe1=reads['pe1'], pe2=reads['pe2'])
    if 'se' in reads:
        cmd += ' -s {se}'.format(se=reads['se'])

    ret = run_command(cmd, timing_log)
    if ret != 0:
        return None

    #remove intermediate sub-folder
    #run_command('rm -rf ' + os.path.join(path_out,'corrected'))
    run_command('rm -rf ' + os.path.join(path_out,'misc'))
    run_command('rm -rf ' + os.path.join(path_out,'tmp'))    

    # Read in list of contigs
    contigs = list(SeqIO.parse(os.path.join(path_out, 'contigs.fasta'), "fasta"))
    contigs = sorted(contigs, key=len, reverse=True)
    logger.info("Read in {} contigs".format(len(contigs)))
    #assembly_file = os.path.join(path_out, prefix_name + '_contigs.fasta')
    with gzip.open(assembly_file, 'wt') as f:
        for i, contig in enumerate(contigs):
            contig.id =prefix_name+'_C'+str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")
    return assembly_file


def assemble_skesa(prefix_name, reads,base_dir = '.', threads=0, memory=50, overwrite=False, timing_log=None, **kargs):
    if threads == 0:
        threads = NUM_CORES_DEFAULT
    path_out = os.path.join(base_dir, prefix_name + '_skesa')
    
    assembly_file = os.path.join(path_out, prefix_name + '_contigs.fasta.gz')
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    elif os.path.isfile(assembly_file) and (not overwrite):
        logger.info(f'Assembly file {assembly_file} found, skip assembling')
        return assembly_file

    cmd = 'skesa --memory {memory} --cores {threads} --fastq '.format(
        memory=int(memory), threads=threads)
    if 'pe1' in reads and 'pe2' in reads:
        cmd += '{pe1} {pe2}'.format(pe1=reads['pe1'], pe2=reads['pe2'])
    if 'se' in reads:
        cmd += '{se}'.format(se=reads['se'])
    assembly_file_raw= os.path.join(path_out, 'contigs.fasta')
    cmd+=' >{output}'.format(output=assembly_file_raw)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        return None

    # Read in list of contigs
    contigs = list(SeqIO.parse(os.path.join(path_out, 'contigs.fasta'), "fasta"))
    contigs = sorted(contigs, key=len, reverse=True)
    logger.info("Read in {} contigs".format(len(contigs)))
    assembly_file = os.path.join(path_out, prefix_name + '_contigs.fasta.gz')
    with gzip.open(assembly_file, 'wt') as f:
        for i, contig in enumerate(contigs):
            contig.id =prefix_name+'_C'+str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")
    return assembly_file

def assemble_shovill(prefix_name, reads,base_dir, trim=False, threads=4, memory=50, overwrite=False, timing_log=None,gsize=None):
    """
        Run assembly process for pair-end input using shovill
        :param prefix_name: sample name, sample id, etc
        :param read: reads
        :param base_dir: working directory
        :return: path to assembly file (normailized)
    """

    path_out = os.path.join(base_dir, prefix_name + '_shovill')
    assembly_file = os.path.join(path_out, prefix_name + '_contigs.fasta.gz')

    if not os.path.exists(path_out):
        os.makedirs(path_out)
    elif os.path.isfile(assembly_file) and (not overwrite):
        return assembly_file

    cmd = 'shovill  --ram {memory} --cpus {threads} --outdir {path_out}'.format(
        memory=int(memory), threads=threads, path_out=path_out)
    if trim:
        cmd += ' --trim --depth 200'
    else:
        cmd += ' --depth 120'
    if 'pe1' in reads and 'pe2' in reads:
        cmd += ' --force  --R1 {pe1} --R2 {pe2}'.format(pe1=reads['pe1'], pe2=reads['pe2'])
    else:
       raise Exception('Only support pair-end reads!')
    if gsize:
        cmd +=' --gsize {gsize}'.format(gsize=gsize)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running shovill!')

    # Read in list of contigs
    contigs = list(SeqIO.parse(os.path.join(path_out, 'contigs.fa'), "fasta"))
    contigs = sorted(contigs, key=len, reverse=True)
    with gzip.open(assembly_file, 'wt') as f:
        for i, contig in enumerate(contigs):
            contig.id =prefix_name+'_C'+str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")
    run_command('rm -f ' + os.path.join(path_out, 'spades.fasta'))
    run_command('rm -f ' + os.path.join(path_out, 'contigs.fa'))
    run_command('rm -f ' + os.path.join(path_out, 'contigs.gfa'))
    run_command('gzip ' + os.path.join(path_out, 'shovill.log'))
    return assembly_file
def get_assembly(prefix_name,assembly, base_dir, overwrite=False):
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

    path_out = os.path.join(base_dir, prefix_name + '_assembly')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    open_func = get_open_func(assembly)
    with open_func(assembly, 'rt') as fn:
        contigs = list(SeqIO.parse(fn, 'fasta'))

    assembly_file = os.path.join(path_out,prefix_name + '_contigs.fasta.gz')
    if os.path.isfile(assembly_file) and (not overwrite):
        logger.info(f'Not overwrite {assembly_file}')
        return assembly_file

    contigs = sorted(contigs, key=len, reverse=True)
    with gzip.open(assembly_file, 'wt') as f:
        for i, contig in enumerate(contigs):
            contig.id = prefix_name + '_C' + str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")


    return assembly_file
def get_assembly_from_gff(prefix_name,gff, base_dir, overwrite=False):
    """
    Get the assembly from user input

    Parameters
    ----------
    gff:
        gff file with assembly
    sample_dir: str
        the directory of the sample
    Returns
    -------
        path to assembly file
    """

    path_out = os.path.join(base_dir, prefix_name + '_assembly')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    open_func = get_open_func(gff)
    temp_asm=open(os.path.join(path_out,prefix_name+'_temp.fasta'),'w')
    with open_func(gff, 'rt') as fn:
        start_fasta=False
        for line in fn:
            if start_fasta:
                temp_asm.write(line)
            if line.startswith('##FASTA'):
                #Done reading gff, move on to reading fasta
                start_fasta=True

    temp_asm.close()
    with open(os.path.join(path_out,prefix_name+'_temp.fasta'), 'rt') as fn:
        contigs = list(SeqIO.parse(fn, 'fasta'))

    assembly_file = os.path.join(path_out,prefix_name + '_contigs.fasta.gz')
    if os.path.isfile(assembly_file) and (not overwrite):
        logger.info(f'Not overwrite {assembly_file}')
        return assembly_file

    contigs = sorted(contigs, key=len, reverse=True)
    with gzip.open(assembly_file, 'wt') as f:
        for i, contig in enumerate(contigs):
            #contig.id = prefix_name + '_C' + str(i)
            #contig.description = ''
            SeqIO.write(contig, f, "fasta")


    return assembly_file
def assemble_flye(prefix_name, reads, input_type, base_dir, threads=4, overwrite=False, timing_log=None, gsize=None):
    """
        Run assembly process for long reads input using flye
        :param prefix_name: sample name, sample id, etc
        :param read: reads
        :param input_type: type of long reads
        :param base_dir: working directory
        :return: path to assembly file (normailized)
    """

    path_out = os.path.join(base_dir, prefix_name + '_flye')
    assembly_file = os.path.join(path_out, prefix_name + '_contigs.fasta.gz')

    if not os.path.exists(path_out):
        os.makedirs(path_out)
    elif os.path.isfile(assembly_file) and (not overwrite):
        logger.info(f'Assembly file {assembly_file} found, skip assembling')
        return assembly_file

    cmd = 'flye --threads {threads} --out-dir {path_out} --{input_type} {reads}'.format(
        threads=threads, path_out=path_out, input_type=input_type, reads=reads['long-read'])
    if gsize:
        cmd +=' -g {gsize}'.format(gsize=gsize)

    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running flye!')

    # Read in list of contigs
    contigs = list(SeqIO.parse(os.path.join(path_out, 'assembly.fasta'), "fasta"))
    contigs = sorted(contigs, key=len, reverse=True)
    with gzip.open(assembly_file, 'wt') as f:
        for i, contig in enumerate(contigs):
            contig.id =prefix_name+'_C'+str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")

    # clean
    run_command('rm -rf ' + os.path.join(path_out, '00-assembly'))
    run_command('rm -rf ' + os.path.join(path_out, '10-consensus'))
    run_command('rm -rf ' + os.path.join(path_out, '20-repeat'))
    run_command('rm -rf ' + os.path.join(path_out, '30-contigger'))
    run_command('rm -rf ' + os.path.join(path_out, '40-polishing'))
    run_command('rm -f ' + os.path.join(path_out, 'params.json'))
    run_command('rm -f ' + os.path.join(path_out, 'assembly.fasta'))
    run_command('gzip ' + os.path.join(path_out, 'assembly_graph.gfa'))
    run_command('gzip ' + os.path.join(path_out, 'assembly_graph.gv'))
    run_command('gzip ' + os.path.join(path_out, 'assembly_info.txt'))
    run_command('gzip ' + os.path.join(path_out, 'flye.log'))

    return assembly_file
