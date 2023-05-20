import os
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
def trim_pe_trimmomatic(prefix_name, reads,threads=0, base_dir='.', timing_log=None, **kargs):
    """
    read_data is a dictionary with field `sample_id`
    :param read_data:
    :param threads:
    :param base_dir:
    :param kargs:
    :return:
    """
    if threads <= 0:
        threads = NUM_CORES_DEFAULT

    out_dir = os.path.join(base_dir, 'trimmomatic_pe')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    out_p1 = os.path.join(out_dir, prefix_name + '_R1.fastq.gz')
    out_p2 = os.path.join(out_dir, prefix_name + '_R2.fastq.gz')

    out_s1 = os.path.join(out_dir, prefix_name + '_S1.fastq')
    out_s2 = os.path.join(out_dir, prefix_name + '_S2.fastq')

    cmd = 'trimmomatic.sh PE -threads {threads}'.format(threads=threads)
    cmd += ' {in_p1} {in_p2} {out_p1} {out_s1} {out_p2} {out_s2}'.format(
        in_p1=reads['pe1'], in_p2=reads['pe2'],
        out_p1=out_p1, out_s1=out_s1, out_p2=out_p2, out_s2=out_s2)
    ret = run_command(cmd, timing_log)
    # Combine single-ended reads into one
    out_s = os.path.join(out_dir, prefix_name+ '_S.fastq')
    with open(out_s, 'w') as fn:
        for seq in bioseq.read_sequence_file(out_s1):
            fn.write(seq.format_fastq())
        for seq in bioseq.read_sequence_file(out_s2):
            fn.write(seq.format_fastq())

    if ret == 0:

#        read_data['se'] = out_s
        return out_p1,out_p2

###NGS assembly using SPAdes
def assemble_spades(prefix_name,reads, base_dir = '.', threads=0, memory=50,overwrite=False, timing_log=None,trim=False,gsize=None, **kargs):
    if threads == 0:
        threads = NUM_CORES_DEFAULT


    path_out = os.path.join(base_dir, prefix_name + '_spades')
    assembly_file = os.path.join(path_out, prefix_name + '_contigs.fasta.gz')
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    elif os.path.isfile(assembly_file) and (not overwrite):
        return assembly_file
    #TODO: add trimmomatic if needed, add seqtk step
    if reads['pe1'].endswith('fastq') or reads['pe1'].endswith('fastq.gz') or reads['pe1'].endswith('fasta') or reads['pe1'].endswith('fasta.gz'):
        pass
    else:
        if get_compress_type(reads['pe1'])=='gzip':
            os.rename(reads['pe1'], reads['pe1']+'.fastq.gz')
            reads['pe1']= reads['pe1']+'.fastq.gz'
        else:
            os.rename(reads['pe1'], reads['pe1']+'.fastq')
            reads['pe1']= reads['pe1']+'.fastq'
    if reads['pe2'].endswith('fastq') or reads['pe2'].endswith('fastq.gz') or reads['pe2'].endswith('fasta') or reads['pe2'].endswith('fasta.gz'):
        pass
    else:
        if get_compress_type(reads['pe2'])=='gzip':
            os.rename(reads['pe2'], reads['pe2']+'.fastq.gz')
            reads['pe2']= reads['pe2']+'.fastq.gz'
        else:
            os.rename(reads['pe2'], reads['pe2']+'.fastq')
            reads['pe2']= reads['pe2']+'.fastq'
    cmd = 'spades.py -m {memory} -t {threads} -k 77,99,127 --careful -o {path_out}'.format(
        memory=int(memory), threads=threads, path_out=path_out)
    if 'pe1' in reads and 'pe2' in reads:
        cmd += ' --pe1-1 {pe1} --pe1-2 {pe2}'.format(pe1=reads['pe1'], pe2=reads['pe2'])
    if 'se' in reads:
        cmd += ' --s1 {se}'.format(se=reads['se'])

    ret = run_command(cmd, timing_log)
    if ret != 0:
        return None

    #Read in list of contigs
    contigs = list(bioseq.read_sequence_file(os.path.join(path_out, 'contigs.fasta')))
    #TODO: filter based on coverage

    contigs = sorted(contigs, key=len, reverse=True)
    logger.info("Read in {} contigs".format(len(contigs)))
    #assembly_file = os.path.join(path_out, prefix_name + '_contigs.fasta')
    with gzip.open(assembly_file, 'wt') as f:
        for i, contig in enumerate(contigs):
            contig.set_desc(contig.get_name())
            contig.set_name(prefix_name + '_C' + str(i+1))
            f.write(contig.format_fasta())
    return assembly_file
def assemble_skesa(prefix_name, reads,base_dir = '.', threads=0, memory=50, timing_log=None, **kargs):
    if threads == 0:
        threads = NUM_CORES_DEFAULT
    path_out = os.path.join(base_dir, prefix_name + '_skesa')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    cmd = 'skesa --memory {memory} --cores {threads} --fastq '.format(
        memory=int(memory), threads=threads)
    if 'pe1' in reads and 'pe2' in reads:
        cmd += '{pe1} {pe2}'.format(pe1=read_data['pe1'], pe2=read_data['pe2'])
    if 'se' in reads:
        cmd += '{se}'.format(se=reads['se'])
    assembly_file_raw= os.path.join(path_out, 'contigs.fasta')
    cmd+=' >{output}'.format(output=assembly_file_raw)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        return None

    #Read in list of contigs
    contigs = list(bioseq.read_sequence_file(assembly_file_raw))
    #TODO: filter based on coverage
    assembly_file = os.path.join(path_out, prefix_name + '_contigs.fasta')
    contigs = sorted(contigs, key=len, reverse=True)
    logger.info("Read in {} contigs".format(len(contigs)))

    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            contig.set_desc(contig.get_name())
            contig.set_name(prefix_name + '_C' + str(i+1))
            f.write(contig.format_fasta())
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
