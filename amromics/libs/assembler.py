import os

import logging
import multiprocessing
from Bio import SeqIO
import amromics.libs.bioseq as bioseq
from amromics.utils.command import run_command
from amromics.utils.utils import get_open_func
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
logger = logging.getLogger(__name__)

###NGS assembly using SPAdes
def assemble_spades(sample_id,reads, assembly_file, base_dir = '.', threads=4, memory=50, timing_log=None, **kargs):
    
    path_out = os.path.join(base_dir, sample_id + '_spades')    
    if not os.path.exists(path_out):
        os.makedirs(path_out)    
    
    #cmd = 'spades.py -m {memory} -t {threads} -k 77,99,127 --isolate --disable-gzip-output -o {path_out}'.format(
    cmd = 'spades.py -m {memory} -t {threads} --only-assembler --isolate --cov-cutoff auto -o {path_out}'.format(
        memory=int(memory), threads=threads, path_out=path_out)
    if 'pe1' in reads and 'pe2' in reads:
        cmd += ' -1 {pe1} -2 {pe2}'.format(pe1=reads['pe1'], pe2=reads['pe2'])
    if 'se' in reads:
        cmd += ' -s {se}'.format(se=reads['se'])
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception(f'Fail to assemble sample {sample_id} with spades  ({ret})')

    
    # Read in list of contigs
    contigs = list(SeqIO.parse(os.path.join(path_out, 'contigs.fasta'), "fasta"))
    contigs = sorted(contigs, key=len, reverse=True)
    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            contig.id =sample_id+'_C'+str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")    
    
    #remove intermediate sub-folder
    run_command('rm -rf ' + os.path.join(path_out,'corrected'))
    run_command('rm -rf ' + os.path.join(path_out,'misc'))
    run_command('rm -rf ' + os.path.join(path_out,'tmp'))

    if 'pe1' in reads and 'pe2' in reads:
        os.remove(reads['pe1'])
        os.remove(reads['pe2'])
    if 'se' in reads:
        os.remove(reads['se'])
    #clean up



def assemble_skesa(sample_id, reads, assembly_file, base_dir = '.', threads=4, memory=50, timing_log=None, **kargs):
    
    path_out = os.path.join(base_dir, sample_id + '_skesa')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

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
        raise Exception(f'Fail to assemble sample {sample_id} with SKESA ({ret})')

    # Read in list of contigs
    contigs = list(SeqIO.parse(os.path.join(path_out, 'contigs.fasta'), "fasta"))
    contigs = sorted(contigs, key=len, reverse=True)
    logger.info("Read in {} contigs".format(len(contigs)))    
    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            contig.id =sample_id+'_C'+str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")
    
    if 'pe1' in reads and 'pe2' in reads:
        os.remove(reads['pe1'])
        os.remove(reads['pe2'])
    if 'se' in reads:
        os.remove(reads['se'])

    #return assembly_file
            

def assemble_shovill(prefix_name, reads,base_dir, trim=False, threads=4, memory=50, overwrite=False, timing_log=None,gsize=None):
    """
        Run assembly process for pair-end input using shovill
        :param prefix_name: sample name, sample id, etc
        :param read: reads
        :param base_dir: working directory
        :return: path to assembly file (normailized)
    """

    path_out = os.path.join(base_dir, prefix_name + '_shovill')
    assembly_file = os.path.join(path_out, prefix_name + '_contigs.fasta')

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
    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            contig.id =prefix_name+'_C'+str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")
    run_command('rm -f ' + os.path.join(path_out, 'spades.fasta'))
    run_command('rm -f ' + os.path.join(path_out, 'contigs.fa'))
    run_command('rm -f ' + os.path.join(path_out, 'contigs.gfa'))
    run_command('gzip ' + os.path.join(path_out, 'shovill.log'))
    return assembly_file

def get_assembly(sample_id, assembly_input, assembly_file):
    """
    Get the assembly from the input assembly to the target assembly

    Parameters
    ----------
    input_assembly:
        
    assembly: str
        the directory of the sample
    Returns
    -------
        None
    """
    
    logger.info(f'Get assembly {assembly_file} from {assembly_input} ')
    open_func = get_open_func(assembly_input)
    with open_func(assembly_input, 'rt') as fn:
        contigs = list(SeqIO.parse(fn, 'fasta'))    
    
    contigs = sorted(contigs, key=len, reverse=True)
    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            contig.id = sample_id + '_C' + str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")    
    
def get_assembly_from_gff(sample_id, gff, assembly_file, base_dir='.', overwrite=False):
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

    path_out = os.path.join(base_dir, sample_id + '_assembly')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    open_func = get_open_func(gff)
    
    #TODO: no write to temp file
    temp_assembly = os.path.join(path_out,sample_id+'_temp.fasta')
    with open_func(gff, 'rt') as fn, open(temp_assembly,'w') as t_out:
        start_fasta=False
        for line in fn:
            if start_fasta:
                t_out.write(line)
            if line.startswith('##FASTA'):
                #Done reading gff, move on to reading fasta
                start_fasta=True
    
    with open(temp_assembly) as fn:
        contigs = list(SeqIO.parse(fn, 'fasta'))
    
    contigs = sorted(contigs, key=len, reverse=True)
    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            #contig.id = prefix_name + '_C' + str(i)
            #contig.description = ''
            SeqIO.write(contig, f, "fasta")
    os.remove(temp_assembly)
    

def assemble_flye(sample_id, reads, input_type, assembly_file, base_dir='.', threads=4, overwrite=False, timing_log=None, gsize=None):
    """
        Run assembly process for long reads input using flye
        :param sample_id: sample name, sample id, etc
        :param read: reads
        :param input_type: type of long reads
        :assembly_file
        :param base_dir: working directory
        
    """

    path_out = os.path.join(base_dir, sample_id + '_flye')
    #assembly_file = os.path.join(path_out, sample_id + '_contigs.fasta')

    #if not os.path.exists(path_out):
    #    os.makedirs(path_out)
    #elif os.path.isfile(assembly_file) and (not overwrite):
    #    logger.info(f'Assembly file {assembly_file} found, skip assembling')
    #    return assembly_file

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
    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            contig.id =sample_id+'_C'+str(i)
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
    os.remove(reads['long-read'])
    