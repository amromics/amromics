import os,shutil
import logging
import multiprocessing
from Bio import SeqIO
import csv
import pandas as pd
from datetime import datetime

from amromics.utils.command import run_command
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def run_roary(samples, overwrite=False,threads=0, base_dir='.', timing_log=None):
    """
        Run roay make pangeome analysis (using prokka results in previous step)
        :param read_data: result holder
        :param base_dir: working directory
        :param threads: number of core CPU
        :return:
    """

    roary_folder=os.path.join(base_dir,'pangenome/roary')
    #temp_folder = os.path.join(base_dir, 'pangenome/temp_roary')
    roary_output = os.path.join(roary_folder, 'summary_statistics.txt')
    if os.path.isfile(roary_output) and (not overwrite):
        logger.info('roary has run and the input has not changed, skip roarying')
        return roary_folder
    
    gff_list = []    
    for sample in samples:
        gff_list.append(sample['annotation_gff'])
    
    # Make sure the directory is not there or roary will add timestamp
    if os.path.isfile(roary_folder):
        os.remove(roary_folder)
    if os.path.exists(roary_folder):
        shutil.rmtree(roary_folder)
    cmd = f'roary -p {threads} -f {roary_folder} -v ' + ' '.join(gff_list)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('roary fail to run!')    
    return roary_folder

def run_panta_cmd(samples, AL=0.8, rate_coverage=0.15, dont_split=False, overwrite=False, progressive=False, threads=0, base_dir='.', timing_log=None):
    starttime = datetime.now()
    sample_gffs = []
    for sample in samples:
        sample_gffs.append(sample['annotation_gff'])

    out_folder=os.path.join(base_dir,'pangenome/panta')
    outfile= os.path.join(out_folder, 'gene_presence_absence.csv')
    if os.path.isfile(outfile) and (not overwrite) and (not progressive):
        logger.info('panta has run and the input has not changed, skip run panta')
        return out_folder
    
    logger.info(f"progressive={progressive}  at {os.path.isfile(outfile)}")
    if os.path.isfile(outfile) and progressive:
        cmd = f'panta add -s -c {out_folder} -t {threads} --AL {AL} -a protein -r  {rate_coverage} -g '
        cmd += ' '.join(sample_gffs)
        logger.info('Run panta with progressive mode')
    else:        
        cmd = f'panta main -s -o {out_folder} -t {threads} --AL {AL} -a protein -r {rate_coverage} -g '
        cmd += ' '.join(sample_gffs)
        logger.info('Run panta with normal mode')
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('panta fail to run!')
    elapsed = datetime.now() - starttime
    logger.info(f'Done -- time taken {str(elapsed)}')
    return out_folder
