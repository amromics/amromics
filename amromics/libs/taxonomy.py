
import os
import shutil
import csv
import logging
import multiprocessing
from amromics.utils.command import run_command
import amromics.libs.mlst as mlst
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def species_identification_kraken(prefix_name,assembly, db='db/kraken2/k2std', base_dir='.', timing_log=None,threads=0):
    """
        Run kraken2 for identifying species
        :param read_data: result holder
        :param base_dir: working directory
        :return: path to output file in result holder
    """
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir,   prefix_name+'_kraken2')
    kraken2_report = os.path.join(path_out, prefix_name + '_kraken2.tsv')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    cmd = 'kraken2 --db {db} --use-names --threads {threads} --report {report} {asm}'.format(
        db=db,
        threads=threads,
        report=kraken2_report,
        asm=assembly

    )



    cmd = "bash -c '{}'".format(cmd)
    if run_command(cmd, timing_log) != 0:
        return None
    return kraken2_report
###Sequence typing using mlst


def detect_mlst(sample):
    m = mlst.find_mlst(sample['assembly'])
    with open(sample['mlst'], 'w') as f:
        f.write("%s\t%s\t%s"%(m['file'],m['scheme'],m['st']))
        for gene in m['profile']:
            f.write("\t%s"%gene)
        f.write("\n")
        
