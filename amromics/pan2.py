import os
import shutil
import re
import json
import gzip
import csv
from datetime import datetime
import logging
from Bio import SeqIO
import pandas as pd
from amromics.utils import run_command
from amromics.pan_genome import extract_proteins, combine_proteins

logger = logging.getLogger(__name__)







def run_add_samples(report, collection_dir='.', threads=8, overwrite=False, timing_log=None):
    pan_genome_folder = os.path.join(collection_dir, 'pan_genome')
    temp_dir = os.path.join(collection_dir, 'temp')
    report['pan_genome'] = pan_genome_folder
    report['temp_dir'] = temp_dir
    if not os.path.exists(pan_genome_folder):
        os.makedirs(pan_genome_folder)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    report = extract_proteins(report, timing_log=timing_log)
    report = combine_proteins(report, timing_log=timing_log)


if __name__ == '__main__':
    report = json.load(open('dev/report2.json', 'r'))
    report = run_add_samples(
        report, 
        collection_dir='dev/pan2', 
        threads=4, 
        overwrite=True
    )
    json.dump(report, open('dev/pan2/report2_output.json', 'w'), indent=4, sort_keys=True)