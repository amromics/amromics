import os
import shutil
import csv
import logging
import multiprocessing
import gzip
from amromics.utils.command import run_command
from amromics.libs.bioseq import read_sequence_file
from amromics.utils.utils import get_compress_type
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()

def trim_trimmomatic(prefix_name, reads,threads=0, base_dir='.', overwrite=False, timing_log=None, **kargs):
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

    out_dir = os.path.join(base_dir, 'trimmomatic')
    out_p1 = os.path.join(out_dir, prefix_name + '_R1.fastq.gz')
    out_p2 = os.path.join(out_dir, prefix_name + '_R2.fastq.gz')
    out_s = os.path.join(out_dir, prefix_name + '_S.fastq.gz')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    elif not overwrite: 
        if os.path.isfile(out_p1) and os.path.isfile(out_p2):
            return out_p1,out_p2
        elif os.path.isfile(out_s):
            return out_s

    TRIMOPT='ILLUMINACLIP:db/trimmomatic.fa:1:30:11 LEADING:3 TRAILING:3 MINLEN:30 TOPHRED33'
    if 'pe1' in reads and 'pe2' in reads:
        cmd = 'trimmomatic PE -threads {threads} -phred33'.format(threads=threads)
        cmd += ' {in_p1} {in_p2} {out_p1} /dev/null {out_p2} /dev/null {opt}'.format(
            in_p1=reads['pe1'], in_p2=reads['pe2'],
            out_p1=out_p1, out_p2=out_p2, opt=TRIMOPT)
        ret = run_command(cmd, timing_log)

        if ret == 0:
            return out_p1,out_p2

    elif 'se' in reads:
        cmd = 'trimmomatic SE -threads {threads} -phred33'.format(threads=threads)
        cmd += ' {in_s} /dev/null {opt}'.format(
            in_s=reads['se'], opt=TRIMOPT)
        ret = run_command(cmd, timing_log)

        if ret == 0:
            return out_s

    else:
        raise Exception('ERROR: Trimmomatic only for Illumina PE or SE!')


def trim_fastp(prefix_name, reads,threads=0, base_dir='.', overwrite=False, timing_log=None, **kargs):
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

    out_dir = os.path.join(base_dir, 'fastp')
    out_p1 = os.path.join(out_dir, prefix_name + '_R1.fastq.gz')
    out_p2 = os.path.join(out_dir, prefix_name + '_R2.fastq.gz')
    out_s = os.path.join(out_dir, prefix_name + '_S.fastq.gz')
    report_json=os.path.join(out_dir,prefix_name+".json")
    report_html=os.path.join(out_dir,prefix_name+".html")

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    elif not overwrite: 
        if os.path.isfile(out_p1) and os.path.isfile(out_p2):
            return out_p1,out_p2
        elif os.path.isfile(out_s):
            return out_s
    #fastp -w 8 --in1 SRR1616936_1.fastq.gz --in2 SRR1616936_2.fastq.gz --out1 fastp/SRR1616936_1.fastq.gz --out2 fastp/SRR1616936_2.fastq.gz -j fastp/report.json -h fastp/report.html
    if 'pe1' in reads and 'pe2' in reads:
        
        cmd = 'fastp -w {threads} --in1 {in_p1} --in2 {in_p2} --out1 {out_p1} --out2 {out_p2} -j {report_json} -h {report_html}'.format(
            threads=threads, 
            in_p1=reads['pe1'], in_p2=reads['pe2'],
            out_p1=out_p1, out_p2=out_p2, 
            report_json=report_json, report_html=report_html)
        
        ret = run_command(cmd, timing_log)

        if ret == 0:
            #return out_p1,out_p2,report_json,report_html
            return out_p1,out_p2

    elif 'se' in reads:
        cmd = 'fastp -w {threads} --in1 {in_s} --out1 {out_s} -j {report_json} -h {report_html}'.format(
            threads=threads, 
            in_s=reads['se'],
            out_s=out_s,
            report_json=report_json, report_html=report_html)
        ret = run_command(cmd, timing_log)

        if ret == 0:
            #return out_s,None,report_json,report_html
            return out_s

    else:
        raise Exception('ERROR: Fastp only for Illumina PE or SE!')

def estimate_gsize_mash(prefix_name, reads, overwrite=False, threads=0, base_dir='.', timing_log=None,**kargs):
    ##estimate genome size
    #gsize=$(mash sketch -p 8 -o /dev/null -k 21 -m 5 -r fastp/SRR1616936_1.fastq.gz 2>&1 | awk -F: '/Estimated genome size/{printf("%d",$2)}')
    
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    mash_log = os.path.join(base_dir, prefix_name + '_mash.log')
    
    ret=0
    if os.path.isfile(mash_log) and (not overwrite):
        print("mash results exist for {}! skip running mash sketch...".format(prefix_name))
    else:
        mash_out = os.path.join(base_dir, prefix_name + '_mash.out')
        if 'pe1' in reads and 'pe2' in reads:
            in1 = reads['pe1']
        elif 'se' in reads:
            in1 = reads['se']
        else:
            raise Exception('ERROR: Mash (gsize est.) only for Illumina PE or SE!')

        cmd='mash sketch -p {threads} -o {mash_out} -k 21 -m 5 -r {in1} > {mash_log} 2>&1'.format(
            threads=threads,
            mash_out=mash_out, in1=in1, mash_log=mash_log)
        cmd+=' && rm -f {mash_out}'.format(mash_out=mash_out)
        ret=run_command(cmd,timing_log)

    if ret != 0:
        print("ERROR: seqtk for {} return non-zero code! skip subsampling by seqtk...".format(prefix_name))
    else:
        with open (mash_log, 'r') as f:
            for num, line in enumerate(f, 1):
                if 'Estimated genome size' in line:
                    return int(float(line.split(':')[1]))
    return None

def rename_reads(reads):
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


def subsample_seqtk(prefix_name, reads, expect_depth=100, threads=8, base_dir='.', overwrite=False, timing_log=None, gsize=None,**kargs):

    ##orig_depth=total_basas/gsize
    #expect_depth=100 #100X
    #factor=$(awk -v gsize=${gsize} -v edepth=${expect_depth} -F: '$1~/after_filtering/{flag=1}$1~/total_bases/&&flag{printf("%.2f",edepth/($2/gsize)); exit}' fastp/fastp.json)

    #if (( $(echo "$factor < 1" |bc -l) )) ; then
    #	    seqtk sample fastp/SRR1616936_1.fastq.gz $factor | pigz --fast -c -p 8 > subsample/SRR1616936_1.fastq.gz
    #	    seqtk sample fastp/SRR1616936_2.fastq.gz $factor | pigz --fast -c -p 8 > subsample/SRR1616936_2.fastq.gz
    #else
    #	    ln -s $(realpath fastp/SRR1616936_1.fastq.gz) subsample
    #	    ln -s $(realpath fastp/SRR1616936_2.fastq.gz) subsample
    #fi
    if gsize==None:
        return reads
    totBases=0
    if 'pe1' in reads:
        for seq in read_sequence_file(reads['pe1']):
            totBases+=seq.length()
    if 'pe2' in reads:
        for seq in read_sequence_file(reads['pe2']):
            totBases+=seq.length()
    if 'se' in reads:
        for seq in read_sequence_file(reads['se']):
            totBases+=seq.length()
    if totBases==0:
        logger.info("Raw data of {} is neither PE or SE! skip subsampling by seqtk...".format(prefix_name))
        return reads

    depth=totBases/gsize
    if depth < expect_depth:
        logger.info(f"{prefix_name} has less than 100x! skip subsampling by seqtk... ({gsize} {totBases} {depth})")
        return reads
    else:
        factor=expect_depth/depth
        path_out=os.path.join(base_dir, prefix_name + '_subsample')
        out_p1 = os.path.join(path_out, prefix_name + '_R1.fastq.gz')
        out_p2 = os.path.join(path_out, prefix_name + '_R2.fastq.gz')
        out_s = os.path.join(path_out, prefix_name + '_S.fastq.gz')
        if not os.path.exists(path_out):
            os.makedirs(path_out)
        elif ((os.path.isfile(out_p1) and os.path.isfile(out_p2)) or os.path.isfile(out_s)) and (not overwrite):
            logger.info("{} subsampling has been done! skip subsampling by seqtk...".format(prefix_name))
            return reads
        
        if 'pe1' in reads:
            cmd='seqtk sample {in_p1} {factor} | pigz --fast -c -p {threads} > {out_p1}'.format(
                    in_p1=reads['pe1'], factor=factor, threads=threads, out_p1=out_p1)
        if 'pe2' in reads:
            cmd+='&&seqtk sample {in_p2} {factor} | pigz --fast -c -p {threads} > {out_p2}'.format(
                    in_p2=reads['pe2'], factor=factor, threads=threads, out_p2=out_p2)
        if 'se' in reads:
            cmd='seqtk sample {in_s} {factor} | pigz --fast -c -p {threads} > {out_s}'.format(
                    in_s=reads['se'], factor=factor, threads=threads, out_s=out_s)
        logger.info("Running subsampling command: {}".format(cmd)) 
        ret=run_command(cmd,timing_log)

        if ret != 0:
            logger.info("ERROR: seqtk for {} return non-zero code! skip subsampling by seqtk...".format(prefix_name))
        elif 'pe1' in reads and 'pe2' in reads:
            reads['pe1']=out_p1
            reads['pe2']=out_p2
        elif 'se' in reads:
            reads['se']=out_p1

        return reads

    
