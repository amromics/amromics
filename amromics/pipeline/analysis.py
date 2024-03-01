import json
import os
import shutil

from amromics.pipeline import wrapper
def copy_file(source_file, dest_dir):
    #try:
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    dest_file = os.path.join(dest_dir, os.path.basename(source_file))
    shutil.copyfile(source_file, dest_file)
    return dest_file
def prepareDataCollectionAnalysis(report,collection_dir):
    temp_folder = os.path.join(collection_dir, 'temp_data')
    #clean up the folder
    if not os.path.isdir(temp_folder):
        os.makedirs(temp_folder)

    gff_dir=os.path.join(temp_folder,'gff')
    #collect ffn files: annotation files should be named as <sample_id>.ffn
    ffn_dir=os.path.join(temp_folder,'ffn')
    faa_dir=os.path.join(temp_folder,'faa')

    for sample in report['samples']:
        #print(sample)
        copy_file(sample['annotation_gff'],gff_dir)
        copy_file(sample['annotation_faa'],faa_dir)
        copy_file(sample['annotation_ffn'],ffn_dir)

    return temp_folder,gff_dir,ffn_dir,faa_dir

def single_genome_analysis(samples, work_dir, assembly_method='spades', overwrite=False, threads=0, memory=8, timing_log=None):
    # TODO: move the analysis in single_genome_analysis_func here
    for sample in samples:
        sample_id = sample['id']
        sample_dir = os.path.join(work_dir, 'samples', sample_id)
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
            
        wrapper.run_single_sample(
            sample, sample_dir=sample_dir, assembly_method=assembly_method, threads=threads, overwrite=overwrite,
            memory=memory,trim=sample['trim'], timing_log=timing_log)

        with open(os.path.join(sample_dir, sample_id + '_dump.json'), 'w') as fn:
            json.dump(sample, fn)
    return samples


def pan_genome_analysis(samples, work_dir, collection_id, collection_name=None, progressive=False,overwrite=False, threads=0, memory=8, timing_log=None,method='panta',rate_coverage=0.15,  genetree=True, tree='fasttree'):

    report = {'collection_id': collection_id,
              'collection_name': collection_name,
              'samples': samples}
    collection_dir = os.path.join(work_dir, 'collections', collection_id)

    # Check if the pan-genone analysis need to be updated
    dataset_sample_ids = []
    for sample in report['samples']:
        sample_id = sample['id']
        dataset_sample_ids.append(sample_id)
        overwrite = overwrite or sample['updated']

    if not os.path.exists(collection_dir):
        os.makedirs(collection_dir)

    # to check if the set of samples has not changed
    sample_set_file = os.path.join(collection_dir, 'sample_set.json')
    if os.path.isfile(sample_set_file):
        with open(sample_set_file) as fn:
            sample_set = json.load(fn)

        if set(sample_set) != set(dataset_sample_ids):
            overwrite = True

    if overwrite:
        roary_folder = os.path.join(collection_dir, 'roary')
        if os.path.exists(roary_folder):
            shutil.rmtree(roary_folder)

        phylogeny_folder = os.path.join(collection_dir, 'phylogeny')
        if os.path.exists(phylogeny_folder):
            shutil.rmtree(phylogeny_folder)

        # Note: Check for existing alignment is done within
    # Write the set of sample IDs
    with open(sample_set_file, 'w') as fn:
        json.dump(dataset_sample_ids, fn)
    #report,genome_dir,gff_dir,ffn_dir,reference, base_dir='.', threads=0, memory=50
    temp_folder,gff_dir,ffn_dir,faa_dir=prepareDataCollectionAnalysis(report,collection_dir)

    report = wrapper.run_collection(report,gff_dir,ffn_dir,faa_dir,overwrite=overwrite,base_dir=collection_dir,progressive=progressive, threads=threads,timing_log=timing_log,method=method, rate_coverage=rate_coverage, genetree=genetree, tree=tree)
    with open(os.path.join(collection_dir, collection_id + '_dump.json'), 'w') as fn:
        json.dump(report, fn)

    # # TODO clean up tmp files
    # if os.path.exists(collection_dir + "/temp"):
    #     shutil.rmtree(collection_dir + "/temp")
    # extract_json.export_json(work_dir, webapp_data_dir, collection_id, collection_name)
    return report
