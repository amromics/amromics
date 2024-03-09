# AMRomics   


## Welcome to AMRomics

**AMRomics** is a software package to analyze microbial genomics data collections.
It provides a pipeline that bundles the current best practice for 
multiple aspects of AMR analyses. The pipeline analysis results can be 
represented and visualized via a web application. The web application also 
provides efficient data management.
 
AMRomics is written in python, it includes the followings dependencies:
 * blast (known to work with 2.10.1+)
 * samtools (1.11)
 * trimmomatic (0.39)
 * spades (3.15.2)
 * shovill (1.1.0)
 * flye (2.9)
 * prokka (1.14.6)
 * mlst (2.19.6)
 * abricate (1.0.1 | Database: vfdb ecoli_vf ecoh card megares resfinder argannot ncbi plasmidfinder)
 * roary (3.13.0) 
 * iqtree (2.1.2)
 * fasttree

## Installation
### Conda
The simplest method is installed via conda:

0. Make sure a conda version is installed in the computer. If not done already, download and install the appropriate conda, such as anaconda from 
   https://repo.anaconda.com/archive/
   
1. Create a conda environment with all the necessary dependencies: From the repository directory run
```bash
git clone --recursive https://github.com/amromics/amromics.git
cd amromics
conda create -y -c conda-forge -c defaults --name amromics python=3.10 mamba
source activate amromics
mamba install -y -c conda-forge -c bioconda -c anaconda -c etetoolkit -c rpetit3 -c defaults --file requirements.txt
pip install panta
```

2. Setup MLST database: AMRomics requires a copy of pubMLST database set up on the folder that AMRomics pipeline is run from. We make available the database in the accompanied file `db.tar.gz` updated on Feb 1, 2024. Just unzip the tarball 
```bash
tar zxvf db.tar.gz
```
Alternatively, you can update the latest database by running the following command line. Note the running make take some time depending on the bandwidth network and the responsiveness of pubMLST server
```
./amr-analysis.py download_db
```
### Docker
The pipeline can be installed via Docker as well.
```bash
git clone --recursive https://github.com/amromics/amromics.git
cd amromics
docker built -t amromics .
```
The working directory from the container is `/tmp/amromics`, user can run amromics commands by mounting the host working directory (where the git cloned into, e.g. `~/workspace/amromics`) into this destination (by using `-v`). For example, if user want to update the latest database: 
```bash
chmod 777 ~/workspace/amromics
docker run -v ~/workspace/amromics/:/tmp/amromics/ amromics amr-analysis.py download_db
```
Note that from the container working directory, a `db` is already available by unzipping the file `db.tar.gz` before hand.
## Usage

### Input preparation

AMRomics takes in as input a collection of bacterial samples in various data types: sequencing reads (fastq format), assembly (fasta format) or assembly with annotations (gff3 where annotations are followed by assembly in fasta). The list of samples in the collection is summarized in a tsv file with the following column headers:
- `sample_id`: An unique ID for each sample. sample_id has to be one word i.e., no space and contains only alphanumerical, _, and - characters
- `sample_name`: A description of the sample
- `input_type`: can take in one of `gff` (for annotations in gff3 format, expecting the assembly at the end of the gff file as in the input from prokka), `asm`, `assembly` (for assembly in fasta format), `Illumina` (for Illumina sequencing reads in fastq format), `pacbio-raw`, `pacbio-hifi`, `pacbio-corr` (for Pacbio sequencing reads in fastq format), `nano-raw`, `nano-hq`, `nano-corr` (for Nanopore sequencing reads in fastq format).
- `files`: path to the input data file(s). For Illumina input type, two fastq files can be supplied if paired-end sequencing is used, and they are separated by a semicolon. For all other input type, only one input file is expected.
- `trim`: if the value is `True` or `yes`, the read data will be trimmed - only applied for Illumina sequencing reads
- `genus`: the Genus of the sample, eg `Escherichia`
- `species`: the species of the sample, eg. `coli`
- `gsize`: an estimate of the genome size, in number, eg 5000000, needed only for subsample Illumina sequencing data to 100x
- `metadata`: any relevant information for the sample, in the format `key1:value1; key2:value2`. For example: `Geographic Location:Houston,USA;Insert Date:8/8/2017;Host Name:Human, Homo sapiens;ampicillin:Resistant;aztreonam:Resistant;ciprofloxacin:Resistant;gentamicin:Susceptible;tetracycline:Susceptible`

We provide several examples in the folder `examples` - See below.


AMRomics pipeline can be invoked with one command line `./amr-analysis.py`. Its usage is as follows:

```bash
./amr-analysis.py  pg -h
usage: amromics pg [-h] [-t THREADS] [-m MEMORY] -c COLLECTION_ID [-n COLLECTION_NAME] -i INPUT [--work-dir WORK_DIR] [--time-log TIME_LOG]
                   [--method METHOD] [--genetree GENETREE] [--progressive PROGRESSIVE] [--tree TREE] [--overwrite OVERWRITE] [--initdb {True,False}]

Pan-genome analysis of a collection

options:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use, 0 for all (default: 0)
  -m MEMORY, --memory MEMORY
                        Amount of memory in Gb to use (default: 30)
  -c COLLECTION_ID, --collection-id COLLECTION_ID
                        Collection ID (default: None)
  -n COLLECTION_NAME, --collection-name COLLECTION_NAME
                        Collection name (default: )
  -i INPUT, --input INPUT
                        Input file (default: None)
  --work-dir WORK_DIR   Working directory (default: data/work)
  --time-log TIME_LOG   Time log file (default: None)
  --method METHOD       panta or roary (default: panta)
  --genetree GENETREE   Run phylogenty for each gene cluster or not (default: False)
  --progressive PROGRESSIVE
                        Run pangenome in progressive mode (default: False)
  --tree TREE           fasttree or iqtree (default: fasttree)
  --overwrite OVERWRITE
                        Force overwrite exist results (default: False)
  --initdb {True,False}
                        Init full database (default: False)

```

#### Examples

We prepare several collections from public data. To download the raw data for these collections

To download the miniature dataset that consists of 5 samples: 1 in sequence assembly, 2 with Illumina sequencing reads, 1 with Nanopore and 1 with Pacbio sequencing reads.
```bash
cd examples/Kp24/raw
./download_kp4.sh
cd ../../
```

To download the 24 sample collection
```bash
cd examples/Kp24/raw
./download_kp24.sh
cd ../../
```

To download the 90 sample collection
```bash
cd examples/Kp24/raw
./download_kp24.sh
cd ../../
```

The following command will run that 24 samples through the pipeline, and import the results
to the web-app for visualization:

```bash
./amr-analysis.py pg --time-log k24_time.log  -t 7 -m 25 -c KpClinicalGRBZ -i examples/Kp24/Kp24.tsv --work-dir data/work  -n "Collection of 24 clinical isolates from Greek and Brazil"
```
or if using docker image `amromics` built as above:
```bash
docker run -v ~/workspace/amromics/examples/Kp24:/tmp/amromics/examples/Kp24 amromics amr-analysis.py pg --time-log examples/Kp24/k24_time.log  -t 16 -m 25 -c KpClinicalGRBZ -i examples/Kp24/config_Kp24.tsv --work-dir examples/Kp24/data/work  -n "Collection of 24 clinical isolates from Greek and Brazil"
```
#### Run with progressive mode:
```bash
./amr-analysis.py pg --time-log k24_time.log  -t 7 -m 25 -c KpClinicalGRBZ --progressive True -i examples/Kp89/Kp89.tsv --work-dir data/work  -n "Collection of 24+89 clinical isolates from Greek and Brazil"
```
### Output
Output from the pipeline is generated under the directory specified by `--work-dir`. 
The results include 2 output sub-folders corresponding to 2 stages of the pipeline: `samples/` for isolate analysis of each individual sample and `collections/` for pan-genome analysis results of the whole collection.
```
work-dir/
├── samples/
│   ├── sample1/
│   ├── sample2/
│   ├── ...
├── collections/<col_name>/
│   ├── alignments/
│   ├── pangenome/
│   ├── phylogeny/
│   ├── VCFs/
│   ├── sample_set.json
```
#### Individual samples
Results for each sample are written in its dedicated sub-folder, e.g. `sample1/` as below.
It stores final output from various modules, such as assembly, MLST, annotation, resistome/virolome detection using different databases:
```
sample1/
├── sample1_assembly.fasta
├── sample1_dump.json
├── sample1.faa
├── sample1.ffn
├── sample1.gff
├── sample1_mlst.tsv
├── sample1_plasmid.tsv
├── sample1_resistome_argannot.tsv
├── sample1_resistome_card.tsv
├── sample1_resistome_ecoh.tsv
├── sample1_resistome_megares.tsv
├── sample1_resistome_ncbi.tsv
├── sample1_resistome_resfinder.tsv
├── sample1_resistome.tsv
├── sample1_virulome.tsv
```
#### Collection pangenome
Analysis for the collection from AMRomics returns output in sub-folders below:
- `pangenome/`: pangenome results from `panta` or `roary` (gene clusters, representative sequences, gene presence/absence matrix...)
- `alignments/`: alignments of each gene cluster from the generated pangenome
- `phylogeny/`: core gene alignment and the corresponding phylogenetic tree
- `VCFs/`: variant calling for each sample from the pangenome (DNA and protein).
A pangenome reference, `pangenome_reference.fasta`, is built by concatenation of representative sequences for each and every gene clusters. This reference is used to generate variant profile for each sample, which is the concatenation of the variations of all its genes, in a VCF file.


## License

## Citations

## Authors
- Duc Quang Le <quangld@huce.edu.vn>
- Minh Duc Cao <minhduc.cao@gmail.com>
- Son Hoang Nguyen <dr.sonhoangnguyen@gmail.com>

