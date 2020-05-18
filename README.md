# Motif Enrichment Pipeline - Lixing Yang Lab

## Introduction:

This pipeline is intended to automate the motif enrichment analysis around the breakpoints of structural variants (SVs). The program first extracts SV breakpoints from the input bedpe files: for each SV breakpoint, the program first checks whether it falls within the mappable regions of the chromosome, and then generates random positions within that same chromosome's mappable regions. The pipeline then extracts the 100 bp length sequence from either side of the breakpoint from the reference genome for each real and random SV breakpoint. The program then runs both the AME and FIMO algorithm to determine whether a specific sequence motif is enriched in the real SVs breakpoint sites versus the randomly generated ones.  
  
The Python script should be run as follows:  
**`$ python runner_analysis.py <options>`**  

## Requirements:

Python 3: https://www.python.org/downloads/release/python-373/  
Bedtools: https://bedtools.readthedocs.io/en/latest/  
Pandas: https://pandas.pydata.org/  
Matplotlib: https://matplotlib.org/  
MEME Suite: http://meme-suite.org/  
Parsl: https://parsl.readthedocs.io/en/stable/index.html  
Scipy: https://www.scipy.org/  

*On Gardner, the modules to load are: gcc/6.2.0 perl/5.24.0 zlib/1.2.8 meme/5.0.5 bedtools/2.27.1 R/3.6.1 miniconda3/4.7.10*  
*The user must then activate a previously created conda virtual environment with Python 3, Parsl, Pandas, Scipy and Matplotlib installed, e.g.*  
`$ module load gcc/6.2.0 perl/5.24.0 zlib/1.2.8 meme/5.0.5 bedtools/2.27.1 R/3.6.1 miniconda3/4.7.10`  
`$ source activate motifenrichment_env`  
`(motifenrichment_env)$ python runner_analysis.py <options>`  

## Options:  

**`-h` or `--help` :** displays the list of options and exits the program  
**`-i` or `--input_dir <path/to/directory>` :** specify the directory path containing all the input bedpe files (REQUIRED)  
**`-o` or `--output_dir <path/to/directory>` :** specify the directory for all output files (REQUIRED)  
**`-f` or `--genome_fasta <path/to/file>` :** specify the genome fasta file path (REQUIRED)  
**`-l` or `--genome_len <path/to/file>` :** specify the chromosome length file path (REQUIRED)  
**`-e` or `--genome_include <path/to/file>` :** specify the file path of the genome regions to include (REQUIRED)  
**`-m` or `--motif_path <path/to/file>` :** specify the path of the motif file in MEME format (REQUIRED)  
**`-t` or `--sampleinfo_table <path/to/file>` :** specify the sample attribute table file path *[default: all samples in the input directory will be selected]*  
**`-a` or `--sample_attr <attribute column:specific attribute,attribute column:specific attribute>` :** specify the sample attribute(s) to select for in the sample info table *[default: all samples in the input directory will be selected]*  
**`-s` or `--SV_types <tra,inv,del,dup>` :** specify the SV types to select for *[default: tra,inv,del,dup]*  
**`-r` or `--rand_sv_ratio <int>` :** specify the ratio of random SVs to real SVs *[default: 3:1 random to real SVs, e.g. 3]*  
**`-F` or `--FIMO_thresh <float>` :** specify the p-value threshold for the FIMO algorithm *[default: 0.0001]*  
**`-A` or `AME_scoring <avg|max>` :** specify the scoring method for AME *[default: max]*  
**`c` or `--config <filename>` :** specify the config file to use for Parsl parallel processing *[default: local]*  
**_NOTE: The paths for all files and directories in the command line should always be absolute paths._**  

## Input:

* Input directory: should contain all necessary bedpe files, with names formated as such `<samplename>.*.bedpe`. Chromosome names should be without "chr", e.g. 1 or X.
* Output directory: if the output directory already contains a `/bed_files` directory, only bed files missing from this directory will be generated while running the program.
* Genome fasta file: again, the fasta headers should be chromosome names without "chr".
* Genome length file: should be a tsv file with columns "chrom" and "length". Again, chromosome names should not start with "chr".
* Genome include file: should be a bed file. Here, the chromosome names start with "chr", e.g. chr1 or chrX. Please see k100 single read mappability files from Umap project (https://bismap.hoffmanlab.org/).
* Motif file: should be in MEME format (see: http://meme-suite.org/doc/meme-format.html?man_type=web)
* Sample info table: should be a csv file. The first column will be the list of sample names, and the remaining columns should be attributes that can be selected for with `--sample_attr`. The first row should be the column names ("attribute_column") and the remaining rows the attributes corresponding to each sample ("specific_attribute").
* Config file: this file should be located in the `/configs` directory of this program. See `gardner.py` or `local.py` for an example, and the Parsl documentation for additional information. When specifying the config in the command line arguments, ignore the `.py` file extension (i.e. `local` not `local.py`).

## Ouput:

* `bed_files/` directory: contains all the randomly generated or real breakpoint bed files listing all the 100bp sequences to analyze (with names formated as such `<sample_name>_<sv_type>_sv.bed` for real SVs, and `<samplename>_<SV type>_rand<random to real ratio>.bed`).
* Results directory (name format: `<attribute column>-<specific attribute>+<attribute column>-<specific attribute>_<SV type>+<SV type>_<random to real ratio>_AME-<AME scoring method>_FIMO-<FIMO threshold>`): contains the following files:
    * `*_results_summary.txt` file: contains all summary of all the statistical analyses.
    * `AME/` directory: contains the raw output of the AME algorithm.
    * `FIMO_(rand|sv)/` directory: contains the raw output of the FIMO algorithm for the real or random SV breakpoints.
    * `*_(rand|sv).fasta` and `*_(rand|sv).bed` files: bed and fasta files containing all breakpoint locations or sequences respectively for random and real SVs.
    * `*_AME_results.csv` file: contains the scoring of all SV breakpoints during AME analysis (used to create the plot).
    * `*_plot_(rand|sv).pdf` file: contains the histogram of the relative rank of AME scores for real or random SV breakpoints.
