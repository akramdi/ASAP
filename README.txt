# Introduction
===============

ASAP: ATAC-Seq data Analysis Pipeline.

ASAP is a flexible bioinformatic pipeline that performs a full ATAC-seq data analysis; starting from raw sequencing reads to peak calling. Outputs include raw and filtered mapping files, coverage files (reads coverage and/or tn5 insertion coverage	), plus peak calling results. 

Fulls steps:

I.Mapping
II.Post-mapping processing (mark duplicate pairs) and filtering
 II.1. Select reads that fall into predefined selected regions
 II.2. Select reads that do not carry more than minMismatch
 II.3. Mark duplicated pairs
 II.4. Select concordant, non-duplicated pairs. 
 II.5. Shift read pairs as described in Buenorestro et al.,2013
III. Compute coverage
 III.1. Compute read coverage
 III.2. Compute insertion events coverage
IV. Peak calling

ASAP is:

*User-friendly: requires a single configuration file. Thus, only one option is required when running the command line (-c configFile.conf). 
*Flexible: provides the possibility to skip a given step(s) and target specific post-processing step(s).

For question/suggestions: kramdi@biologie.ens.fr

# Prerequisites
===============

bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): >= version 2.2.9 
MACS2 (https://github.com/taoliu/MACS): >= versiobn 2.1.0.20151222
samtools (http://samtools.sourceforge.net/): >= version: 1.2
genomeCoverageBed (bedtools: http://bedtools.readthedocs.io/en/latest/index.html): >= v2.20.1
igvtools (https://software.broadinstitute.org/software/igv/igvtools):  IGV Version >= 2.3.72

# Run example: 
==============

Parameters must be provided in a form of a configuration file which is the only input accepted in order to run the pipeline. A configuration file example can be found in: /kingdoms/a2e/programs/a2e-ak-scripts/NGSpipelines/ATAC/ASAP_configFile_example.conf

Command example:
bash /kingdoms/a2e/programs/a2e-ak-scripts/NGSpipelines/ATAC/ASAP.sh -c myConfig_file.conf

IMPORANT: on Bioclust, the script must be called from its directory.


# About the configuration file:
==============================
The configuration file gathers the parameters of every step. Note that, when running the pipeline, only the "turned on" steps will performed. 


Content of the configuration file:

--> General options section: general information option about the run. Must be always filled. 
OUTDIR: Main output directory where results are written. OUTDIR is created if does not exist.
sampleName: Name of the processed sample. No space is allowed: use _/- to mimic space if needed.
path* : indicate the full path to the different dependencies. Default values are set for the Bioclust environment.


--> Mapping options: 
map: set to "yes/no".If mapping is skipped (map=no), a BAM file must be provided to proceed. (see post-mapping steps).
FASTQ1: absolute path to fastq file (R1). File can be gizpped.
FASTQ2: absolute path to fastq file (R2). File can be gizpped.
bowtieIndex: absolute path to prefix of bowtie2 indexes
mappingParameters: bowtie2 mapping parameters

--> post mapping steps: if one is interested only in post-mapping steps, it is possible to skip the mapping step (map=no) and perform filtering and/or peak calling and/or coverage computation. To do so, aligned reads must be provided in a BAM file. If map=yes, the "turned on" post-mapping steps will performed on the internal mapping results.

BAM: absolute path to a BAM file

------> Filtering step:
filter: set to "yes"/"no". If map=yes, filtering will be performed internal mapping results, if map=no, filtering will be performed on the provided bam file in BAM option. 
maxMis: Maximum number of mismatches allowed per read.
selectedRegions: regions used to select reads.(tab-delimited format: <Chr name><start><end>)

------> Coverage step:
readCoverage: set to "yes/no" if read coverage should be computed or not
ieventsCoverage: set to "yes/no" if Tn5 insertion event coverage should be computed or not
GENOME: genome assembly (tair10)
CHRLEN: absolute path to chromosome info (tab-delimited format: <Chr name><chr length>)


------> Peak calling step: 
callpeak: set to "yes/no" if peak calling should be computed or not.
control: absolute path to a control bam file. Note that peak calling can be performed without a control, however, one can provide a control such as ATAC on genomic DNA. Leave option empty if no control is used.

MODE: peak calling mode: <broad/narrow>. Default=broad
fdr: cutoff for peak detection. Default=0.01
gsize: effective genome size of tair10 (gsize=10e7)



# How to fill the configuration file:
====================================

- Do not change the name of the option.
- No space is allowed between "=" and option value
- BAM option is taken into account only when mapping step is skipped (map=no)


# Output files
==============
ASAP outputs mapping files, coverage files and MACS2 peak calling results.

----> Mapping results:
*mapped.sorted.bam: Contains mapped reads (bowtie2 raw mapping results)

----> Filtering/post-processing results:
*mis.mkdup.f3F1024.masked.shifted.bam: Contains the selected set of reads after filtering. Ideally, accessible peaks are called using this output. 
*.csv: Summery of filtering step is CSV format

----> Coverage results
*mis.mkdup.f3F1024.masked.shifted.tdf:genome-wide coverage of ATAC reads 
*mis.mkdup.f3F1024.masked.shifted.ievent.tdf: genome-wide coverage of Tn5 insertion events

(Extra output: *mis.mkdup.f3F1024.masked.shifted.ievent.bam: contains Tn5 insertion events. Basically, instead of showing reads, only the position corresponding to Tn5 insertion event are shown)

----> Peak calling output: 
Output are stored in an directory: peak_calling_<sampleName>. List of MACS2 output and description: https://github.com/taoliu/MACS#output-files


# Memory requirement (approximately, tested on mapping of ~ 200M reads)
====================
*Disk memory: >= 25G
*Whole pipeline: 10cpu, ~7G memory
*Only mapping step: 10 cpu, 3G memory
*Mapping step skipped: 1cpu, ~7G memory


# Contact 
==========
Questions, suggestions: kramdi@biologie.ens.fr



