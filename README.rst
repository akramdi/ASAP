===================================
README for ASAP (v0.1)
===================================
ASAP version 0.1

Introduction 
============

ASAP is a flexible bioinformatic pipeline for ATAC-seq data analysis. Starting from raw ATAC-seq sequencing reads, ASAP outputs raw and filtered mapping files, coverage files (reads coverage ; tn5 insertion events coverage), fragment length distribution, read exraction based on fragment length, and peak calling results. 

Overview of major steps 
==========================

- Mapping 
- Post-mapping processing and filtering:

  - Filter (or not) reads that fall into user-defined blacklisted regions
  - Select reads that do not carry more than minMismatch
  - Mark duplicated pairs
  - Select concordant, non-duplicated pairs. 
  - Shift read pairs as described in Buenorestro et al.,2013
- Compute read coverage
- Compute insertion events coverage
- Fragment length distribution
- Extract reads pairs based on a fragment length range
- Peak calling

ASAP is:
 **User-friendly**: requires a single configuration file. Thus, only one option is required when running the command line (see `Usage of ASAP`_)
 

 **Flexible**: provides the possibility to skip a given step(s) and target specific post-processing step(s).


Dependencies
============

* `Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_   version >= 2.2.9 
* `MACS2 <https://github.com/taoliu/MACS>`_ 
* `Samtools <http://samtools.sourceforge.net/>`_ version >=1.2
*  `GenomeCoverageBed <http://bedtools.readthedocs.io/en/latest/index.html>`_  version >= v2.20.1
* `igvtools <https://software.broadinstitute.org/software/igv/igvtools>`_  IGV Version >= 2.3.72


Usage of ASAP
=============
A configuration file required to execute the pipeline. 

::
 
 bash ASAP.sh [-h] [-v] [-c]



Options
--------

-c CONFIGFILE
```````````````
This is the only REQUIRED parameter for ASAP. The configuration file is a text file that gathers the full set of parameters required to execute the pipeline. (check the example *ASAP_configFile_example.conf* in distribution)

-h/-v 
``````
Print out the help/current version


About the configuration file:
=============================

The configuration file gathers the parameters of each step. Note that, when running the pipeline, only the "turned on" steps will performed. A step is turned on by a *yes/no* argument.

Here we list the different set of parameters to be filled in the configuration file: 








General parameters
------------------
General information option about the run. Must be always filled. 



:OUTDIR:              Main output directory where results are written. OUTDIR is created if does not exist
:sampleName:          Name of the processed sample. No space is allowed: use _ or - to mimic space if needed
:CHRLEN:              Chromosome info file (tab-delimited format: *<Chr name><chr length>*)
:path:                Full path to the different dependencies, if not already added to $PATH


Mapping step parameters
-----------------------
:map:                         Set to "*yes/no*". If mapping is skipped (*map=no*), a BAM file must be provided to proceed. 
                              (see post-mapping steps).
:FASTQ1:                      fastq file (R1). File can be gizpped
:FASTQ2:                      fastq file (R2). File can be gizpped
:bowtieIndex:                 Prefix of bowtie2 indexes
:mappingParameters: Bowtie2  mapping parameters. Default: *--very-sensitive -X 2000 -p 10*

 
Post mapping steps 
-------------------
It is possible to skip the mapping step (*map=no*) and perform any of the post-mapping steps. To do so, aligned reads must be provided in a BAM file. If* map=yes*, the "turned on" post-mapping steps will performed on the internal mapping results.

:BAM: aligment file in BAM format


Filtering parameters
---------------------

:filter:                     Set to "*yes/no*". If *map=yes*, filtering will be performed on internal mapping results, 
                             if *map=no*, filtering will be performed on the provided alignment file in BAM option. 
                             
:maxMis:                     Maximum number of mismatches allowed per read.
:blacklist:                  Set to "yes/no" if reads should be filtred based on a list of blacklisted regions. 
                             If "*blacklist=yes*", blacklisted regions must be provided in the next parameter. 

:blacklistedRegions:         Regions used to filter reads.(tab-delimited format: <Chr name><start><end>)

:shift:                      Set to "yes"/"no". If *shift=yes*, reads are shifted by 4bp so that read starts reflect the center of the Tn5 transposition event

Coverage
---------
:readCoverage:                Set to "*yes/no*" if read coverage should be computed or not
:ieventsCoverage:             Set to "*yes/no*" if Tn5 insertion events coverage should be computed or not
:GENOME:                      Genome assembly. Use same genomes names as igvtools (tair10, hg38, mm9..)

Fragment length
---------------
:fragDist:                    Set to "*yes/no*" if fragment length distribution should be computed or not


Read extraction
---------------
:extractReads:                Set to "*yes/no*" if read pairs should be extracted based on a given range of fragment length
:lowBoundary:                 Lower boundery of the range: [lowBoundary,upBoundary]. Default=100
:upBoundary:				  Upper boundery of the range: [lowBoundary,upBoundary]. Default=250          


Peak calling
------------
:callpeak:                     Set to "yes/no" if peak calling should be computed or not.
:control:                      Control bam file. Note that peak calling can be performed without a control, however, one can                            provide a control such as ATAC-seq on genomic DNA. Leave option empty if no control is used.
:MODE:                         Peak calling mode: *<broad/narrow>*. Default=broad
:fdr:                          Cutoff for peak detection. Default=0.01
:gsize:                        Effective genome size of tair10 (gsize=10e7)



Output files
============

ASAP outputs mapping files, coverage files, fragments distribution table/plot and MACS2 peak calling results.
Mapping output
---------------

:*.mapped.sorted.bam:                Contains mapped reads (bowtie2 raw mapping results)

Filtering/post-processing outputs
---------------------------------

:*.masked.shifted.bam: Contains the selected set of reads after filtering. Ideally, accessible peaks are called using this file. 

:*.filter.stats.csv: Summary of filtering step is CSV format

Coverage outputs
----------------
:*.masked.shifted.tdf: Genome-wide coverage of ATAC reads 
:*.masked.shifted.ievent.tdf: Genome-wide coverage of Tn5 insertion events

:*.masked.shifted.ievent.bam: Contains Tn5 insertion events. Basically, instead of showing reads, only the position corresponding to Tn5 insertion event are shown)


Fragment length distribution
----------------------------
:TLEN.{sampleName}.f66.txt: Counts/frequencies of fragments length
:TLEN.{sampleName}.f66.txt: Plot of fragment length distribution

Read extraction
---------------
:*.subReads.f3.frag*.bam: Contains the set of extracted reads based on the given rage of fragment length


Peak calling outputs 
--------------------
Output are stored in an directory: *peak_calling_<sampleName>*. Check `MACS2 output list <https://github.com/taoliu/MACS#output-files>`_


