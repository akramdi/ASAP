#!/bin/bash 


:<<'readme'
15/05/2017
kramdi@biologie.ens.fr

This script is written in an effort to gather the main analysis steps of ATAC-seq data.Here are the major steps:

I.Mapping
II.Post-mapping processing (mark duplicate pairs) and filtering
 II.1. Select reads that fall into predefined selected regions
 II.2. Select reads that do not carry more than minMismatch
 II.3. Select concordant, non-duplicated pairs. 
 II.4. Shift read pairs as described in Buenorestro et al.,2013
III. Compute coverage
 III.1. Compute read coverage
 III.2. Compute insertion events coverage
IV. Compute fragment length distribution
 V. Peak calling

readme

#Here, we get the directory of this script in order to display a proper help
#SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#===============================================================================
# Help setion
#===============================================================================
usage() { 
cat << EOF
ASAP: ATAC-Seq Analysis Pipeline

Flexible processing of ATAC-seq data from mapping; filtering, to coverage and peak calling.

USAGE: bash $0 -c myConfig_file.conf
       
OPTIONS:
	-c <Parameters file>: Full set of parameters. Format and example file in : $SCRIPTDIR 
       

A.KRAMDI - a2eTeam - May 2017
		 
EOF

}

 

while getopts "vhc:" OPTION
do
     case $OPTION in
		 v)
             echo "ASAP : beta version 0.9.1"
             exit 1
             ;;

         h)
             usage
             exit 1
             ;;
        c)
            CONFIG=$OPTARG
			if [[ ! -s $CONFIG ]]; then >&2 echo "[ERROR]: Empty/non-existing parameter file or parameter file not correctly provided (use -c option).Exit."; exit 1; fi
			STEP="all"
             ;;
	       ?)
          usage
          exit
          ;;
     esac
done

if [  $# -eq 0 ]; then usage ; exit 1 ; fi 

#===============================================================================
# Create timestamp for log
#===============================================================================
date=`date "+[%b-%d-%Y %T]"`
stamp="@$date"


#===============================================================================
# Create tmp dir and LOGS
#===============================================================================

LOCALTMP=/localtmp/tmp$RANDOM
if [ ! -d $LOCALTMP ]; then mkdir $LOCALTMP ; fi 

function finish {
rm -rf $LOCALTMP
}
trap finish EXIT 
#LOG file is created after getting samplename/projcetName

#===============================================================================
# Parse configuration file to extract parameters 
#===============================================================================

#Initialize argument array
typeset -A config # init array
config=( # set default values in config array
[OUTDIR]=""
[sampleName]=""
[FASTQ1]=""
[FASTQ2]=""
[bowtieIndex]=""
[mappingParameters]=""
[maxMis]=""
[gsize]=""
[fdr]=""
[MODE]=""
[BAM]=""
[map]=""
[filter]=""
[readCoverage]=""
[ieventsCoverage]=""
[callpeak]=""
[fragDist]=""
[control]=""
[CHRLEN]=""
[GENOME]=""
[pathSamtools]=""
[pathIgvTools]=""
[pathGenomeCoverageBed]=""
[pathMACS2]=""
[pathBowtie2]=""
[pathToJava]=""
[pathToPicardJar]=""
)


while read line
do
    if echo $line | grep -F = &>/dev/null
    then
        varname=$(echo "$line" | cut -d '=' -f 1)
        config[$varname]=$(echo "$line" | cut -d '=' -f 2-)
    fi
done < $CONFIG


#Create variables
OUTDIR=${config[OUTDIR]}
ID=${config[sampleName]} #changed variable name
FASTQ1=${config[FASTQ1]}
FASTQ2=${config[FASTQ2]}
bowtieIndex=${config[bowtieIndex]}
mappingParameters=${config[mappingParameters]}
maxMis=${config[maxMis]}
SELECTEDREGIONS=${config[selectedRegions]} #changed variable name
gsize=${config[gsize]}
fdr=${config[fdr]}
MODE=${config[MODE]}
BAM=${config[BAM]}
map=${config[map]}
filter=${config[filter]}
CHRLEN=${config[CHRLEN]}
GENOME=${config[GENOME]}
readCoverage=${config[readCoverage]}
ieventsCoverage=${config[ieventsCoverage]}
callpeak=${config[callpeak]}
fragDist=${config[fragDist]}
peakCallingControl=${config[control]} #changed variable name
pathSamtools=${config[pathSamtools]}
pathIgvTools=${config[pathIgvTools]}
pathGenomeCoverageBed=${config[pathGenomeCoverageBed]}
pathMACS2=${config[pathMACS2]}
pathBowtie2=${config[pathBowtie2]}
pathToJava=${config[pathToJava]}
pathToPicardJar=${config[pathToPicardJar]}

#================================== Check some parameters and issue warnings/errors

#check yes/no parameters
if [[ ( "$map" != "yes" && "$map" != "no" ) || ( "$filter" != "yes" && "$filter" != "no" ) || ( "$readCoverage" != "yes" && "$readCoverage" != "no" )  || ( "$ieventsCoverage" != "yes" && "$ieventsCoverage" != "no" ) || ( "$fragDist" != "yes" && "$fragDist" != "no" ) || ( "$callpeak" != "yes" && "$callpeak" != "no" ) ]] ; then 
	>&2 echo -e "\n[ERROR]: Possible values to turn on a section: yes/no. Exit." ; exit 1 
fi

if [ -d $OUTDIR ]; then >&2 echo "[WARNING]: Output directory ($OUTDIR) already exists. Results will be overwitten." ;fi

#check if BAM option is filled when map=no. If BAM is filled, check if a valid BAM file was given 
if [[ "$map" == "no" && -z $BAM ]]; then 
	>&2 echo -e "\n[ERROR]: Mapping step is set to 'no', yet no valid BAM file is provided. Please provide a path to an input BAM file. Exit." ; exit 1 
fi
if [[ "$map" == "no" && ! -s $BAM ]]; then 
	>&2 echo -e "\n[ERROR]: Mapping step is set to 'no', yet no valid BAM file is provided. Please provide a path to an input BAM file. Exit." ; exit 1 
fi
# chech in input option is valid (if not filled the peak calling will be performed without a control)
if [[ ! -z $peakCallingControl && ! -s $peakCallingControl ]]; then 
	>&2 echo -e "\n[ERROR]: Error due to control file for peak calling: $peakCallingControl. File does not exist or has zero size. Please provide a correct control file or leave option empty. Exit." ; exit 1 
fi

# Create a coverage directory if any of the coverage steps are turned on: 
COVDIR=$OUTDIR/Coverage_${ID}
if [[ ( "$readCoverage" == "yes" || "$ieventsCoverage" == "yes" ) && ! -d $COVDIR ]]; then mkdir $COVDIR ;fi 

#Create fragment dist directory
FRAGDIR=$OUTDIR/Fragment_distribution_${ID}
if [[  "$fragDist" == "yes" && ! -d $FRAGDIR ]]; then mkdir $FRAGDIR ;fi 



#print out parameters and info
echo -e "\t=================================================================================\n\tUSER: $USER\n\tSample: $ID\n\tASAP start time: `date`\n\t================================================================================="

echo "++++++++++++++++++"
echo -e "PARAMETERS SUMMARY:"
echo "++++++++++++++++++"

#write general parameters
cat > $LOCALTMP/tmp.conf <<EOL
# Output directory: $OUTDIR
# Sample: $ID
# Mapping step: $map
EOL

#Write mapping section + options
if [[ "$map" == "yes" ]]; then 
cat >> $LOCALTMP/tmp.conf <<EOL
- Bowtie2: $pathBowtie2
- Index: $bowtieIndex
- Mapping parameters: $mappingParameters
EOL
fi

#write filtering section + options
cat >> $LOCALTMP/tmp.conf <<EOL
# Filtering step: $filter 
EOL
if [[ "$filter" == "yes" ]]; then 
cat >> $LOCALTMP/tmp.conf <<EOL
- Select reads that fall oustide blacklisted regions
- Maximum mismatches per read: $maxMis
- Select concordant, non-duplicated pairs
[ Final set of reads: concordant, non-duplicated read pairs that carry no more than $maxMis mismatch and falls outside of blacklisted regions ]
EOL
fi


#write coverage section/fragment distribution section/peak calling section + options
cat >> $LOCALTMP/tmp.conf <<EOL
# Compute read coverage: $readCoverage
# Compute insertion events coverage: $ieventsCoverage
# Compute fragment length distribution: $fragDist
# Peak calling: $callpeak
EOL


#Check if an control file was provided
if [ -z $peakCallingControl ]; then
	controlInfo="Control: no control"
else
	controlInfo="Control: $peakCallingControl"
fi

#peak calling parameters
if [[ "$callpeak" == "yes" ]]; then
cat >> $LOCALTMP/tmp.conf <<EOL
- MACS2: $pathMACS2
- $controlInfo
- Treatement: $BAM
- Mode: $MODE
- Cutoff for peak detection (fdr): $fdr
- $GENOME effective genome size: 10e7
- Other parameters: --keep-dup all --nomodel
EOL
fi

cat $LOCALTMP/tmp.conf

echo "++++++++++++++++++++++"
echo -e "END PARAMETERS SUMMARY"
echo -e "++++++++++++++++++++++\n"


#IMPORTANT: Create log file where all commands are kept. Comes in handy in case of debugging
LOGDIR=$OUTDIR/logs_${ID} ; if [ ! -d $LOGDIR ] ; then mkdir $LOGDIR ;fi 
LOG=$LOGDIR/commands_ASAP_run_${ID}.log
touch $LOG
echo -e "++++++++++++++++++++++++++++++ List of main commands and logs ++++++++++++++++++++++++++++++++++\nUSER: $USER\nSample: $ID\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n" > $LOG


########################################################################################### MAIN #################################################################################"
cd $OUTDIR


if [[ $map == "yes" ]]; then
#===============================================================================
# Mapping step
#===============================================================================
echo -e "\n######################## Mappinng step" | tee -a $LOG

echo "$stamp: Mapping..." | tee -a $LOG
echo "#$pathBowtie2 -x $bowtieIndex $mappingParameters -1 $FASTQ1 -2 $FASTQ2 1> ${ID}.mapped.sam" | tee -a $LOG
eval "$pathBowtie2 -x $bowtieIndex $mappingParameters -1 $FASTQ1 -2 $FASTQ2 1> ${ID}.mapped.sam "


echo "$stamp: Convert mapping outpur to bam..." | tee -a $LOG
echo "#$pathSamtools view -bhS ${ID}.mapped.sam -o ${ID}.mapped.bam" >> $LOG
$pathSamtools view -bhS ${ID}.mapped.sam -o ${ID}.mapped.bam
 

echo "$stamp: Sort bam..." | tee -a $LOG
echo "#$pathSamtools sort ${ID}.mapped.bam -o ${ID}.sorted.bam -O bam -T $LOCALTMP" >> $LOG
$pathSamtools sort  ${ID}.mapped.bam -o ${ID}.mapped.sorted.bam -O bam -T $LOCALTMP

if [ -s ${ID}.mapped.sorted.bam ]; then rm -f ${ID}.mapped.sam ${ID}.mapped.bam ; fi 

echo "$stamp: Index sorted bam..." | tee -a $LOG
echo "#$pathSamtools index ${ID}.mapped.sorted.bam" >> $LOG
$pathSamtools index ${ID}.mapped.sorted.bam

#if mapping step is skipped. Index the provided BAM file
else

	echo "Input BAM: $BAM"
	if [ ! -s $BAM.bai ]; then echo "$stamp: Index input BAM..." ; $pathSamtools index $BAM ; fi
fi




if [[ "$filter" == "yes" ]]; then 
if [[ "$map" == "yes" ]]; then BAM=${ID}.mapped.sorted.bam ;fi 
#===============================================================================
# Post-processing and filtering
#===============================================================================
echo -e "\n######################## Post-processing and filtering step (using BAM:$BAM)" | tee -a $LOG

#this file will hold the number of reads before and after the mismatch filter for informative purpose.(from the start,we only consider the reads that fell outside the blacklisted regions)
IDnreads=$LOCALTMP/$ID.nreads.txt

#echo "** Moved to $FINALBAMDIR" >> $LOG

#Select reads that fall into our selected regions + N mis + markduplicates + get concordant-non-duplicates 
echo "$stamp: Filter by blacklisted regions and allowed mismatches..." | tee -a $LOG
$pathSamtools view -h -F 4 $BAM -L $SELECTEDREGIONS | awk -v n=$maxMis -v ID=$ID -v IDnreads=$IDnreads  '{ a++ ; if($0 ~ /^@/) {print $0} else { mis=match($0,"XM:i:") ; nmis=substr($0,mis+5,2) ;if (int(nmis)<=n) {print $0; b++} } }END{printf("selected-regions\t%d\nreads-selected-by-nmis\t%d\n",a,b) > IDnreads }'  | $pathSamtools sort - -o $ID.${maxMis}mis.masked.sorted.bam -O bam -T $LOCALTMP



echo "$stamp: Mark duplicates..." | tee -a $LOG
echo "#$pathToJava -Xmx5g -Djava.io.tmpdir=$LOCALTMP \
-jar $pathToPicardJar MarkDuplicates \
METRICS_FILE=$LOGDIR/MKDUP.$ID.${maxMis}mis.masked.log  \
INPUT=$ID.${maxMis}mis.masked.sorted.bam \
OUTPUT=$ID.${maxMis}mis.mkdup.masked.bam 2> $log/MKDUP.$ID.${maxMis}mis.masked.err" >> $LOG


$pathToJava -Xmx5g -Djava.io.tmpdir=$LOCALTMP \
-jar $pathToPicardJar MarkDuplicates \
METRICS_FILE=$LOGDIR/MKDUP.$ID.${maxMis}mis.masked.log  \
INPUT=$ID.${maxMis}mis.masked.sorted.bam \
OUTPUT=$ID.${maxMis}mis.mkdup.masked.bam 2> $LOGDIR/MKDUP.$ID.${maxMis}mis.masked.err


#Get number of duplicates (from the reads falling into our selected regions and carrying no more than nmis) 
echo "$stamp: Get stats on duplicates..." 
nDup=`$pathSamtools view -f 1024 $ID.${maxMis}mis.mkdup.masked.bam -c`

echo "$stamp: Select concordant, non-duplicated reads..."   | tee -a $LOG
echo "#$pathSamtools view -hb -f 3 -F 1024 $ID.${maxMis}mis.mkdup.masked.bam >  $ID.${maxMis}mis.mkdup.f3F1024.masked.bam  " >> $LOG
 
$pathSamtools view -hb -f 3 -F 1024 $ID.${maxMis}mis.mkdup.masked.bam >  $ID.${maxMis}mis.mkdup.f3F1024.masked.bam 

#remove intermediate files
if [ -s  $ID.${maxMis}mis.mkdup.f3F1024.masked.bam ] ; then rm -f $ID.${maxMis}mis.mkdup.masked.bam $ID.${maxMis}mis.masked.sorted.bam ; fi

#========4. Fix mates!
echo "$stamp: Shift reads..." 
$pathSamtools sort -n $ID.${maxMis}mis.mkdup.f3F1024.masked.bam -T $LOCALTMP -O bam | $pathSamtools fixmate - $ID.tmp.bam

#========5. Shift reads! and sort by coordinates for indexing 
$pathSamtools view -h -F 16 $ID.tmp.bam |  awk 'BEGIN{OFS="\t"} { if ($1 ~ /^@/) {print} else { Rlen=length($10) ; $4=$4+4 ; $8=$8-5 ; if ($8<0 && $8==0) {$8=1} ;  print $0   }}' > $ID.tmp.shifted.sam
$pathSamtools view -f 16 $ID.tmp.bam | awk 'BEGIN{OFS="\t"} { $4=$4-5 ;if ($4<0 && $4==0) {$4=1} ;  $8=$8+4 ;  print $0 }' >> $ID.tmp.shifted.sam
$pathSamtools view -bh $ID.tmp.shifted.sam | $pathSamtools sort - -o $ID.${maxMis}mis.mkdup.f3F1024.masked.shifted.bam -O bam -T $LOCALTMP

# index file
echo "$stamp: Index output bam file..."
$pathSamtools index $ID.${maxMis}mis.mkdup.f3F1024.masked.shifted.bam

#remove tmp bam file
if [ -s $ID.${maxMis}mis.mkdup.f3F1024.masked.shifted.bam ]; then rm -f $ID.tmp.bam $ID.tmp.shifted.sam  $ID.${maxMis}mis.mkdup.f3F1024.masked.bam ; fi 

echo "$stamp: Extract filtering stats in CSV file..."

Rawstats=$LOCALTMP/samtools.stats.inBAM.$ID.txt ; touch $Rawstats ; echo "# RAW SAMTOOLS STATS" >> $LOG ; cat $Rawstats >> $LOG
$pathSamtools flagstat $BAM | awk 'BEGIN{OFS="\t"}{print "#flagstat", $0}' > $Rawstats
$pathSamtools idxstats $BAM | awk 'BEGIN{OFS="\t"}{print "#idx", $0}'  >> $Rawstats

totalFastqReads=`cat $Rawstats | grep "#flagstat"  | awk '{print $2}' | awk '{FS="+"} NR ==1 {print $1}'`
echo "$stamp: Total reads in $ID fastq: $totalFastqReads" 

#total unmapped/mapped (no mismatche filtering)
totalMa=`cat $Rawstats | grep "#flagstat" | awk '{print $2}' | awk '{FS="+"} NR ==5 {print $1}'`
totalUnma=$((totalFastqReads-totalMa)) 
echo "$stamp: Unmapped reads: $totalUnma"
echo "$stamp: Mapped reads: $totalMa"


#To get the number of reads falling into blacklisted regions, we subtract the number of reads that DID not fall into the black listed regions - total number of reads
#BLreads=`echo "$totalFastqReads-`grep selected $IDnreads | cut -f2`" | bc`

SLreads=`grep regions $IDnreads | cut -f2`
echo "$stamp: Reads falling into selected regions: $SLreads "

SLreadsMis=`grep nmis $IDnreads | cut -f2`
echo "$stamp: Reads carrying selected by $maxMis (within selected regions): $SLreadsMis"


#----------- TOTAL READS THAT MAPPED TO CHRC/M

#--- mapped on ChrM/C
nMaChrM=`grep "#id" $Rawstats | awk '{print $2,$3,$4,$5}' | awk '$1=="ChrM" {print $3}'`
nMaChrC=`grep "#id" $Rawstats | awk '{print $2,$3,$4,$5}' | awk '$1=="ChrC" {print $3}'`
nMaChrMC=$((nMaChrM+nMaChrC))
echo "$stamp: Mapped reads on ChrM+ChrC: $nMaChrMC (ChrM: $nMaChrM ChrC: $nMaChrC )"

#----------- TOTAL READS THAT MAPPED TO CHR1-5
#--- mapped on  Chr1-5
nMaChr15=`grep "#id" $Rawstats | awk '{print $2,$3,$4,$5}' |awk '$1!="ChrC" && $1!="ChrM" {s+=$3}END{print s}'`
echo "$stamp: Mapped reads on Chr1-->Chr5: $nMaChr15"


#----------- TOTAL READS SELECTED READS

FinalReads=`$pathSamtools view $ID.${maxMis}mis.mkdup.f3F1024.masked.shifted.bam -c`
echo "$stamp: Final number of selected reads: $FinalReads"

#----------------------------------------------- print all
title="SAMPLE,Total_reads_fastq,Total_mapped_reads,Mapped_on_ChM,Mapped_on_ChC,Mapped_on_ChM-ChrC,Mapped_on_Chr1-5,Total_mapped_selected_reads,Total_mapped_selected_reads_${maxMis}mis,Total_mapped_selected_reads_${maxMis}mis_dup,Final_reads"

outcsv=$OUTDIR/$ID.filter.stats.csv

cat > $outcsv <<EOF
#Total_reads_fastq: total sequenced read
#Total_mapped_reads: total mapped reads
#Mapped_on_ChM: number of reads that mapped to ChrM
#Mapped_on_ChC:number of reads that mapped to ChrC
#Mapped_on_ChM-ChrC: sum of number of reads that mapped to ChrM and ChC
#Mapped_on_Chr1-5: number of reads that mapped to nuclear DNA (Chr1->Chr5)
#Total_mapped_selected_reads: total mapped reads that fall outside of pre-defined blackisted regions (blackilisted regions include ChrM-ChrC)
#Total_mapped_selected_reads_${maxMis}mis: total mapped reads that fall outside of pre-defined blackisted regions and carry no more than $maxMis mismatches
#Total_mapped_selected_reads_${maxMis}mis_dup: total mapped duplicate reads that fall outside of pre-defined blackisted regions and carry no more than $maxMis mismatches
#Final_reads: Number of reads selected after filtering steps based on (i) blacklisted regions (ii) mismatches (iii) non-duplicated concordant pairs
,,
$title
$ID,$totalFastqReads,$totalMa,$nMaChrM,$nMaChrC,$nMaChrMC,$nMaChr15,$SLreads,$SLreadsMis,$nDup,$FinalReads
EOF
echo "$stamp: Wrote filtering stat in CSV format: $outcsv"
fi



if [[ $readCoverage == "yes" ]]; then 
	if [[ "$map" == "yes" && "$filter"=="no" ]]; then BAM=${ID}.mapped.sorted.bam ; fi
	if [[ "$filter" == "yes" ]]; then BAM=$ID.${maxMis}mis.mkdup.f3F1024.masked.shifted.bam ; fi
#===============================================================================
# Compute read coverage 
#===============================================================================
echo -e "\n######################## Compute read coverage (using BAM:$BAM)" | tee -a $LOG

#outputs
BEDGRAPH=${BAM/.bam/.bedgraph}
TDF=${BAM/.bam/.tdf}


echo "$stamp: Get bedgraph..."
$pathGenomeCoverageBed -ibam $BAM -g $CHRLEN -bga -trackline > $BEDGRAPH
echo "$stamp: Convert bedgraph to tdf..."
$pathIgvTools toTDF $BEDGRAPH $TDF $GENOME 1>> $LOG 2>&1
#remove bedgraph if tdf was correctly created: -s means file exists and has a non zero size
if [ -s $TDF ]; then rm -f $BEDGRAPH ; fi 
rm -f igv.log

#mv output to Coverage dir
mv $TDF $COVDIR
echo "$stamp: Wrote coverage file: $COVDIR/$TDF "

fi



if [[ $ieventsCoverage == "yes" ]]; then 
if [[ "$map" == "yes" && "$filter" == "no" ]]; then BAM=${ID}.mapped.sorted.bam ; fi
if [[ "$filter" == "yes" ]]; then BAM=$ID.${maxMis}mis.mkdup.f3F1024.masked.shifted.bam ; fi

#===============================================================================
# Compute insertion event coverage
#===============================================================================
echo -e "\n######################## Create insertion event BAM files and compute ievent coverage (using BAM:$BAM)" | tee -a $LOG


#outputs
OUTBAM=${BAM/.bam/.ievent.bam}
BEDGRAPH=${BAM/.bam/.ievent.bedgraph}
TDF=${BAM/.bam/.ievent.tdf}

echo "$stamp: Extract insertion events..." | tee -a $LOG

$pathSamtools view -h $BAM -F 16 |awk 'BEGIN{OFS="\t"}{ if ($1 ~ /^@/) {print ;} else {Rlen=length($10) ; $6=1"M"; $8=$8+Rlen-1 ; $10=substr($10,1,1);$11=substr($11,1,1); print $0 }}' > $ID.tmp.ievent.sam
$pathSamtools view $BAM -f 16 | awk 'BEGIN{OFS="\t"}{Rlen=length($10) ; $6=1"M"; $4=$4+Rlen-1; $10=substr($10,Rlen,1) ; $11=substr($11,Rlen,1) ; print $0}' >> $ID.tmp.ievent.sam


$pathSamtools view $ID.tmp.ievent.sam -bh | $pathSamtools sort - -o $OUTBAM -O bam -T $LOCALTMP
$pathSamtools index $OUTBAM
if [ -s $OUTBAM ]; then rm -f $ID.tmp.ievent.sam ; fi 

echo "$stamp: Wrote insertion events in bam file: $OUTBAM "


#then, get the coverage of insertion events
#cd $IEDIR
echo "$stamp: Compute ievents coverage..." | tee -a $LOG
echo "$stamp: Get bedgraph..."
$pathGenomeCoverageBed -ibam $OUTBAM -g $CHRLEN -bga -trackline > $BEDGRAPH
echo "$stamp: Convert bedgraph to tdf... "
$pathIgvTools toTDF $BEDGRAPH $TDF $GENOME 1>> $LOG 2>&1
if [ -s $TDF ]; then rm -f $BEDGRAPH ; fi  
rm -f igv.log

#mv output to Coverage dir
mv $TDF $COVDIR
echo "$stamp: Wrote coverage of insertion events: $COVDIR/$TDF "

fi



if [[ "$fragDist" == "yes" ]]; then 

#get correct BAM file (user-provided or created)
if [[ "$map" == "yes" && "$filter" == "no" ]]; then BAM=${ID}.mapped.sorted.bam ; fi
if [[ "$filter" == "yes" ]]; then BAM=$ID.${maxMis}mis.mkdup.f3F1024.masked.shifted.bam ; fi

#===============================================================================
# Compute fragment length distribution
#===============================================================================
echo -e "\n######################## Compute fragment length distribution (using BAM:$BAM)" | tee -a $LOG

#create dir for fragment length
if [ ! -d $OUTDIR/Fragment_distribution_${ID} ]; then mkdir $OUTDIR/Fragment_distribution_${ID} ; fi 

echo "$stamp: Extraction of fragment length... " | tee -a $LOG
#meaning of flag 66 -> 64 (first in pair) + 2 (read mapped in proper pair. Because we still have reads with no 2nd mate (*/=) and fixmate puts a TLEN=0 in this case)
echo "samtools view $BAM -f 66 |  awk -v tot=$LOCALTMP/tot.tmp ' function abs(v) {return v < 0 ? -v : v}  {s++; print abs($9)} END {print s > tot }' | sort - -T $LOCALTMP | uniq -c |sort -k2 -g > TLEN.$ID.f66.txt " >> $LOG
samtools view $BAM -f 66 |  awk -v tot=$LOCALTMP/tot.tmp ' function abs(v) {return v < 0 ? -v : v}  {s++; print abs($9)} END {print s > tot }' | sort - -T $LOCALTMP | uniq -c |sort -k2 -g > TLEN.$ID.f66.txt 

#Convert counts to frequencies:
tot=`awk '{print $1}' $LOCALTMP/tot.tmp`
awk -v tot=$tot 'BEGIN{OFS="\t"; print "#TLEN\tcount\tFreq"} {printf ("%d\t%d\t%f\n",$2,$1,$1/tot)}'  TLEN.$ID.f66.txt > tmp ; mv tmp TLEN.$ID.f66.txt

echo "$stamp: Plot fragment length distribution... "
gnuplot << EOF
set term pngcairo size 1024,768
set output "TLEN.$ID.f66.png"
set nokey

set datafile commentschars "#"
#set autoscale
 
set title "Fragment length distribution ($ID)"
set xlabel "Fragment length"
set ylabel "Frequencies"
set grid
 
set xrange [0:600]
set xtics 50
 
set grid mxtics  # draw lines for each xtics and mytics
set mxtics 10    # set the spacing for the mytics
set grid         # enable the grid


#set arrow from 2,graph(0,0) to 2,graph(1,1) nohead
#plot 'TLEN.$ID.f66.txt' title "" with points pt 7 ps 0.1

plot 'TLEN.$ID.f66.txt' using 1:3 with lines
EOF


#move outputs to directory
mv  TLEN.$ID.f66.txt TLEN.$ID.f66.png $FRAGDIR
echo "$stamp: Wrote $FRAGDIR/TLEN.$ID.f66.txt ; $FRAGDIR/TLEN.$ID.f66.png"
fi





if [[ $callpeak == "yes" ]]; then 
#get correct BAM file (user-provided or created)
if [[ "$map" == "yes" && "$filter" == "no" ]]; then BAM=${ID}.mapped.sorted.bam ; fi
if [[ "$filter" == "yes" ]]; then BAM=$ID.${maxMis}mis.mkdup.f3F1024.masked.shifted.bam ; fi
#===============================================================================
# Peak calling
#===============================================================================
OUTPC=$OUTDIR/Peak_calling_${ID}
if [ ! -d $OUTPC ]; then mkdir  $OUTPC ; fi
echo -e "\n######################## MACS2 Peak calling (using BAM:$BAM)" | tee -a $LOG
FORMAT=BAM 

#commands with and without control for both broad and narrow mode
if [ ! -z $peakCallingControl ]; then
	if [[ "$MODE" == "broad" ]]; then 

echo "$pathMACS2 callpeak -c $peakCallingControl -t $BAM -n "${ID}_${MODE}_nomodel" -g $gsize  -f $FORMAT  --outdir $OUTPC --keep-dup all --nomodel --broad -B --trackline -q $fdr --broad-cutoff $fdr" >> $LOG
$pathMACS2 callpeak -c $peakCallingControl -t $BAM -n "${ID}_${MODE}_nomodel" -g $gsize  -f $FORMAT  --outdir $OUTPC --keep-dup all --nomodel --broad -B --trackline -q $fdr --broad-cutoff $fdr

	elif [[ "$MODE" == "narrow" ]]; then 

echo "$pathMACS2 callpeak -c $peakCallingControl -t $BAM -n "${ID}_${MODE}_nomodel" -g $gsize  -f $FORMAT  --outdir $OUTPC --keep-dup all --nomodel --call-summits -B --trackline -q $fdr" >> $LOG
$pathMACS2 callpeak -c $peakCallingControl -t $BAM -n "${ID}_${MODE}_nomodel" -g $gsize  -f $FORMAT  --outdir $OUTPC --keep-dup all --nomodel --call-summits -B --trackline -q $fdr 
	fi

#no control
else
	if [[ "$MODE" == "broad" ]]; then  

echo "$pathMACS2 callpeak -t $BAM -n "${ID}_${MODE}_nomodel" -g $gsize  -f $FORMAT  --outdir $OUTPC --keep-dup all --nomodel --broad -B --trackline -q $fdr --broad-cutoff $fdr " >> $LOG
$pathMACS2 callpeak -t $BAM -n "${ID}_${MODE}_nomodel" -g $gsize  -f $FORMAT  --outdir $OUTPC --keep-dup all --nomodel --broad -B --trackline -q $fdr --broad-cutoff $fdr 

	elif [[ "$MODE" == "narrow" ]]; then 

echo "$pathMACS2 callpeak  -t $BAM -n "${ID}_${MODE}_nomodel" -g $gsize  -f $FORMAT  --outdir $OUTPC --keep-dup all --nomodel --call-summits -B --trackline -q $fdr " >> $LOG
$pathMACS2 callpeak  -t $BAM -n "${ID}_${MODE}_nomodel" -g $gsize  -f $FORMAT  --outdir $OUTPC --keep-dup all --nomodel --call-summits -B --trackline -q $fdr 
	fi
fi


fi
#============================================================================

echo -e "\t=================================================================================\n\tUSER: $USER\n\tSample: $ID\n\tASAP end time: `date`\n\t================================================================================="

