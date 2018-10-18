#!/bin/bash 


:<<'readme'
11/10/2018
kramdi@biologie.ens.fr
kramdi.a@gmail.com

This script is written in an effort to gather the main analysis steps of ATAC-seq data.Here are the major steps:

I.Mapping
II.Post-mapping processing (mark duplicate pairs) and filtering
 II.1. Filter (or not) reads that fall into user-defined blacklisted regions
 II.2. Select reads that do not carry more than minMismatch. Filter by mapping quality (MAPQ)
 II.3. Select concordant, non-duplicated pairs. 
 II.4. Shift read by 4bp as described in Schep et al.,2015
III. Compute coverage
 III.1. Compute read coverage
 III.2. Compute insertion events coverage
IV. Compute fragment length distribution
 V. Extract read pairs based on a given range of fragment length and compute arcs between fragments extremities ( protection visualization)
VI. Peak calling

readme

#Here, we get the directory of this script in order to display a proper help
#SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
VERSION="v2.0"

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

 
while getopts "hvc:" OPTION
do
     case $OPTION in
		v)
		echo "ASAP - version $VERSION"; exit 0
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
stamp() {
date=`date "+[%b-%d-%Y %T]"`
echo "@$date"
}



#===============================================================================
# Parse configuration file to extract parameters 
#===============================================================================

#Initialize argument array
typeset -A config # init array
config=( # set default values in config array
[OUTDIR]=""
[LOCALTMP]=""
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
[mapQual]=""
[readCoverage]=""
[ieventsCoverage]=""
[callpeak]=""
[fragDist]=""
[control]=""
[CHRLEN]=""
[GENOME]=""
[blacklist]=""
[blacklistedRegions]=""
[shift]=""
[pathSamtools]=""
[pathIgvTools]=""
[pathGenomeCoverageBed]=""
[pathMACS2]=""
[pathBowtie2]=""
[pathToJava]=""
[pathToPicardJar]=""
[extractReads]=""
[lowBoundary]=""
[upBoundary]=""
[arcs]=""
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
LOCALTMP=${config[TMPDIR]} #changed variable name
ID=${config[sampleName]} #changed variable name
FASTQ1=${config[FASTQ1]}
FASTQ2=${config[FASTQ2]}
bowtieIndex=${config[bowtieIndex]}
mappingParameters=${config[mappingParameters]}
maxMis=${config[maxMis]}
blacklist=${config[blacklist]}
blacklistedRegions=${config[blacklistedRegions]}
shift=${config[shift]}
gsize=${config[gsize]}
fdr=${config[fdr]}
MODE=${config[MODE]}
BAM=${config[BAM]}
map=${config[map]}
filter=${config[filter]}
mapQual=${config[mapQual]}
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
extractReads=${config[extractReads]}
lowBoundary=${config[lowBoundary]}
upBoundary=${config[upBoundary]}
arcs=${config[arcs]}

#===============================================================================
# Create tmp dir and LOG file
#===============================================================================

#remove possible tailing "/"
LOCALTMP=${LOCALTMP%\/}
if [ ! -d $LOCALTMP ]; then >&2 echo -e "\n[ERROR]: TMPDIR must be a directory. Exit." ; exit 1 ;fi

#create tmp dir
LOCALTMP=$LOCALTMP/tmp$RANDOM
if [ ! -d $LOCALTMP ]; then mkdir $LOCALTMP ; fi 

function finish {
rm -rf $LOCALTMP
}
trap finish EXIT 

#Create log file where all commands are kept. Comes in handy in case of debugging
if [ -d $OUTDIR ]; then >&2 echo "[WARNING]: Output directory ($OUTDIR) already exists. Results will be overwitten." ; else mkdir $OUTDIR; fi 
LOGDIR=$OUTDIR/logs_${ID} ; if [ ! -d $LOGDIR ] ; then mkdir $LOGDIR ;fi 
LOG=$LOGDIR/commands_ASAP_run_${ID}.log
touch $LOG

#================================== Check some parameters and issue warnings/errors

#check yes/no parameters
if [[ ( "$map" != "yes" && "$map" != "no" ) || ( "$filter" != "yes" && "$filter" != "no" ) || ( "$readCoverage" != "yes" && "$readCoverage" != "no" )  || ( "$ieventsCoverage" != "yes" && "$ieventsCoverage" != "no" ) || ( "$fragDist" != "yes" && "$fragDist" != "no" ) || ( "$callpeak" != "yes" && "$callpeak" != "no" ) || ( "$extractReads" != "yes" && "$extractReads" != "no" ) || ( "$blacklist" != "yes" && "$blacklist" != "no" ) || ( "$shift" != "yes" && "$shift" != "no" ) ]] ; then 
	>&2 echo -e "\n[ERROR]: Possible values to turn on a section/option: yes/no. Exit." ; exit 1 
fi



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
#COVDIR=$OUTDIR/Coverage_${ID}
#if [[ ( "$readCoverage" == "yes" || "$ieventsCoverage" == "yes" ) && ! -d $COVDIR ]]; then mkdir $COVDIR ;fi 

#Create fragment dist directory
#FRAGDIR=$OUTDIR/Fragment_distribution_${ID}
#if [[  "$fragDist" == "yes" && ! -d $FRAGDIR ]]; then mkdir $FRAGDIR ;fi 



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
- Use blacklisted regions to filter reads: $blacklist
EOL
if [[ "$blacklist" == "yes" ]]; then
cat >> $LOCALTMP/tmp.conf <<EOL
 - Blacklisted regions: $blacklistedRegions
EOL
fi
cat >> $LOCALTMP/tmp.conf <<EOL
- Maximum mismatches per read: $maxMis
- Miinimum MAPQ:$mapQual
- Shift reads by 4bp: $shift
- Select concordant, non-duplicated pairs
EOL
fi


#write coverage section/fragment distribution section/peak calling section + options
cat >> $LOCALTMP/tmp.conf <<EOL
# Compute read coverage: $readCoverage
# Compute insertion events coverage: $ieventsCoverage
# Compute fragment length distribution: $fragDist
# Extract read pairs based on fragment length: $extractReads
EOL

if [[ "$extractReads" == "yes" ]]; then 
cat >> $LOCALTMP/tmp.conf <<EOL
- Lower fragment length boundary: $lowBoundary nt
- Upper fragment length boundary: $upBoundary nt
- Compute arcs between fragment extremities: $arcs
EOL
fi

cat >> $LOCALTMP/tmp.conf <<EOL
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


########################################################################################### MAIN #################################################################################"
cd $OUTDIR


if [[ $map == "yes" ]]; then
if [ ! -s $pathBowtie2 ]; then >&2 echo `file $pathBowtie2`  ; exit 1 ; fi
#===============================================================================
# Mapping step
#===============================================================================
echo -e "\n######################## Mappinng step"
echo -e "\n######################## Mappinng step" > $LOG

echo "`stamp`: Mapping..." | tee -a $LOG
echo "`stamp`: Mapping command: #$pathBowtie2 -x $bowtieIndex $mappingParameters -1 $FASTQ1 -2 $FASTQ2 1> ${ID}.mapped.sam" | tee -a $LOG
eval "$pathBowtie2 -x $bowtieIndex $mappingParameters -1 $FASTQ1 -2 $FASTQ2 1> ${ID}.mapped.sam "


echo "`stamp`: Convert mapping output to bam format..." | tee -a $LOG
echo "#$pathSamtools view -bhS ${ID}.mapped.sam -o ${ID}.mapped.bam" >> $LOG
$pathSamtools view -bhS ${ID}.mapped.sam -o ${ID}.mapped.bam
 

echo "`stamp`: Sort bam..." | tee -a $LOG
echo "#$pathSamtools sort ${ID}.mapped.bam -o ${ID}.sorted.bam -O bam -T $LOCALTMP" >> $LOG
$pathSamtools sort  ${ID}.mapped.bam -o ${ID}.mapped.sorted.bam -O bam -T $LOCALTMP

if [ -s ${ID}.mapped.sorted.bam ]; then rm -f ${ID}.mapped.sam ${ID}.mapped.bam ; fi 

echo "`stamp`: Index sorted bam..." | tee -a $LOG
echo "#$pathSamtools index ${ID}.mapped.sorted.bam" >> $LOG
$pathSamtools index ${ID}.mapped.sorted.bam

#if mapping step is skipped. Index the provided BAM file
else

	echo "Input BAM: $BAM"
	if [ ! -s $BAM.bai ]; then echo "`stamp`: Index input BAM..." ; $pathSamtools index $BAM ; fi
fi




if [[ "$filter" == "yes" ]]; then 
if [ ! -s $pathSamtools ]; then >&2 echo `file $pathSamtools`  ; exit 1 ; fi
if [ ! -s $pathToJava ]; then >&2 echo `file $pathToJava`  ; exit 1 ; fi
if [ ! -s $pathToPicardJar ]; then >&2 echo `file $pathToPicardJar`  ; exit 1 ; fi

if [[ "$map" == "yes" ]]; then BAM=${ID}.mapped.sorted.bam ;fi 
#===============================================================================
# Post-processing and filtering
#===============================================================================
echo -e "\n######################## Post-processing and filtering step (using BAM:$BAM)" | tee -a $LOG

#this file will hold the number of reads before and after the mismatch filter for informative purposes
IDnreads=$LOCALTMP/$ID.nreads.txt

# Remove reads based on blacklisted regions (if provided) + filter by max N mismatches
if [[ "$blacklist" == "yes" ]]; then
##create complement of the blacklisted regions for easy selection using samtools -L
bedtools complement -i <(sort -k1,1 -k2,2n $blacklistedRegions) -g <(grep -v 'ChrC\|ChrM' $CHRLEN) | sort -k1,1 -k2,2n >  $LOCALTMP/selected_ATACseq_regions.bed

#Select reads that fall into our selected regions + filter by N mis
mask="masked"
echo "`stamp`: Filter by blacklisted regions and allowed mismatches..." | tee -a $LOG
$pathSamtools view -h -F 4 $BAM -L $LOCALTMP/selected_ATACseq_regions.bed | awk -v n=$maxMis -v ID=$ID -v IDnreads=$IDnreads  '{if($0 ~ /^@/) {print $0} else { a++ ; mis=match($0,"XM:i:") ; nmis=substr($0,mis+5,2) ;if (int(nmis)<=n) {print $0; b++} } }END{printf("selected-regions\t%d\nreads-selected-by-nmis\t%d\n",a,b) > IDnreads }'  | $pathSamtools sort - -o $ID.${maxMis}mis.$mask.sorted.bam -O bam -T $LOCALTMP
else
#if no blacklist is provided, do: filter by N mis
mask="unmasked"
echo "`stamp`: Filter by allowed mismatches..." | tee -a $LOG
$pathSamtools view -h -F 4 $BAM | awk -v n=$maxMis -v ID=$ID -v IDnreads=$IDnreads  '{if($0 ~ /^@/) {print $0} else {a++ ;  mis=match($0,"XM:i:") ; nmis=substr($0,mis+5,2) ;if (int(nmis)<=n) {print $0; b++} } }END{printf("selected-regions\t%d\nreads-selected-by-nmis\t%d\n",a,b) > IDnreads }'  | $pathSamtools sort - -o $ID.${maxMis}mis.$mask.sorted.bam -O bam -T $LOCALTMP
fi

#check if duplicates are already marked: saves time when filtering is performed on already marked reads.. 
if [ `samtools view $ID.${maxMis}mis.$mask.sorted.bam -H | grep "MarkDuplicates" -c` -gt 0 ]; then 
#skip marking duplicates and just rename file
echo "`stamp`: BAM header says that duplicate reads are already marked (skip step)..." | tee -a $LOG
mv $ID.${maxMis}mis.$mask.sorted.bam $ID.${maxMis}mis.mkdup.$mask.bam 

else

echo "`stamp`: Mark duplicates..." | tee -a $LOG
echo "#$pathToJava -Xmx5g -Djava.io.tmpdir=$LOCALTMP \
-jar $pathToPicardJar MarkDuplicates \
METRICS_FILE=$LOGDIR/MKDUP.$ID.${maxMis}mis.$mask.log  \
INPUT=$ID.${maxMis}mis.$mask.sorted.bam \
OUTPUT=$ID.${maxMis}mis.mkdup.$mask.bam 2> $log/MKDUP.$ID.${maxMis}mis.$mask.err" >> $LOG


$pathToJava -Xmx5g -Djava.io.tmpdir=$LOCALTMP \
-jar $pathToPicardJar MarkDuplicates \
METRICS_FILE=$LOGDIR/MKDUP.$ID.${maxMis}mis.$mask.log  \
INPUT=$ID.${maxMis}mis.$mask.sorted.bam \
OUTPUT=$ID.${maxMis}mis.mkdup.$mask.bam 2> $LOGDIR/MKDUP.$ID.${maxMis}mis.$mask.err
fi


#Get number of duplicates (from the reads falling into our selected regions and carrying no more than nmis) 
#echo "`stamp`: Get stats on duplicates..." 
#nDup=`$pathSamtools view -f 1024 $ID.${maxMis}mis.mkdup.$mask.bam -c`

echo "`stamp`: Select concordant, non-duplicated reads and fix flags..."   | tee -a $LOG
echo "#$pathSamtools view -hb -f 3 -F 1024 -q $mapQual $ID.${maxMis}mis.mkdup.$mask.bam >  $ID.${maxMis}mis.mkdup.f3F1024.$mask.bam  " >> $LOG
 
$pathSamtools view -hb -f 3 -F 1024 -q $mapQual $ID.${maxMis}mis.mkdup.$mask.bam >  $ID.${maxMis}mis.mkdup.f3F1024.$mask.bam  #get concordant, non-dup, mapping quality >= mapQual

#remove intermediate files
if [ -s  $ID.${maxMis}mis.mkdup.f3F1024.$mask.bam ] ; then rm -f $ID.${maxMis}mis.mkdup.$mask.bam $ID.${maxMis}mis.$mask.sorted.bam ; fi

#========4. Fix mates!
$pathSamtools sort -n $ID.${maxMis}mis.mkdup.f3F1024.$mask.bam -T $LOCALTMP -O bam | $pathSamtools fixmate - $ID.tmp.bam

#========5. Shift reads (if yes), adjust fragment length (-8) and sort by coordinates for indexing; otherwise, carry on 
shiftTag="shifted"
if [[ $shift == "yes" ]]; then 
echo "`stamp`: Shift reads..." 
#use -f3 because reads are fixed (fixmate) we have to extract the new concordant reads
$pathSamtools view -h -f 3 -F 16 $ID.tmp.bam | awk '  function abs(v) {return v < 0 ? -v : v} BEGIN{OFS="\t"} { if ($1 ~ /^@/) {print} else { Rlen=length($10) ; $4=$4+4 ; $8=$8-4 ; if ($8<0 && $8==0) {$8=1} ;  if ($9>0) {i=1} else if ($9<0) {i=-1} else {i=0}; $9=(abs($9)-8)*i;  print $0   }}' > $ID.tmp.$shiftTag.sam
$pathSamtools view -f 19 $ID.tmp.bam | awk ' function abs(v) {return v < 0 ? -v : v}  BEGIN{OFS="\t"} { $4=$4-4 ;if ($4<0 || $4==0) {$4=1} ;  $8=$8+4 ;  if ($9>0) {i=1} else if ($9<0) {i=-1} else {i=0}; $9=(abs($9)-8)*i; print $0 }' >> $ID.tmp.$shiftTag.sam
$pathSamtools view -bh $ID.tmp.$shiftTag.sam | $pathSamtools sort - -o $ID.${maxMis}mis.mkdup.f3F1024.$mask.$shiftTag.bam -O bam -T $LOCALTMP
else
echo "`stamp`: Skip shiftting reads (shift=$shift), sort reads..." 
shiftTag="unshifted"
$pathSamtools view -f 3 $ID.tmp.bam  -h | $pathSamtools sort - -o $ID.${maxMis}mis.mkdup.f3F1024.$mask.$shiftTag.bam -O bam -T $LOCALTMP
fi

# index file
echo "`stamp`: Index output bam file..."
$pathSamtools index $ID.${maxMis}mis.mkdup.f3F1024.$mask.$shiftTag.bam

#remove tmp bam file
if [ -s $ID.${maxMis}mis.mkdup.f3F1024.$mask.$shiftTag.bam ]; then rm -f $ID.tmp.bam $ID.tmp.$shiftTag.sam  $ID.${maxMis}mis.mkdup.f3F1024.$mask.bam ; fi 

echo "`stamp`: Extract filtering stats in CSV file..."

Rawstats=$LOCALTMP/samtools.stats.inBAM.$ID.txt ; touch $Rawstats 
$pathSamtools flagstat $BAM | awk 'BEGIN{OFS="\t"}{print "#flagstat", $0}' > $Rawstats
$pathSamtools idxstats $BAM | awk 'BEGIN{OFS="\t"}{print "#idx", $0}'  >> $Rawstats
#write raw samtools stats in LOG 
echo "# RAW SAMTOOLS STATS" >> $LOG ; cat $Rawstats >> $LOG

totalFastqReads=`cat $Rawstats | grep "#flagstat"  | awk '{print $2}' | awk '{FS="+"} NR ==1 {print $1}'`
echo "`stamp`: Total reads in $ID fastq: $totalFastqReads" 

#total unmapped/mapped (no mismatche filtering)
totalMa=`cat $Rawstats | grep "#flagstat" | awk '{print $2}' | awk '{FS="+"} NR ==5 {print $1}'`
totalUnma=$((totalFastqReads-totalMa)) 
echo "`stamp`: Unmapped reads: $totalUnma"
echo "`stamp`: Mapped reads: $totalMa"


#To get the number of reads falling into blacklisted regions, we subtract the number of reads that DID not fall into the black listed regions - total number of reads
#BLreads=`echo "$totalFastqReads-`grep selected $IDnreads | cut -f2`" | bc`

#SLreads=`grep regions $IDnreads | cut -f2`
#echo "`stamp`: Reads falling into selected regions: $SLreads "

#SLreadsMis=`grep nmis $IDnreads | cut -f2`
#echo "`stamp`: Reads selected by $maxMis mismatch(es) (among selected reads by blacklisted regions): $SLreadsMis"


#----------- TOTAL READS THAT MAPPED TO CHRC/M

#--- mapped on ChrM/C
nMaChrM=`grep "#id" $Rawstats | awk '{print $2,$3,$4,$5}' | awk '$1=="ChrM" {print $3}'`
nMaChrC=`grep "#id" $Rawstats | awk '{print $2,$3,$4,$5}' | awk '$1=="ChrC" {print $3}'`
nMaChrMC=$((nMaChrM+nMaChrC))
echo "`stamp`: Mapped reads on ChrM+ChrC: $nMaChrMC (ChrM: $nMaChrM ChrC: $nMaChrC )"

#----------- TOTAL READS THAT MAPPED TO CHR1-5
#--- mapped on  Chr1-5
nMaChr15=`grep "#id" $Rawstats | awk '{print $2,$3,$4,$5}' |awk '$1!="ChrC" && $1!="ChrM" {s+=$3}END{print s}'`
echo "`stamp`: Mapped reads on Chr1-->Chr5: $nMaChr15"


#----------- TOTAL READS SELECTED READS

FinalReads=`$pathSamtools view $ID.${maxMis}mis.mkdup.f3F1024.$mask.$shiftTag.bam -c`
echo "`stamp`: Final number of selected reads: $FinalReads"

#----------------------------------------------- print all
#title="SAMPLE,Total_reads_fastq,Total_mapped_reads,Mapped_on_ChM,Mapped_on_ChC,Mapped_on_ChM-ChrC,Mapped_on_Chr1-5,Total_mapped_selected_reads,Total_mapped_selected_reads_${maxMis}mis,Total_mapped_selected_reads_${maxMis}mis_dup,Final_reads"

title="SAMPLE,Total_reads_fastq,Total_mapped_reads,Mapped_on_ChM,Mapped_on_ChrC,Mapped_on_ChM-ChrC,Mapped_on_Chr1-5,selected_reads"

outcsv=$OUTDIR/$ID.filter.stats.$mask.$shiftTag.csv

cat > $outcsv <<EOF
#Mapping parameters: bowtie2 $mappingParameters
#Blacklist: $blacklist
,,
#Total_reads_fastq: total sequenced read
#Total_mapped_reads: total mapped reads
#Mapped_on_ChM: number of reads that mapped to ChrM
#Mapped_on_ChrC:number of reads that mapped to ChrC
#Mapped_on_ChM-ChrC: sum of number of reads that mapped to ChrM and ChrC
#Mapped_on_Chr1-5: number of reads that mapped to nuclear DNA (Chr1->Chr5)
#selected_reads: Number of reads selected after filtering steps based on (i) blacklisted regions if set (ii) mismatches, MAPQ (iii) non-duplicated concordant pairs
,,
$title
$ID,$totalFastqReads,$totalMa,$nMaChrM,$nMaChrC,$nMaChrMC,$nMaChr15,$FinalReads
EOF
echo "`stamp`: Wrote filtering stat in CSV format: $outcsv"
fi


#removed categories
#Total_mapped_selected_reads: Number of selected reads after blacklist filtering. (Equals the initial number of reads if no blacklist was provided)
#Total_mapped_selected_reads_${maxMis}mis: Number of selected reads based on maximum $maxMis mismatches
#Total_mapped_selected_reads_${maxMis}mis_dup: Number of duplicates (reads are counted after the two previous filtering steps

#$SLreads,$SLreadsMis,$nDup,


#===============================================================================

if [[ $readCoverage == "yes" ]]; then 

if [ ! -s $pathIgvTools ]; then >&2 echo `file $pathIgvTools`  ; exit 1 ; fi

	if [[ "$map" == "yes" && "$filter"=="no" ]]; then BAM=${ID}.mapped.sorted.bam ; fi
	if [[ "$filter" == "yes" ]]; then BAM=$ID.${maxMis}mis.mkdup.f3F1024.$mask.$shiftTag.bam ; fi
#===============================================================================
# Compute read coverage 
#===============================================================================
echo -e "\n######################## Compute read coverage (using BAM:$BAM)" | tee -a $LOG

#outputs
base=`basename $BAM`
BEDGRAPH=${base/.bam/.bedgraph}
TDF=${base/.bam/.tdf}


echo "`stamp`: Get bedgraph..."
$pathGenomeCoverageBed -ibam $BAM -bga -trackline > $BEDGRAPH
echo "`stamp`: Convert bedgraph to tdf..."
$pathIgvTools toTDF $BEDGRAPH $TDF $GENOME 1>> $LOG 2>&1
#remove bedgraph if tdf was correctly created: -s means file exists and has a non zero size
if [ -s $TDF ]; then rm -f $BEDGRAPH ; fi 
rm -f igv.log

#mv output to Coverage dir
#mv $TDF $COVDIR
#echo "`stamp`: Wrote coverage file: $COVDIR/$TDF "
echo "`stamp`: Wrote coverage file: $TDF "


fi

#===============================================================================

if [[ $ieventsCoverage == "yes" ]]; then 

if [ ! -s $pathIgvTools ]; then >&2 echo `file $pathIgvTools`  ; exit 1 ; fi

if [[ "$map" == "yes" && "$filter" == "no" ]]; then BAM=${ID}.mapped.sorted.bam ; fi
if [[ "$filter" == "yes" ]]; then BAM=$ID.${maxMis}mis.mkdup.f3F1024.$mask.$shiftTag.bam ; fi

#===============================================================================
# Compute insertion event coverage
#===============================================================================
echo -e "\n######################## Create insertion event BAM files and compute ievent coverage (using BAM:$BAM)" | tee -a $LOG


#outputs
base=`basename $BAM`
OUTBAM=${base/.bam/.ievent.bam}
BEDGRAPH=${base/.bam/.ievent.bedgraph}
TDF=${base/.bam/.ievent.tdf}

echo "`stamp`: Extract insertion events..." | tee -a $LOG

$pathSamtools view -h $BAM -F 16 |awk 'BEGIN{OFS="\t"}{ if ($1 ~ /^@/) {print ;} else {Rlen=length($10) ; $6=1"M"; $8=$8+Rlen-1 ; $10=substr($10,1,1);$11=substr($11,1,1); print $0 }}' > $ID.tmp.ievent.sam
$pathSamtools view $BAM -f 16 | awk 'BEGIN{OFS="\t"}{Rlen=length($10) ; $6=1"M"; $4=$4+Rlen-1; $10=substr($10,Rlen,1) ; $11=substr($11,Rlen,1) ; print $0}' >> $ID.tmp.ievent.sam


$pathSamtools view $ID.tmp.ievent.sam -bh | $pathSamtools sort - -o $OUTBAM -O bam -T $LOCALTMP
$pathSamtools index $OUTBAM
if [ -s $OUTBAM ]; then rm -f $ID.tmp.ievent.sam ; fi 

echo "`stamp`: Wrote insertion events in bam file: $OUTBAM "


#then, get the coverage of insertion events
#cd $IEDIR
echo "`stamp`: Compute ievents coverage..." | tee -a $LOG
echo "`stamp`: Get bedgraph..."
$pathGenomeCoverageBed -ibam $OUTBAM -bga -trackline > $BEDGRAPH
echo "`stamp`: Convert bedgraph to tdf... "
$pathIgvTools toTDF $BEDGRAPH $TDF $GENOME 1>> $LOG 2>&1
if [ -s $TDF ]; then rm -f $BEDGRAPH ; fi  
rm -f igv.log

#mv output to Coverage dir
#mv $TDF $COVDIR
echo "`stamp`: Wrote coverage of insertion events: $TDF "

fi

#===============================================================================

if [[ "$fragDist" == "yes" ]]; then 
if [ ! -s $pathSamtools ]; then >&2 echo `file $pathSamtools`  ; exit 1 ; fi
#get correct BAM file (user-provided or created)
	if [[ "$map" == "yes" && "$filter" == "no" ]]; then BAM=${ID}.mapped.sorted.bam ; fi
	if [[ "$filter" == "yes" ]]; then BAM=$ID.${maxMis}mis.mkdup.f3F1024.$mask.$shiftTag.bam ; fi

#===============================================================================
# Compute fragment length distribution
#===============================================================================
echo -e "\n######################## Compute fragment length distribution (using BAM:$BAM)" | tee -a $LOG

#create dir for fragment length
#if [ ! -d $OUTDIR/Fragment_distribution_${ID} ]; then mkdir $OUTDIR/Fragment_distribution_${ID} ; fi 

echo "`stamp`: Extraction of fragment length... "

#names for fragment dist txt file/png file
base=`basename $BAM`
TLENTXT=${base/.bam/.TLEN.f3F16.txt}
TLENPNG=${base/.bam/.TLEN.f3F16.png}


# meaning of flag -f 3 -F16 -> -f 3 (get pairs mapped in propoer pairs to avoid reads that do not have mates) + -F 16 (get forward reads out these selected pairs) 
# Get the absolute value of TLEN

echo "$pathSamtools view $BAM -f 3 -F 16 |  awk -v tot=$LOCALTMP/tot.tmp ' function abs(v) {return v < 0 ? -v : v}  {s++; print abs($9)} END {print s > tot }' | sort - -g -T $LOCALTMP | uniq -c |sort -k2 -g > $TLENTXT 
" >> $LOG
$pathSamtools view $BAM -f 3 -F 16 |  awk -v tot=$LOCALTMP/tot.tmp ' function abs(v) {return v < 0 ? -v : v}  {s++; print abs($9)} END {print s > tot }' | sort - -g -T $LOCALTMP | uniq -c |sort -k2 -g > $TLENTXT 

#Convert counts to frequencies:
tot=`awk '{print $1}' $LOCALTMP/tot.tmp`
awk -v tot=$tot 'BEGIN{OFS="\t"; print "#TLEN\tcount\tFreq"} {printf ("%d\t%d\t%f\n",$2,$1,$1/tot)}'  $TLENTXT > tmp ; mv tmp $TLENTXT


echo "`stamp`: Plot fragment length distribution... "
gnuplot << EOF
set term pngcairo size 1024,768
set output "$TLENPNG"
set nokey

set datafile commentschars "#"
#set autoscale
 
set title "Fragment length distribution ($ID)"
set xlabel "Fragment length"
set ylabel "Frequencies"
set grid
 
#set xrange [0:600]
set xtics 100
 
set grid mxtics  # draw lines for each xtics and mytics
set mxtics 10    # set the spacing for the mytics
set grid         # enable the grid

set xtics rotate by 315 offset -1,0

#set arrow from 2,graph(0,0) to 2,graph(1,1) nohead
#plot 'TLEN.$ID.f66.txt' title "" with points pt 7 ps 0.1

plot "$TLENTXT" using 1:3 with lines
EOF
#move outputs to directory
#mv  TLEN.$ID.f3F16.txt TLEN.$ID.f3F16.png $FRAGDIR
#echo "`stamp`: Wrote $FRAGDIR/TLEN.$ID.f3F16.txt ; $FRAGDIR/TLEN.$ID.f3F16.png"
echo "`stamp`: Wrote $TLENTXT ; $TLENPNG"
fi

#===============================================================================



if [[ "$extractReads" == "yes" ]]; then 
if [ ! -s $pathSamtools ]; then >&2 echo `file $pathSamtools`  ; exit 1 ; fi
#get correct BAM file (user-provided or created)
	if [[ "$map" == "yes" && "$filter" == "no" ]]; then BAM=${ID}.mapped.sorted.bam ; fi
	if [[ "$filter" == "yes" ]]; then BAM=$ID.${maxMis}mis.mkdup.f3F1024.$mask.$shiftTag.bam ; fi

#===============================================================================
# Extract read pairs based on a range of fragment length 
#===============================================================================
echo -e "\n######################## Extract reads based on fragment length range (using BAM:$BAM)" | tee -a $LOG
# Create output bam name:
base=`basename $BAM`
baseNoex=${base%.bam}

subSetSam=$LOCALTMP/$baseNoex.subReads.f3.frag-${lowBoundary}-${upBoundary}.sam
subSetBam=$baseNoex.subReads.f3.frag-${lowBoundary}-${upBoundary}.bam

header=$LOCALTMP/header.txt ; touch $header

echo "`stamp`: Extract properly paired reads (-f3) with fragment lengths [ $lowBoundary , $upBoundary ]..." | tee -a $LOG
echo "$pathSamtools view $BAM -H >  $LOCALTMP/header.txt" >> $LOG
echo "$pathSamtools view -f 3 $BAM | awk -v up=$upBoundary -v low=$lowBoundary -v out=$subSetSam ' function abs(v) {return v < 0 ? -v : v} $9 ~ /^[-0-9]+$/ {
if (abs($9) >= low && abs($9) <= up) {print $0 > out }}'" >> $LOG

$pathSamtools view $BAM -H >  $header
$pathSamtools view -f 3 $BAM | awk -v up=$upBoundary -v low=$lowBoundary -v out=$subSetSam ' function abs(v) {return v < 0 ? -v : v} $9 ~ /^[-0-9]+$/ {
if (abs($9) >= low && abs($9) <= up) {print $0 > out }}'

echo "`stamp`: Sort and index extracted reads..."
cat $header $subSetSam | $pathSamtools view -bh - | $pathSamtools sort - -O bam -o $subSetBam -T $LOCALTMP 
$pathSamtools index $subSetBam

N=`$pathSamtools view $subSetBam -c`
echo "`stamp`: Extracted $N reads in $subSetBam"

#compute reads coverage of extracted bam file
BEDGRAPH=${subSetBam/.bam/.bedgraph}
TDF=${subSetBam/.bam/.tdf}

echo "`stamp`: Get coverage file (.tdf)..." | tee -a $LOG
echo "$pathGenomeCoverageBed -ibam $subSetBam -bga -trackline > $BEDGRAPH" >> $LOG
echo "$pathIgvTools toTDF $BEDGRAPH $TDF $GENOME" >> $LOG

$pathGenomeCoverageBed -ibam $subSetBam -bga -trackline > $BEDGRAPH
$pathIgvTools toTDF $BEDGRAPH $TDF $GENOME 1>> $LOG 2>&1
if [ -s $TDF ]; then rm -f $BEDGRAPH ; fi  
rm -f igv.log


#compute insertion coverage of extracted bam file
iBEDGRAPH=${subSetBam/.bam/.ievent.bedgraph}
iTDF=${subSetBam/.bam/.ievent.tdf}

echo "`stamp`: Get insertion coverage file (.tdf)..." | tee -a $LOG
$pathGenomeCoverageBed -ibam  $subSetBam -bga -5 > $iBEDGRAPH
$pathIgvTools toTDF $iBEDGRAPH $iTDF $GENOME 1>> $LOG 2>&1
if [ -s $iTDF ]; then rm -f $iBEDGRAPH ; fi  
rm -f igv.log

#compute arcs between fragment extremities: useful for visualization of protection (nucleosomes) on IGV
if [[ "$arcs" == "yes" ]]; then 
arcFile=$baseNoex.subReads.f3.frag-${lowBoundary}-${upBoundary}.arcs.bed
echo "`stamp`: Compute arcs between fragment extremities..." | tee -a $LOG
echo 'track name=junctions description="fragment links from '$baseNoex'"' > $arcFile
#about the format: https://software.broadinstitute.org/software/igv/splice_junctions
#format: [seqname] [start] [end] [id] [score] [strand] [thickStart] [thickEnd] [r,g,b] [block_count] [block_sizes] [block_locations] #"frag"NR
samtools view -F 16 -f 3 $subSetBam | awk 'BEGIN{OFS="\t"} {Rlen=length($10); start=$4-1 ; end=$8+Rlen ; print $3,start,end, $1, 1,"+" }' | sort -k1,1 -k2n,2 -T $LOCALTMP  >> $arcFile
fi 


fi

#===============================================================================

if [[ $callpeak == "yes" ]]; then

if [ ! -s $pathMACS2 ]; then >&2 echo `file $pathMACS2`  ; exit 1 ; fi

#get correct BAM file (user-provided or created)
if [[ "$map" == "yes" && "$filter" == "no" ]]; then BAM=${ID}.mapped.sorted.bam ; fi
if [[ "$filter" == "yes" ]]; then BAM=$ID.${maxMis}mis.mkdup.f3F1024.$mask.$shiftTag.bam ; fi
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
