#!/bin/bash -l
##By: Hans Vasquez-Gross
set -e
set -u
set -o pipefail

##Written specifically for the tetraploid project and the structure in which our sample data is downloaded from BGI sequencing facility
##This script can be modified to run a similar pipeline for any project
##Assumes PE reads are in the following directory structure:  CaptureSet/DeMupltiplexed Sample Number/reads_1.fq.gz or reads_2.fq.gz

usage() { echo "Usage: $0 [-p <base_path> (no trailing slash)] inputcapturename{s}" 1>&2; exit 1; }
##Set default BasePath
base_path=/group/dubcovskygrp/projects/sequence/BGI_data

while getopts "hp:" opt; do
    case $opt in
        p) base_path=$OPTARG;;
        h) usage;;
    esac
done
shift $(( $OPTIND-1 ))
(($# >= 1)) || { usage; }

#####
##Trimming Configuration
####
OUT_BASEDIR=/group/dubcovskygrp/projects/sequence/processed/trimmed-reads
ADAPTORS_DIR=/group/dubcovskygrp/projects/sequence/BGI_data
## presets
PRIOR=0.4
QUAL_THRESH=20
export TMP_DIR=/home/hansvgdub/scratch/$SLURM_JOBID
##############
####END
#############

#########################
##BWA Mapping Configuration
########################
SCRIPTS_DIR=/group/dubcovskygrp/scripts
REF_DIR=/group/dubcovskygrp/projects/reference
REF=$REF_DIR/IWGSC_CSS_AB-TGAC_UCW_v1.fa
today=$(date +%Y%m%d)
#today="20140511" ##if rerunning an old pipeline another day
export TMP_DIR=/home/hansvgdub/scratch/$SLURM_JOBID
####################
##END CONFIGURATION$
####################

##Log Start time
date
hostname
module load sickle
module load scythe
module load seqqs 
module load bwa
module load samtools
module load picardtools

##Setup a function
echoerr() { echo "$@" 1>&2; }

if [ ! -e $TMP_DIR ]; then
    echo "Making $TMP_DIR"
    mkdir $TMP_DIR
fi
echo "created and assigned TMP_DIR to $TMP_DIR";


for capture in "$@"
do

    IN_DIR=$base_path/$capture
    IN_DIRFILES=$IN_DIR/*
    samp_ids_str=""

    OUT_RESULTS=/group/dubcovskygrp/projects/sequence/processed/mapped/results/genome-CS-AB-ref/$capture
    TRIMED_IN_DIR=/group/dubcovskygrp/projects/sequence/processed/trimmed-reads
    OUT_STATS=/group/dubcovskygrp/projects/sequence/processed/mapped/stats/genome-CS-ABU-ref/$capture
    ##Rmdup Configuration
    IN_DIR_RMDUP=/group/dubcovskygrp3/projects/sequence/processed/mapped/results/genome-CS-ABU-ref/$capture
    OUT_RESULTS_RMDUP=/group/dubcovskygrp3/projects/sequence/processed/mapped/rmdup/genome-CS-ABU-ref/$capture

    total_reads=$(cat $IN_DIR/*/*.adapter.stat | grep "total_reads" | cut -d' ' -f2 | awk '{total += $1}END{ print total}')
    echo "Total Reads in $capture: $total_reads"


    if (( $total_reads >= 180000000 )); then

        echo $IN_DIR
        if [ ! -e "$IN_DIR/done.flag" ]; then

            echo "Processing $capture in directory $IN_DIR"
            # time trimming process
            T="$(date +%s)"

            ##Find the adapters file
            ADAPTERS=$ADAPTORS_DIR/Bioo_adaptors.fasta
            echo "Using adaptor file: $ADAPTERS"

            for FILE in $IN_DIRFILES
            do

                ##For each folder in that sub-directory
                if [ -d $FILE  ]; then
                    basename=$(basename $FILE)
                    SAMPLE_NUMBER=$(echo $basename | awk -F_ '{print $NF}')
                    SAMPLE_NAME="Kronos$SAMPLE_NUMBER"
                    sampleoutdir="$OUT_BASEDIR/$SAMPLE_NAME/"
                    samp_ids_str="${samp_ids_str}${SAMPLE_NUMBER} " ##used later for looping through samples in a particular capture
                    echo "Sample Input directory: $FILE"
                    echo "SampleName: $SAMPLE_NAME"
                    echo "Output Directory: $sampleoutdir"

                    ##Find input files
                    IN1=$FILE/*_1.fq.gz
                    IN2=$FILE/*_2.fq.gz
                    echo "Input PE1 for sample $SAMPLE_NAME: $IN1"
                    echo "Input PE2 for sample $SAMPLE_NAME: $IN2"
                    zcat1="zcat $IN1 > $TMP_DIR/${SAMPLE_NAME}_raw_1.fq &"
                    echo $zcat1
                    eval $zcat1
                    zcat2="zcat $IN2 > $TMP_DIR/${SAMPLE_NAME}_raw_2.fq &"

                    echo $zcat2
                    eval $zcat2
                    echo "Waiting for zcat to finish"
                    wait

                    echo "zcat finished"
                    IN1=$TMP_DIR/${SAMPLE_NAME}_raw_1.fq
                    IN2=$TMP_DIR/${SAMPLE_NAME}_raw_2.fq

                    ##Make sure the output directory is made
                    if [ ! -d $sampleoutdir  ]; then
                        mkdir $sampleoutdir
                        echo "Created $sampleoutdir"
                    fi

                    ###Use Sickle for quality trimming and Scythe for adapter trimming
                    echo
                    trimcommand="time sickle pe -t illumina -q $QUAL_THRESH \
                    -f <(seqqs -e -p ${sampleoutdir}raw_${SAMPLE_NAME}_R1 $IN1 | scythe -a $ADAPTERS -p $PRIOR - 2> ${sampleoutdir}${SAMPLE_NAME}_R1_scythe.stderr) \
                    -r <(seqqs -e -p ${sampleoutdir}raw_${SAMPLE_NAME}_R2 $IN2 | scythe -a $ADAPTERS -p $PRIOR - 2> ${sampleoutdir}${SAMPLE_NAME}_R2_scythe.stderr) \
                    -o >(seqqs -e -p ${sampleoutdir}trimmed_${SAMPLE_NAME}_R1 - | gzip > ${sampleoutdir}${SAMPLE_NAME}_R1_trimmed.fq.gz) \
                    -p >(seqqs -e -p ${sampleoutdir}trimmed_${SAMPLE_NAME}_R2 - | gzip > ${sampleoutdir}${SAMPLE_NAME}_R2_trimmed.fq.gz) \
                    -s >(seqqs -e -p ${sampleoutdir}trimmed_${SAMPLE_NAME}_singles - | gzip > ${sampleoutdir}${SAMPLE_NAME}_singles_trimmed.fq.gz) 2> ${sampleoutdir}${SAMPLE_NAME}_sickle.stderr"

                    ##Run trim command
                    echo $trimcommand
                    eval $trimcommand
                    echo
                    echo "touch $FILE/done.flag"
                    touch $FILE/done.flag

                fi
            done

        fi
    else
        echo "Stopping script because $capture did not pass total_reads filter criteria. $total_reads < 180,000,000 reads"
        exit
    fi


    ####FINISHED TRIMMING READS -> NOW START MAPPING
    ###BWA Map Section
    if [ ! -e $OUT_RESULTS ]; then
        echo "Making $OUT_RESULTS"
        mkdir -p $OUT_RESULTS
    fi

    if [ ! -e $OUT_STATS ]; then
        echo "Making $OUT_STATS"
        mkdir -p $OUT_STATS
    fi

    if [ ! -e $OUT_RESULTS_RMDUP ]; then
        echo "Making $OUT_RESULTS_RMDUP"
        mkdir -p $OUT_RESULTS_RMDUP
    fi

    echo "[map-bwa.sh] $today"

    ###Fix samp_ids variable
    samp_ids=$(echo $samp_ids_str | sed 's/ *$//')

    for sample in $samp_ids
    do

        FW_READS=$TRIMED_IN_DIR/Kronos$sample/"Kronos"$sample"_R1_trimmed.fq.gz"
        REV_READS=$TRIMED_IN_DIR/Kronos$sample/"Kronos"$sample"_R2_trimmed.fq.gz"

        echo "[map-bwa.sh] temporarily uncompressing capture $capture sample Kronos$sample"

        zcat1="zcat $FW_READS > $TMP_DIR/Kronos${sample}_R1_trimmed.fq &"
        echo "[map-bwa.sh]: running $zcat1"
        eval $zcat1
        echoerr $zcat1
    done
    wait

    for sample in $samp_ids
    do

        FW_READS=$TRIMED_IN_DIR/Kronos$sample/"Kronos"$sample"_R1_trimmed.fq.gz"
        REV_READS=$TRIMED_IN_DIR/Kronos$sample/"Kronos"$sample"_R2_trimmed.fq.gz"

        echo "[map-bwa.sh] temporarily uncompressing capture $capture sample Kronos$sample"

        zcat2="zcat $REV_READS > $TMP_DIR/Kronos${sample}_R2_trimmed.fq &"
        echo "[map-bwa.sh]: running $zcat2"
        echoerr $zcat2
        eval $zcat2
    done
    wait

    for sample in $samp_ids
    do

        FW_READS=$TMP_DIR/"Kronos"$sample"_R1_trimmed.fq"
        REV_READS=$TMP_DIR/"Kronos"$sample"_R2_trimmed.fq"

        echo "[map-bwa.sh] processing capture $capture sample Kronos$sample"
        prefix="Kronos"$sample"_"$capture"_bwa_ABUgenome_"$today

        bwaalncmd1="time bwa aln -t 2 -I $REF $FW_READS > $OUT_RESULTS/$prefix.1.sai 2> $OUT_STATS/$prefix.1.log &"
        echo "[map-bwa.sh]: running $bwaalncmd1"
        echoerr $bwaalncmd1
        eval $bwaalncmd1

    done
    wait

    for sample in $samp_ids
    do

        FW_READS=$TMP_DIR/"Kronos"$sample"_R1_trimmed.fq"
        REV_READS=$TMP_DIR/"Kronos"$sample"_R2_trimmed.fq"

        echo "[map-bwa.sh] processing capture $capture sample Kronos$sample"
        prefix="Kronos"$sample"_"$capture"_bwa_ABUgenome_"$today

        bwaalncmd2="time bwa aln -t 2 -I $REF $REV_READS > $OUT_RESULTS/$prefix.2.sai 2> $OUT_STATS/$prefix.2.log &"
        echo "[map-bwa.sh]: running $bwaalncmd2"
        echoerr $bwaalncmd2
        eval $bwaalncmd2

    done
    wait

    for sample in $samp_ids
    do

        FW_READS=$TMP_DIR/"Kronos"$sample"_R1_trimmed.fq"
        REV_READS=$TMP_DIR/"Kronos"$sample"_R2_trimmed.fq"

        PU=$(head -1 $FW_READS| cut -f 1 -d ":" | sed s/@//)
        #prefix="Kronos"$sample"_"$capture"_bwa_ABgenome_"$today
        SAIFILE=($OUT_RESULTS/Kronos${sample}_${capture}_bwa_ABUgenome_*.1.sai)
        prefix=$(basename $SAIFILE .1.sai)
        echo "SAIFILE is $SAIFILE; Prefix is $prefix"

        bwasampe="time bwa sampe -N 10 -n 10 -r '@RG\tID:$sample\tPU:$PU\tLB:$capture\tSM:Kronos$sample' $REF $OUT_RESULTS/$prefix.1.sai $OUT_RESULTS/$prefix.2.sai $FW_READS $REV_READS > $OUT_RESULTS/$prefix.sam 2> $OUT_STATS/$prefix.log &"
        echo "[map-bwa.sh]: running $bwasampe"
        echoerr $bwasampe
        eval $bwasampe

    done
    wait

    for sample in $samp_ids
    do

        SAMFILE=($OUT_RESULTS/Kronos$sample*.sam)
        prefix=$(basename $SAMFILE .sam)

        #Redirection seems to cause an truncated error; must use the -o argument
        samtoolscmd="samtools view -S -b $OUT_RESULTS/$prefix.sam -o $OUT_RESULTS/$prefix.bam &"
        echo $samtoolscmd
        echoerr $samtoolscmd
        eval $samtoolscmd

    done
    wait

    for sample in $samp_ids
    do

        SAMFILE=($OUT_RESULTS/Kronos$sample*.sam)
        prefix=$(basename $SAMFILE .sam)

        samtoolscmd="samtools sort $OUT_RESULTS/$prefix.bam $OUT_RESULTS/$prefix.sorted &"
        echo $samtoolscmd
        echoerr $samtoolscmd
        eval $samtoolscmd
        echo "[map-bwa.sh]: done sorting bam $prefix"
    done
    wait

    #####START RMDUP
    ####Remove duplicates
    for sample in $samp_ids
    do

        BAMFILE=($IN_DIR_RMDUP/Kronos$sample*sorted.bam)
        prefix=$(basename $BAMFILE .bam)

        echo "[map-bwa.sh] processing capture $capture sample Kronos$sample prefix $prefix"
        echo "Starting PicardTools MarkDuplicates for $prefix"
        rmdupcmd="time /share/apps/picardtools-1.107/MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true I=$IN_DIR_RMDUP/$prefix.bam O=$OUT_RESULTS_RMDUP/$prefix.rmdup.bam M=$OUT_RESULTS_RMDUP/$prefix.rmdup.txt TMP_DIR=$TMP_DIR &"
        echo $rmdupcmd
        eval $rmdupcmd

        echo "[rmdup.sh]: done rmdup $prefix"

    done
    wait

    ###calculate stats
    for sample in $samp_ids
    do

        BAMFILE=($IN_DIR_RMDUP/Kronos$sample*sorted.bam)
        prefix=$(basename $BAMFILE .bam)

        perlcmd="time perl $SCRIPTS_DIR/K-sam_stats-v2.1.pl -f $REF -b $OUT_RESULTS_RMDUP/$prefix.rmdup.bam -o $OUT_RESULTS_RMDUP/$prefix.rmdup.samstats > $OUT_RESULTS_RMDUP/$prefix.rmdup.samstats.log &"
        echo $perlcmd
        eval $perlcmd
        echo "[rmdup.sh]: done samstats $prefix"

    done
    wait


done

##remove temp directory
rm -rf $TMP_DIR
