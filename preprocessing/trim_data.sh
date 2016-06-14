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

##Log Start time
date
hostname
module load sickle
module load scythe
module load seqqs 

if [ ! -e $TMP_DIR ]; then
    echo "Making $TMP_DIR"
    mkdir $TMP_DIR
fi
echo "created and assigned TMP_DIR to $TMP_DIR";


for capture in "$@"
do

    IN_DIR=$base_path/$capture
    IN_DIRFILES=$IN_DIR/*

    OUT_RESULTS=/group/dubcovskygrp/projects/sequence/processed/mapped/results/genome-CS-AB-ref/$capture

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

done

##remove temp directory
rm -rf $TMP_DIR
