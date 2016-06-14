#!/bin/bash -l
mkdir -p logs

date
module load perl
hostname

export input=`cat /group/dubcovskygrp4/mmreads_test/alt-locs/MAPS-calls/individual_list.txt | sed -n "$SLURM_ARRAY_TASK_ID p" | awk '{print $1}' `

bamdir="/home/trhowell/TILLING/rescued_reads/rescue-I-final"
outdir="/group/dubcovskygrp4/mmreads_test/alt-locs"
reference="/group/dubcovskygrp/projects/reference/IWGSC_CSS_AB-TGAC_UCW_v1.fa"

mkdir -p ${outdir}

for i in {3..6}; do
    for j in {2..4}; do
        grep -P "${input}\t" /group/dubcovskygrp4/mmreads_test/alt-locs/MAPS-calls/Kronos.mapspart2.HetMinCov${i}HetMinPer15HomMinCov${j}-All-mm.selfdedup.no3690no3746.tsv > ${outdir}/${input}-HetMinCov${i}HomMinCov${j}.muts &
    done
done
wait

for i in {3..6}; do
    for j in {2..4}; do
        samtools view ${bamdir}/${input}_*.bam | /home/trhowell/bin/39.getContigSets.pl \
        -fa ${reference} \
        -mi ${outdir}/${input}-HetMinCov${i}HomMinCov${j}.muts \
        -mo ${outdir}/${input}-HetMinCov${cov}HomMinCov${j}.muts_plus_alt_locs.txt \
        -s ${outdir}/${input}-HetMinCov${cov}HomMinCov${j}.alt-locs.txt
    done
done

# Log ending time
date

