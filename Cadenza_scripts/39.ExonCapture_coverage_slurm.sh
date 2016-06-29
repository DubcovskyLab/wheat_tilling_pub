#!/bin/bash
##################################
# 39.ExonCapture_coverage_slurm.sh

# Paul Bailey	18.3.2016
##################################

echo "#################"
echo "SCRIPT VARIABLES:"
echo "\$SLURM_JOBID	" $SLURM_JOBID
echo "Number of cpu's requested	" $cpu
echo "Amount of memory requested (MB)	" $mem_mb
echo "\$SLURM_JOB_NAME	" $SLURM_JOB_NAME
echo "\$SLURM_SUBMIT_DIR	" $SLURM_SUBMIT_DIR
echo "List of nodes allocated to the job	" $SLURM_JOB_NODELIST
# Also useful would be: start time of job
# Also useful would be: end time of job
# Also useful would be: total length of time of job
echo "#################"


# Bedfile for calculating the coverage stats:
#bedfilePlusPadding=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/sequence_space_coverage/blast/blast_exon_queries_vs_gcontigs_db_all_hits/exons_0-based_200bp_padding_sort_-k11_-k22n_merge.bed 
#bedfileMinusPadding=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/sequence_space_coverage/blast/blast_exon_queries_vs_gcontigs_db_all_hits/exons_0-based_sort_-k11_-k22n_merge.bed
# The "tophit" bedfile for calculating the stats below: 
bedfilePlusPadding=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/sequence_space_coverage/blast/exons_0-based_200bp_padding_sort_-k11_-k22n_bt_merge.bed
bedfileMinusPadding=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/sequence_space_coverage/blast/exons_0-based_sort_-k11_-k22n_merge.bed


:<<***string_marker_for_commenting_out_code***


# Count the number of reads with map quality >= 20:
$source_production_env; source samtools-0.1.18
samtools view -c -F4 -q 20 *_sampe_sort_markdup_rm.bam > sampe_sort_markdup_rm_st_-F4_-q20_count.log


# Count % reads on-target:
# Samtools:
$source_production_env; source samtools-0.1.18
# Plus padding:
samtools view -c -L $bedfilePlusPadding *sampe_sort_markdup_rm.bam  > sampe_sort_markdup_rm_st_-L_200padding_count.log
# Minus padding:
samtools view -c -L $bedfileMinusPadding *sampe_sort_markdup_rm.bam  > sampe_sort_markdup_rm_st_-L_0padding_count.log

#Bedtools:
$source_production_env; source bedtools-2.17.0
# Plus padding:
bedtools intersect -abam *sampe_sort_markdup_rm.bam \
-b $bedfilePlusPadding \
-f 0.000000001 \
-u | wc -l > sampe_sort_markdup_rm_bt_intersect_200padding_count.log
### NB - THIS METHOD IS NOT WORKING NOW - LAST TIME I SAVEWD TO A FILE - OTHERWISE I CAN'T SEE WHAT'S GOING ON
# -f Minimum overlap required as a fraction of A. Default is 1E-9 OR 0.000000001 (i.e. 1bp).
# NB -wa doesn't have to be specified to get the reads to be output.
# Minus padding:
bedtools intersect -abam *sampe_sort_markdup_rm.bam \
-b $bedfileMinusPadding \
-f 0.000000001 \
-u | wc -l > sampe_sort_markdup_rm_bt_intersect_0padding_count.log


***string_marker_for_commenting_out_code***


# Get average on-target coverage and base coverage >=1, >=6 etc for current sample (RUNTIME: up to 2h; 12GB mem):
$source_production_env; source bedtools-2.17.0
# Plus padding:
coverageBed -abam *sampe_sort_markdup_rm.bam \
-b $bedfilePlusPadding \
-d \
| ~baileyp/bin/35.getCoverageStats_On-targetRegions.pl \
> sampe_sort_markdup_rm_covBed_-d_av_cov_more_1_6_10_15_20_25_100x_covs_200padding.log
# Minus padding:
coverageBed -abam *sampe_sort_markdup_rm.bam \
-b $bedfileMinusPadding \
-d \
| ~baileyp/bin/35.getCoverageStats_On-targetRegions.pl \
> sampe_sort_markdup_rm_covBed_-d_av_cov_more_1_6_10_15_20_25_100x_covs_0padding.log



# Get average coverage accross whole genome, off-target and on-target regions (RUNTIME up to 2 days!)
#bedtools genomecov -d -ibam *sampe_sort_markdup_rm.bam \
#-g $bedT_GenomeFile \
#| ~baileyp/bin/35.getCoverageStats_WholeGenome.pl $numbrBasesInIWGSC_v1 \
#> sampe_sort_markdup_rm_bt_genomecov_-d_more_1_6_10_20x_covs.log



# Run Bedtools coverage -hist. The hist option gives you the coverage for each feature in the gff (in this case I wanted to calculate average coverage per Ensembl gene)
# (although this is split up by exact coverage on different rows - see my Bedtools notes).
# RUNTIME: 2.5 h, 15 GB memory:
# For IWGSC_v1 CSS contigs, require Ensembl release 23 gff file (CDS lines only):
#bedfile=/tgac/workarea/group-cg/baileyp/EnsemblPlantsRelease23/Triticum_aestivum.IWGSP1.23_CDS_lines_sort_-k11_-k44n.gff3
# For new 3B contigs, require Ensembl release 25 gff file (CDS lines only).
# Have merged these with the 23 gff file:
#bedfile=/tgac/workarea/group-cg/baileyp/EnsemblPlantsRelease23/Triticum_aestivum.IWGSP1.23_25_3B_CDS_lines_sort_-k11_-k44n_IWGSC_v1contigs_3Bcontigs_reformatted.gff3
#source production_env; source bedtools-2.17.0
#coverageBed -abam *sampe_sort_markdup_rm.bam -b $bedfile -hist \
#> sampe_sort_markdup_rm_bt_cov_-hist_plus_3B.out