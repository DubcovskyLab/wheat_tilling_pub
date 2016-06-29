#!/bin/bash
############################
# 39.EMS_Preference_slurm.sh

# Paul Bailey	29.1.2016
############################

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



# Cadenza main data set (RUN TIME: 30 minutes, seems to need > 18 mem):
#source python-2.7.3
#python ~baileyp/bin/bases-bordering-mutation-by-type-2.py \
#-f cadenza_1200_unique.maps-part2_-d5_-s3_all_muts_no_low_pc_EMS_noWTs_sort_g.correctedNo_RH.maps.nobins.noindels.tsv \
#-r /tgac/workarea/group-tg/baileyp/39.ExonCapture_PlateH2N/mutation_context_analysis/IWGSC_v2_ChrU_3B_Ids_tidied.fa \
#-o ems-preference_main_ABD.out   


#source python-2.7.3
#python ~baileyp/bin/bases-bordering-mutation-by-type-2.py \
#-f Cadenza_repeat_mutations_noRH_3individuals_in_MAPS_format.txt \
#-r /tgac/workarea/group-tg/baileyp/39.ExonCapture_PlateH2N/mutation_context_analysis/IWGSC_v2_ChrU_3B_Ids_tidied.fa \
#-o ems-preference_repeat_muts_3individuals.out


# Cadenza mm data set (RUN TIME: 20 minutes, seems to need > 18 GB mem): :
#source python-2.7.3
#python ~baileyp/bin/bases-bordering-mutation-by-type-2.py \
#-f /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/ty_mm_VEP_analysis_-d3+_-s2+/TABCDEFGHIJKLMN_poorA-G_poorH-N_-d5+_-s3+_no_low_pc_EMS_noWTs_no_delns.tab \
#-r /tgac/workarea/group-tg/baileyp/39.ExonCapture_PlateH2N/mutation_context_analysis/IWGSC_v2_ChrU_3B_Ids_tidied.fa \
#-o ems-preference_mm.out


# Cadenza - looking at A, B and D genomes separately (RUN TIME: 30 minutes, seems to need > 18 mem):
#source python-2.7.3
#python ~baileyp/bin/bases-bordering-mutation-by-type-2.py \
#-f maps_output_D_genome \
#-r /tgac/workarea/group-tg/baileyp/39.ExonCapture_PlateH2N/mutation_context_analysis/IWGSC_v2_ChrU_3B_Ids_tidied.fa \
#-o ems-preference_main_D_genome.out


# Kronos main data set (RUN TIME: 20 minutes, seems to need > 13GB mem :
source python-2.7.3
python ~baileyp/bin/bases-bordering-mutation-by-type-2.py \
-f Kronos_mutations_unique_only_noRH_in_MAPS_format.tab \
-r /tgac/workarea/collaborators/dubcovskylab/reference/IWGSC_CSS_AB-TGAC_UCW_v1.fa \
-o ems-preference_unique_muts_all.out