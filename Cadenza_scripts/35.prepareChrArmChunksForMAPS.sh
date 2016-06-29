#!/bin/bash

##################################
# 35.prepareChrArmChunksForMAPS.sh
#
# Paul Bailey	5.5.2015
#
##################################

referenceFile=/tgac/workarea/group-cg/baileyp/IWGSC_v2_ChrU_Ref/IWGSC_v2_ChrU.fa					# The full and final reference for the project (IWGSC CSS minus CSS 3B + new 3BSeq + CSS 3B not in 3BSeq + ChrU)

#referenceFile=/tgac/workarea/collaborators/dubcovskylab/reference/IWGSC_CSS_AB-TGAC_UCW_v1.fa		# Reference for Ksenia and Hans' Kronos samples

#referenceFile=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_masked_ref/39.ExonCapture_masked_ref_10_CadWT_samples_v1/masked_contigs.fasta	# Masked reference for Cadenza samples

#referenceFile=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_masked_ref/masking_for_Kronos_samples_merged10sam/masked_contigs.fasta			# Masked reference for Kronos samples 

# Prepare contigs for each chr arm from reference file:
chrArms=(1AL 1AS 1BL 1BS 1DL 1DS 2AL 2AS 2BL 2BS 2DL 2DS 3AL 3AS 3DL 3DS 4AL 4AS 4BL 4BS 4DL 4DS 5AL 5AS 5BL 5BS 5DL 5DS 6AL 6AS 6BL 6BS 6DL 6DS 7AL 7AS 7BL 7BS 7DL 7DS)	# minus 3B and Cadenza which have to be done separately
#chrArms=(1AL 1AS 1BL 1BS 2AL 2AS 2BL 2BS 3AL 3AS 3B 4AL 4AS 4BL 4BS 5AL 5AS 5BL 5BS 6AL 6AS 6BL 6BS 7AL 7AS 7BL 7BS) 	# no D arms
# Also ...grep "ChrU"... separately as a special case:
#chrArms=(TGAC_Cadenza_U)
chrArms=(UCW_Kronos_U)
# ... and 3B (IWGSC_3BSEQ_3B AND IWGSC_CSS_3B) - NB - unlocalised new3BSEQ fasta header lines have effectively 2 ids - so you need to remove it:  more 3B_ids_before_rming_v442_field | awk '{print $1}' >3B_ids   
#chrArms=3B


for chrArm in ${chrArms[@]}; do

	echo $chrArm

	bsub -J $chrArm -oo ${chrArm}_get_ids_fetch_fasta_subset.log "
																   # NB - added masked substitution for the masked ref and needed to eliminate a space char for the non-masked contigIds (this space char is now fixed in the masking script)
	cat $referenceFile | grep \">IWGSC_CSS_$chrArm\" | sed 's/>//' | sed 's/ masked//' | sed 's/ //' > ${chrArm}_ids
	#cat $referenceFile | grep TGAC_Cadenza_U | sed 's/>//' | sed 's/ masked//' | sed 's/ //' > ${chrArm}_ids
	#cat $referenceFile | grep UCW_Kronos_U | sed 's/>//' > ${chrArm}_ids
	#cat $referenceFile | grep 'IWGSC_3BSEQ_3B\|IWGSC_CSS_3B' | sed 's/>//' | sed 's/ masked//' > ${chrArm}_ids
	36.fetch_fasta_subset.pl $referenceFile ${chrArm}_ids > ${chrArm}.fa
	# Get length of each contig and prepare bedfile (start coords are zero-based):
	source exonerate-2.2.0
	fastalength ${chrArm}.fa | awk '{print \$2 \"\t\" 0 \"\t\" \$1}' > ${chrArm}.fa.bed
	"
done
