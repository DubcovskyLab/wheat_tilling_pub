#!/bin/bash

##########################################
# 35.runMAPS_Pipeline_GetSNP_Annotation.sh
#
# Paul Bailey	25.1.2015
#
##########################################

# Remove samples from the maps-part2 Excel results summary file (-d1 --> -d 12) that have these issues:
# (need to load results into Excel first and examine the data)
# 1. low % EMS
# 2. WT's
# 3. samples with normal # hets but low # homs need to be removed from the main set, 
#    just to calculate a better average het/hom ratio


# 1. Remove low % EMS samples and place them at the end of the data:
testSamples54='LIB8438'
plateH='LIB15597'
plateI='LIB15312'
plateN='LIB15385'
:<<***string_marker_for_commenting_out_code***
cat ExonCapture_maps-part2_-dall_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel.txt | \
grep "$testSamples54\|$plateH\|$plateI\|$plateN" \
> ExonCapture_maps-part2_-dall_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel_low_pc_EMS.txt
***string_marker_for_commenting_out_code***

# 2. Remove WT samples and place them at the end of the data:
testSamples54_WT='LIB5782\|LIB8440'
CadWTs='LIB10700\|LIB1070[123456]\|LIB1002[012]\|LIB10495\|LIB1031[345678]\|LIB1100[789]\|LIB1149[012345]'
:<<***string_marker_for_commenting_out_code***
cat ExonCapture_maps-part2_-dall_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel.txt | \
grep "$testSamples54_WT\|$CadWTs" \
> ExonCapture_maps-part2_-dall_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel_WTs_only.txt
***string_marker_for_commenting_out_code***

# 3. Samples with low # homs: 
testSamples54_low_homs='LIB8411'
plateA='LIB10633\|LIB10698'
plateB='LIB9954\|LIB9958'
plateC='LIB11218\|LIB11227\|LIB11178'
plateD='LIB10410\|LIB10436\|LIB10442\|LIB10469'
plateE='LIB10219\|LIB10304'
plateF='LIB10995\|LIB10914\|LIB10937\|LIB10946\|LIB10951'			# 28.9.2015 - NB I had missed LIB10937 because it was a poor capture so not yet done 
plateG='LIB11484\|LIB11485\|LIB11486\|LIB11429\|LIB11410\|LIB11425'	# 28.9.2015 - NB I had missed LIB11484 because it was a poor capture so not yet done
plateJ='LIB16150'
plateK='LIB15880\|LIB15953\|LIB15973'
plateL='LIB16260\|LIB16208'
plateM='LIB16550\|LIB16562'
plateN_low_homs='LIB15443'
plate20='LIB17362'
# First grep these samples and add them to the bottom of the data:
:<<***string_marker_for_commenting_out_code***
cat ExonCapture_maps-part2_-dall_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel.txt | \
grep "$testSamples54_low_homs\|$plateA\|$plateB\|$plateC\|$plateD\|$plateE\|$plateF\|$plateG\|$plateJ\|$plateK\|$plateL\|$plateM\|$plateN_low_homs\|$plate20" \
> ExonCapture_maps-part2_-dall_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel_low_homs_only.txt
***string_marker_for_commenting_out_code***

# Then we want to remove the low % EMS and WT samples plus the low#hom samples from the table so we get an accurate het/hom ratio and SNP number:
:<<***string_marker_for_commenting_out_code***
cat ExonCapture_maps-part2_-dall_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel.txt | \
grep -v "$testSamples54\|$plateH\|$plateI\|$plateN\|$testSamples54_WT\|$CadWTs\|$testSamples54_low_homs\|$plateA\|$plateB\|$plateC\|$plateD\|$plateE\|$plateF\|$plateG\|$plateJ\|$plateK\|$plateL\|$plateM\|$plateN_low_homs\|$plate20" \
> ExonCapture_maps-part2_-dall_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel_no_low_pc_EMS_noWTs_no_low_homs.txt 
***string_marker_for_commenting_out_code***


# The file prepared here is for the VEP analysis and getting the stats:
#Now need to repeat this 'grep -v' on  maps-part2_3x_4x_5x_more6x_lines.txt - in order to remove all low % EMS and WT samples before the VEP analysis:
:<<***string_marker_for_commenting_out_code***
#cat maps-part2_3x_4x_5x_more6x_lines.txt | \
#grep -v "$testSamples54\|$plateH\|$plateI\|$plateN\|$testSamples54_WT\|$CadWTs" \
#> maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs.txt
# Same for other data:
cat maps-part2_-d6_-s4_all_muts.tab | \
grep -v "$testSamples54\|$plateH\|$plateI\|$plateN\|$testSamples54_WT\|$CadWTs" \
> maps-part2_-d6_-s4_all_muts_no_low_pc_EMS_noWTs.tab
# NBNB - Can now get the abolute mutation calling statistics with this file (or the sorted file below) and for the VEP analysis. 
***string_marker_for_commenting_out_code***

# Notes:
# 1. The above data is coming from the individual files for -d2 --> -d 12 and not the file maps-part2_3x_4x_5x_more6x_lines.txt
# 	 This means that the values in the Excel table for -d6,-s3 will be slightly different for the data that is processed for the web pages.
# 2. I didn't remove the indel mutations before counting for the Excel file - but there aren't that many of these and they are still mutations!
# 3. Two samples completely failed from PlateN, LIB15438 and LIB15471 - not enough reads
# 4. Five samples have worked but I forgot to add them to any MAPs run: 
#		PlateL LIB16261	(should have gone into a normal MAPS run)
#		PlateD LIB10467 LIB10468 LIB10469 and LIB10470 (should have gone into one of the poor capture MAPS runs)
#			LIB10469 looks to be a low ~hom line from an early MAPS run. 





# This file - maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs.txt
# was used to calculate the % EMS at -d2 to -d9 and -s3 to -s9:
#Commands used:- e.g.:
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs.txt| awk '($8 ~ /het/ && $10 ==3)' | wc -l 
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs.txt| awk '($8 ~ /het/ && $10 ==3)' | awk '$11 ~ /GA|CT/' | wc -l
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs.txt| awk '($8 ~ /hom/ && $10 ==3)' | wc -l
# To get the values for hetMinCov at 2x coverage, used this file: maps-part2_-lX_-c10_-d2_-p10_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt





# Now this file -  maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs.txt - needs to be sorted here:
# 1. to bring all the chr arms together from the different MAPS runs 
# 2. to see identical SNPs between each line 
# 3. to check that the VEP output is one line per effect
# This file has no header line - OK.
# Also removing blank lines that appear in this file (the tail command when catted the MAPS output leaves a gap at the end of each file section! 29.9.2015 - I think these lines now get removed in the last script - so the grep here is now obsolete)
#sort -g maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs.txt | grep '^TGAC_Cadenza_U_\|^IWGSC_' >  maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt
# For homMinCov2 (NB - found that the grep was still required!):
#sort -g  maps-part2_-d6_-s4_all_muts_no_low_pc_EMS_noWTs.tab | grep '^TGAC_Cadenza_U_\|^IWGSC_' >  maps-part2_-d6_-s4_all_muts_no_low_pc_EMS_noWTs_sort_g.tab







# Also prepare 2 MAPS outputs, one for CHR 1-2; 4-7, TGAC_U and one for chr 3B then process them separately through VEP with different annotation files:
#MAPS_Infile=combined.mapspart2.HetMinCov4HetMinPer15HomMinCov3.corrected.deduped.10kb_bins.RH.byContig.MI.No_RH.maps.nokronos0.noindels.nobins_het4_ONLY.tsv
#cat $MAPS_Infile \
#| grep -v '^IWGSC_3BSEQ_3B' > ${MAPS_Infile}_no_new3B.txt

# To extract new 3B, searching for IWGSC_3BSEQ_3B should get the pseudomolecule and all the unordered contigs.
# Then need to convert the Ids to Ensembl format with sed: 
#cat $MAPS_Infile \
#| grep '^IWGSC_3BSEQ_3B' | sed 's/IWGSC_3BSEQ_3B_traes3bPseudomoleculeV1/3B/' | sed s/IWGSC_3BSEQ_3B_// \
#> ${MAPS_Infile}_new3B.txt
# The above outputs are the the files that need to be converted to vcf.


# Now create a unique list of old IWGSC contig Ids from the "*no_new3B.txt" file
# These ids are used in the Ensembl converter script (AssemblyMapper_NA.pl from Nikolai)) to get the ids on the whole pseudochromosome.
# Info on coordinates :
# coord_system:version:seq_region_name:start:end:strand
#cat ${MAPS_Infile} \
#| grep -v _3B_ | awk '{print $1 "\t" $2}' \
#| sort -u \
#| awk '{print "scaffold:IWGSC2:"$1 ":" $2 ":" $2 ":1"}' \
#> iwgsc_v1_for_convertion_to_iwgsc_v2_coords.txt
# NB - remove for preparing Kronos file: grep -v _3B_ 


# Convert IWGSC_v1 contig ids to IWGSC_v2 here (RUN TIME: 11 h for PlatesTABCDEFGHI; ~3,500 records per second)
# Use Nikolai's modified script, AssemblyMapper_NA.pl
# Had to modify script with a 'lib' to locate the module.
# Copy iwgsc_v1_for_convertion_to_iwgsc_v2_coords.txt into ~baileyp/scratch to run this command because /tgac/workarea/ not visible from software node - e.g.:
#cp -p iwgsc_v1_for_convertion_to_iwgsc_v2_coords.txt ~baileyp/scratch/VEP_Convert_IWGSC_v1_ids_to_IWGSC_v2__15Oct2015
### 18.8.2015 - could also do this for all bases in the whole ref! - I think this would take a much longer time!
:<<***string_marker_for_commenting_out_code***
ssh software
[ need to cd again into ~baileyp/scratch ]
source perl-5.16.2; /usr/users/ga002/baileyp/ProgramFiles/EnsemblAssemblyMapperFromNikolai/AssemblyMapper_NA.pl \
-s triticum_aestivum \
-f iwgsc_v1_for_convertion_to_iwgsc_v2_coords.txt \
--cs_version IWGSC2 \
--cs_name chromosome \
| grep '^scaffold' \
> iwgsc_v1_converted_to_iwgsc_v2_coords.txt 2>iwgsc_v1_converted_to_iwgsc_v2_coords.err &
***string_marker_for_commenting_out_code***
# Output format example:
#	# scaffold:IWGSC2:IWGSC_CSS_1AL_scaff_1000404:4355:4355:1
#	scaffold:IWGSC2:IWGSC_CSS_1AL_scaff_1000404:4355:4355:1,chromosome:IWGSC2:1A:231871519:231871519:1	<-- need to pass this line to 35.prepareVCF_PlusFile.pl 
### next time add this to above cmd:
#cat iwgsc_v1_converted_to_iwgsc_v2_coords.txt | grep '^scaffold' > iwgsc_v1_converted_to_iwgsc_v2_coords_prepare.txt
# Copy back to group-cg area





# Prepare a vcf file from MAPS output for SnpEFF or VEP - vcf e.g.:
#CHROM POS     ID        REF    ALT     QUAL FILTER INFO
#IWGSC_CSS_1AL_scaff_84503       352     .       C       T       .       .       .
#IWGSC_CSS_1AL_scaff_126384      2039    .       G       A       .       .       .
#IWGSC_CSS_1AL_scaff_276808      2169    .       C       T       .       .       .
:<<***string_marker_for_commenting_out_code***
fileList=(maps_mutations_-d4_-s3_het4_ONLY.tab)

for file in ${fileList[@]}; do

	echo "Processing $file"

	# Eliminate astericks and other base calls in columns $5 and $6 for the WT and mut bases respectively e.g. '.+2CT' like this otherwise SnpEff crashes.
	# Leaving in the non-EMS changes here.
	# NBNB - the conditional SHOULD be && not || (!) because we only want to let through rows that have [GATC] in both column 5 and 6!!
	# NBNB - It is important to use column $5 of the MAPS output rather than column $3, then where the reference is wrong, MAPS can still identify a mutation in the cadenza WT.
	# This is because VEP comes back with no annotation for where the ref base is same as the alt allele because the cadenza WT is apparently different to the reference
	awk '$5 !~ /[\.\+\*\d]/ && $6 !~ /[\.\+\*\d]/ {print $1 "\t" $2 "\t" "." "\t" $5 "\t" $6 "\t" "." "\t" "." "\t" "."}' \
	$file \
	> ${file}.vcf

	# Also create the vcf_plus file here.
	# Input files to script:
	# 1. table of info for all samples: ExonCapture_PlateTABCDEFGHIJKLMN.txt [ merge of 35.ExonCapture_300_Samples/35.ExonCapture_PreparingTable.txt ExonCapture_PlateH2N.txt - contains both headers! Don't think it matters ]
	# 2. IWGSC v1 coordinates on the pseudomolecules
	# NB - Run this on the interactive node and check that I am pointing to the correct sample table in the script
	printf '%b' "scaffold\tchr\tlibrary\tline\tposition\tchromosome_position\tref_base\twt_base\talt_base\thet/hom\twt_cov\tmut_cov\tconfidence\thomoeologs_present(1,2or3)\n" > ${file}.vcf_plus
	cat $file \
	| awk '$5 !~ /[\.\+\*\d]/ && $6 !~ /[\.\+\*\d]/' \
	| 35.prepareVCF_PlusFile.pl \
	/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/VEP_analysis_aux_files/ExonCapture_PlateTABCDEFGHIJKLMN.tab \
	iwgsc_v1_converted_to_iwgsc_v2_coords.txt \
	>> ${file}.vcf_plus
	#File for HetMinCov3+, HomMinCov3+: /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/VEP_analysis_aux_files/iwgsc_v1_converted_to_iwgsc_v2_coords.txt \
	# /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/VEP_analysis_aux_files/iwgsc_v1_converted_to_iwgsc_v2_coords.txt \
	
	
	# Also creating a small file of the deletions and insertions:
	# Note that I'm printing out column $5 not $3 as above - this row contains the info on WT insertions:
	# NBNB - the conditional SHOULD be || here because we want to let through any row that has the * characters at either position.
#	awk '$5 ~ /\*/ || $6 ~ /\*/ {print $1 "\t" $2 "\t" "." "\t" $5 "\t" $6 "\t" "." "\t" "." "\t" "."}' \
#	$file \
#	> ${file}_indels_only.vcf

#	printf '%b' "scaffold\tchr\tlibrary\tline\tposition\tchromosome_position\tref_base\twt_base\talt_base\thet/hom\twt_cov\tmut_cov\tconfidence\thomoeologs_present(1,2or3)\n" > ${file}_indels_only.vcf_plus
#	cat $file | \
#	awk '$5 ~ /[\.\+\*\d]/ || $6 ~ /[\.\+\*\d]/' \
#	| 35.prepareVCF_PlusFile.pl \
#	/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/VEP_analysis_aux_files/ExonCapture_PlateTABCDEFGHIJKLMN.txt \
#	/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/VEP_analysis_aux_files/iwgsc_v1_converted_to_iwgsc_v2_coords.txt \
#	>> ${file}_indels_only.vcf_plus
	
	
	
	### 13.10.2014 - do I need these rows below now?
	
	# This just produces the larger indels (NB - just using the above file for annotating all indels together though):
#	#awk '$5 ~ /[\.\+\d]/ || $6 ~ /[\.\+\d]/ {print $1 "\t" $2 "\t" "." "\t" $5 "\t" $6 "\t" "." "\t" "." "\t" "."}' \
#	#${file} \
#	#> ${file}_multiple_insertions_only.vcf
	
	# Also printing out the full MAPS lines for the indel data.
	# Need to grab a MAPS output header from somewhere:
#	head -n 1 PlateA/MAPS_bams_to_use_LIB10604-LIB10701/1AL/maps-part2_-l22_-c10_-d6_-p10_-s3.txt > ${file}_indels_only_MAPS_output.tab
#	awk '$5 ~ /[\.\+\*\d]/ || $6 ~ /[\.\+\*\d]/ {print $0}' \
#	$file \
#	>> ${file}_indels_only_MAPS_output.tab

	# This just produces the larger indels:
#	head -n 1 PlateA/MAPS_bams_to_use_LIB10604-LIB10701/1AL/maps-part2_-l22_-c10_-d6_-p10_-s3.txt > ${file}_multiple_insertions_only_MAPS_output.tab
#	awk '$5 ~ /[\.\+\d]/ || $6 ~ /[\.\+\d]/ {print $0}' \
#	${file} \
#	>> ${file}_multiple_insertions_only_MAPS_output.tab
done
***string_marker_for_commenting_out_code***
# The vcf_plus file has a header line so the line count is one more than the .vcf file.
# NBNB - this file doesn't contain the indels so has less rows than the sort_g_[no]_new3B.txt files!

### possibly need to add the ^ char to the match string as this is also present in the base column - 10.6.2015 - can't find it now...





# VEP:
# Before running the VEP, sort unique the .vcf file to find the number of redundant SNPs:
#sort -u maps-part2_-lX_-c10_-d3_-p10_-s3_platesABCDEFG_some_lines_removed.txt.vcf | wc -l
# Must still use the .vcf and .vcf_plus file prepared above for the VEP and not the non-redundant set revealed by sort -u

# Run the VEP (RUNTIME: cpu=15; mem=1.4 GB; 5h40' for 8.2 million SNPs):
#:<<***string_marker_for_commenting_out_code***
cacheVersion=30			# 23 for no_new3B; 25 for new 3B; 30 for +sift
no_new3B=new3B		# no_new3B OR new3B
indels=					# '_indels_only' OR leave blank

infile=maps_mutations_-d4_-s3_new3B_het4_ONLY.tab.vcf

outfile=VEP_-d4_-s3_${no_new3B}_het4_ONLY.vcf
logfile=VEP_-d4_-s3_${no_new3B}_het4_ONLY.log
sift='--sift s'	# '--sift s'

cpu=15	# Was using 15

#bsub -J VEP_${no_new3B} -q Prod128 -n $cpu -R "rusage [mem=2000] span[hosts=1]" -oo $logfile " \
									  # -t 0-7:00 on the full set
sbatch -J VEP_${no_new3B} -p tgac-long -t 0-7:00 -c $cpu --mem=3000 -o $logfile -e ${logfile}_err --wrap="

# For use with sift in Release30 (NB location of the VEP anno file is still in ensembl-tools-release-78!) (RUNTIME: now seems to take 55 minutes for 6m muts):
# source perl-5.16.2; ~baileyp/ProgramFiles/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \

# For use with sift in Release23 and 25:
#source perl-5.16.2; ~baileyp/ProgramFiles/ensembl-tools-release-78/scripts/variant_effect_predictor/variant_effect_predictor.pl \

source perl-5.16.2; ~baileyp/ProgramFiles/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
--offline \
--cache \
--cache_version $cacheVersion \
--fork $cpu \
--dir_cache ~baileyp/ProgramFiles/ensembl-tools-release-78 \
-i $infile \
--vcf \
-o $outfile \
--force_overwrite \
--species triticum_aestivum \
$sift
"
#***string_marker_for_commenting_out_code***
#--cache_version 25				use the number in the required VEP dir in the triticum_aestivum directory: 23_IWGSC1, 25_IWGSC2
#--no_consequence             	can add an existing vcf file and it will add to it   - don't think this worked         
#--sift b                       try using sift - it might not work for wheat db 
#--fasta                        this may give you the sequence surrounding the SNP!!! 



# Prepare to paste VEF output with the vcf_plus file.
# Need to remove the first 3 header lines from the VEF output file:
#tail -n +4 deln_VEP_-d5_-s3_no_new3B.vcf > deln_VEP_-d5_-s3_no_new3B.vcf_paste_prepared.vcf 
#tail -n +4 VEP_-d3_-s2_new3B.vcf > VEP_-d3_-s2_new3B_paste_prepared.vcf
#tail -n +4 VEP_-s2_PlateK_5DL_6AS_no_new3B.vcf > VEP_-s2_PlateK_5DL_6AS_no_new3B_paste_prepared.vcf 
#tail -n +4 VEP_-s2_PlateK_5DL_6AS_new3B.vcf > VEP_-s2_PlateK_5DL_6AS_new3B_paste_prepared.vcf


# Also need to remove empty lines present in the vcf_plus file:
### 12.6.2015 - NB - I've now removed these lines at the very beginning!!!!!
#hetMinCov=3; grep -v 'No information available for this SNP!' maps-part2_-lX_-c10_-d3_-p10_-s3_platesABCDEFG_some_lines_removed_sort_g.txt.vcf_plus \
#>  maps-part2_-lX_-c10_-d3_-p10_-s3_platesABCDEFG_some_lines_removed_sort_g.txt.vcf_plus_paste_prepared
# NB - 20.3.2015 - this search string comes from my script that creates the vcf_plus file - still required (as of 12.6.2015, most probably not but do anyway)
### 20.8.2015 - yes i think i can miss this little step out now!
### Good still to check that the vcf_plus file and VEP*paste_prepared.vcf have the same number of lines! 

# Glue the lines of the .vcf_plus file and the VEP output using paste.
:<<***string_marker_for_commenting_out_code***
# For no_new3B:
#paste maps-part2_-d3_-s2_all_muts_no_low_pc_EMS_noWTs_sort_g.corrected.10kb_bins.RH.byContig.MI.No_RH.maps_binIds_rmed_het3_hom2.tsv_no_new3B.txt.vcf_plus \
#VEP_-d3_-s2_no_new3B.vcf_paste_prepared.vcf \
#> VEP_-d3_-s2_no_new3B.vcf_plus_info.txt
#paste maps-part2_all_hom==2_PlateK_5DS_6AL_no_new3B.txt.vcf_plus \
#VEP_-s2_PlateK_5DL_6AS_no_new3B_paste_prepared.vcf \
#> VEP_-s2_PlateK_5DL_6AS_no_new3B_plus_info.txt
# For new3B:
#paste maps-part2_-d3_-s2_all_muts_no_low_pc_EMS_noWTs_sort_g.corrected.10kb_bins.RH.byContig.MI.No_RH.maps_binIds_rmed_het3_hom2.tsv_new3B.txt.vcf_plus \
#VEP_-d3_-s2_new3B_paste_prepared.vcf \
#> VEP_-d3_-s2_new3B.vcf_plus_info.txt
#paste maps-part2_all_hom==2_PlateK_5DS_6AL_new3B.txt.vcf_plus \
#VEP_-s2_PlateK_5DL_6AS_new3B_paste_prepared.vcf \
#> VEP_-s2_PlateK_5DL_6AS_new3B_plus_info.txt
***string_marker_for_commenting_out_code***





# Now format the table in ideal format (for VEP):
:<<***string_marker_for_commenting_out_code***
#fileList=(VEP_-d3_-s3_new3B.vcf_plus_info.txt VEP_-d3_-s3_no_new3B.vcf_plus_info.txt)
fileList=(VEP_-s2_PlateK_5DL_6AS_no_new3B_plus_info.txt VEP_-s2_PlateK_5DL_6AS_new3B_plus_info.txt)

for file in ${fileList[@]}; do

	echo "Processing $file"

cat $file | perl -e '

	use strict;
	use warnings;
	while(my $line = <>)	{
		#print $line;
		
		my @fields =  split /\t/, $line;
		
		if($line =~ m/Contig/)	{# Build header line for the SnpEff output
			my @fieldsSubset = @fields[0..13];
			print join("\t", @fieldsSubset);
			print "\t\t\t\tgene\tfeature\tconsequence\tcDNA_position\tCDS_position\tamino_acid\tcodon\t";
			print "\n";
			next;
		}
		
		# 11.6.2015 - Now printing one row for each effect within the loop below.
		#my @fieldsSubset = @fields[0..14];				# Gets all fields from the vcf_plus file plus Contig + position for safety
		#print join("\t", @fieldsSubset);
		#print "\t";

		# Getting the info in the last column - SnpEff info:
		my @eachEffectSplit = split /,/, $fields[22];	# It should be always the zeroth field
		#print $semiColonFields[0];
		
		
		foreach my $eff (@eachEffectSplit)	{

			chomp $eff;
#			print $eff;

			my @innerFields = split /\|/, $eff;
			#print $innerFields[0], "\n";

			my $cDNA_Position = $innerFields[5];
			my $CDS_Position = $innerFields[6];
			my $aaPosition= $innerFields[7];
			my $aaChange=$innerFields[8];
			
			my @aaChangeFields;
			my $aaFrom;
			my $aaTo;
			if($aaChange)	{
				@aaChangeFields = split /\//, $aaChange;
				$aaFrom = $aaChangeFields[0];
				$aaTo = $aaChangeFields[1];
			}
			my $codonChange = $innerFields[9];
# Perl Line 50
			my $geneName = $innerFields[1];						# NEED TO MAKE A COUNTER TO SAY effect1, 2 ETC
			#if($geneName)		{# Gene name will not exist for all effects so gives error in line below

				# Adding a hyperlink to Ensemblplants:
				# 11.6.2015 - removed this because Martin Trick is now making a webpage/database so hyperlink can be picked up from there.
				#$geneName = "=HYPERLINK(\"http://plants.ensembl.org/Triticum_aestivum/Search/Results?species=Triticum%20aestivum;idx=;q=$geneName\",\"$geneName\")";
			#}

			my $Transcript_ID = $innerFields[2];			
			
			my $effect = $innerFields[4];		# CSQ=T||||intergenic_variant||||||||

			
			
			# Printing each effect on each row and also the extra info that goes with it:
       		my @fieldsSubset = @fields[0..14];				# Gets all fields from the vcf_plus file plus Contig + position for safety
			print join("\t", @fieldsSubset);
			print "\t";
				
#			if($effect =~ m/intergenic_variant|downstream_gene_variant|intron_variant|5_prime_UTR_variant|3_prime_UTR_variant|splice_acceptor_variant/)	{# None of the other fields exist so dont print these
#			
#				print "\t\t", $effect, "\t\t\t";
#			}
#			else	{
				### 28.1.2015 - add here no warnings for this print - to silence the null vars.
				### 11.6.2015 - done this but still can t eliminate some of the warnings
				### 9.10.2015 - it would be good to print each variable inside it s own conditional to avoid this issue - i.e. only print if variable has a value 

      			print do {
       				no warnings "all"; # restricted to this do block only
       				$geneName, "\t", $Transcript_ID, "\t", $effect, "\t", $cDNA_Position, "\t", $CDS_Position, "\t", "${aaFrom}${aaPosition}${aaTo}", "\t", $codonChange, "\t";
       			};
				#print $geneName, "\t", $Transcript_ID, "\t", $effect, "\t", "${aaFrom}${aaPosition}${aaTo}", "\t", $codonChange, "\t";
#			}
			print "\n";	
		}		
		#print "\n";	# Moved above to print one effect per line
	}
' \
> ${file}_formatted.txt
done
***string_marker_for_commenting_out_code***





# Now merge the 2 annotation files, no_new3B and new3B:
### NB - none of these file names has a header because I forgot to change search for header line in the above Perl code
### from Contig to scaffold - but this is OK - acually makes it slightly easier to continue: 
#cat  VEP_-d3_-s3_no_new3B.vcf_plus_info.txt_formatted.txt  VEP_-d3_-s3_new3B.vcf_plus_info.txt_formatted.txt  >  VEP_-d3_-s3.vcf_plus_info_formatted_final.txt




# Once the formatted SNPs have been checked, remove the following extra columns from the table,
# currently $14 (homoeolog number) and $15 (space) - 15.10.2015 - except for the mm MAPS run where I will include an extra column $16(?) which should be printed at the very end
#cat VEP_-d3_-s3.vcf_plus_info_formatted_final.txt \
#| awk -F '\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22}' \
#> wheat_6x_mutations.tab
# Can also check the bottom of the file like so:
###tail -n 500  VEP_-d3_-s3.vcf_plus_info_formatted_final.txt | awk -F '\t' '{print $13 "\t" $14 "\t" $15 "\t" $18}'


#########################
# Fields in output table:
#scaffold	chr	library	line	position	chromosome_position	ref_base	wt_base	alt_base	het/hom	wt_cov	mut_cov	confidence	\
#gene	feature	consequence	cDNA_position	CDS_position	amino_acid	codon	
#########################

# Get stats we need for paper Table 1 and other analyses.
# NB - the maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt 
# file used also contains deletion mutations!

# Total number of mutations (high, medium and low confidence):
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /het/ && $10 >= 6) || ($8 ~ /hom/ && $10 >= 3)' | wc -l 
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /het/ && $10 == 5)' | wc -l
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /het/ && ($10 == 4 || $10 == 3) )' | wc -l
# Total number of homozygous mutations:
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /hom/ && $10 >= 3)' | wc -l
# Total number of EMS homozygous mutations: 
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /hom/ && $10 >= 3)' | awk '$11 ~ /GA|CT/' | wc -l
# Total number of heterozygous mutations:
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /het/ && $10 >= 6)' | wc -l
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /het/ && $10 == 5)' | wc -l
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /het/ && $10 == 4)' | wc -l
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /het/ && ($10 == 4 || $10 == 3) )' | wc -l
# Total number of EMS heterozygous mutations:
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /het/ && $10 >= 6)' | awk '$11 ~ /GA|CT/' | wc -l
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /het/ && $10 >= 5)' | awk '$11 ~ /GA|CT/' | wc -l
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /het/ && $10 >= 4)' | awk '$11 ~ /GA|CT/' | wc -l





# 1.4.2016 - redone to remove duplicate libs   (last time done: 22.2.2016)
# Preparing specific data set for paper and website:
# for the paper - HetMonCov5+. HomMinCov3+:
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /het/ && $10 >= 5) || ($8 ~ /hom/ && $10 >= 3)' | wc -l
#Total # mutations up till now: 6,566,073		(last time:	6,574,935)
# First had to count and remove the indeletion mutations:
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt | awk '($8 ~ /het/ && $10 >= 5) || ($8 ~ /hom/ && $10 >= 3)' |  awk '$5 ~ /[\.\+\*\d]/ || $6 ~ /[\.\+\*\d]/' | wc -l
#5863				(last time: 5,869)
# Then added in the 5DS and 6AL data sets previously absent for PlateK LIB15904-LIB15935 MAPS runs + counted the deletions again:
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt ../PlateK/MAPS_bams_to_use_LIB15904-LIB15935/5DS/maps-part2_-l20_-c10_-d3_-p10_-s2.txt ../PlateK/MAPS_bams_to_use_LIB15904-LIB15935/6AL/maps-part2_-l20_-c10_-d3_-p10_-s2.txt \
#| awk '($8 ~ /het/ && $10 >= 5) || ($8 ~ /hom/ && $10 >= 3)' |  awk '$5 ~ /[\.\+\*\d]/ || $6 ~ /[\.\+\*\d]/' | wc -l
#5,868 					(last time: 5,874) (just a few more - makes sense)
# Adding the new 5DS and 6AL data gives a further 3,645 mutations.
# Now REMOVED the deletions and counted the number of mutations in the final file:
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt \
#../PlateK/MAPS_bams_to_use_LIB15904-LIB15935/5DS/maps-part2_-l20_-c10_-d3_-p10_-s2.txt \
#../PlateK/MAPS_bams_to_use_LIB15904-LIB15935/6AL/maps-part2_-l20_-c10_-d3_-p10_-s2.txt \
#| awk '($8 ~ /het/ && $10 >= 5) || ($8 ~ /hom/ && $10 >= 3)' \
#| awk '$5 !~ /[\.\+\*\d]/ && $6 !~ /[\.\+\*\d]/' | wc -l #OR > maps-part2_Het5+_Hom3+_no_low_pc_EMS_noWTs_no_delns.txt 
# Copied the file to here: /tgac/workarea/group-tg/projects/hexaploid_capture/main_mutations_set/maps-part2_Het5+_Hom3+_no_low_pc_EMS_noWTs_no_deltns.txt
# 6,563,850 mutations across 1,209 lines				(last time: 6,572,706 mutations across 1,209 lines)
# NBNB - halt! - 4.4.2016 - I think i have catted 5DS and 6AL twice because they will also be in the maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt file because this is a new file 
# Still need to remove hetzg muts

# Now place the small deletions into a separate file:
# (NB - i.e. removing * from $5, i.e. deletions that are present in the reference and removing '.', '+' and '\d' from the mutant which also implies a deletion in the reference 
# Just noticed that there are a few mutations that have insertions in Cadenza WT and mutant reads w.r.t reference but then there is a base change within the Cadenza 
# insertion - but none are EMS changes! 
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt \
#../PlateK/MAPS_bams_to_use_LIB15904-LIB15935/5DS/maps-part2_-l20_-c10_-d3_-p10_-s2.txt \
#../PlateK/MAPS_bams_to_use_LIB15904-LIB15935/6AL/maps-part2_-l20_-c10_-d3_-p10_-s2.txt \
#| awk '($8 ~ /het/ && $10 >= 5) || ($8 ~ /hom/ && $10 >= 3)' \
#| awk '$5 ~ /[\.\+\d]/ || $6 ~ /*/' | wc -l #OR > maps-part2_Het5+_Hom3+_no_low_pc_EMS_noWTs_delns_ONLY.txt
# 5,857 mutations in 954 libs - average # mutations per lib = 6; median number of mutations per lib = 3
# Clarify whether Cristobal wants this file or the one below minus deletions in multiple libs.  

# Of these mutations, how many are present only in a single library?
#cat maps-part2_Het5+_Hom3+_no_low_pc_EMS_noWTs_delns_ONLY.txt | awk '{print $1 "\t" $2}'|uniq -c | awk '$1 == 1' | wc -l
# 5090 mutations (86.90%)
# Prepare a list file of the 186 mutations in multiple libraries for removal:
#cat maps-part2_Het5+_Hom3+_no_low_pc_EMS_noWTs_delns_ONLY.txt | awk '{print $1 "\t" $2}'|uniq -c | awk '$1 > 1' | awk '{print $2 "\t" $3}' > deln_muts_in_multiple_libs
# Now remove them from the data:
#cat maps-part2_Het5+_Hom3+_no_low_pc_EMS_noWTs_delns_ONLY.txt | grep -v -f deln_muts_in_multiple_libs > maps-part2_Het5+_Hom3+_no_low_pc_EMS_noWTs_delns_ONLY_uniq2singleLib.txt





# for the website - HetMonCov3+. HomMinCov2+:			[ redoing 4.4.2016 to remove duplicate libs ]
# Adding hom mutations with coverage of 2 plus all of the new 5DS and 6AL data. Also need to remove their headers
#The file for the paper (see above) is OK because it has a numerical check on col $10 which fails for the header.
# Old data before removal of suplicate samples../VEP_analysis_all_samples_11_2_2016_hetMinCov2/maps-part2_allmuts_extractedHomMinCovAt2_no_low_pc_EMS_noWTs_sort_g.txt \
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g.txt \
#../VEP_analysis_all_samples_4_4_2016_hetMinCov2/maps-part2_allmuts_extractedHomMinCovAt2_no_low_pc_EMS_noWTs_sort_g.txt \
#../PlateK/MAPS_bams_to_use_LIB15904-LIB15935/5DS/maps-part2_-l20_-c10_-d3_-p10_-s2.txt \
#../PlateK/MAPS_bams_to_use_LIB15904-LIB15935/6AL/maps-part2_-l20_-c10_-d3_-p10_-s2.txt \
#| grep -v ^Chrom \
#| awk '$5 !~ /[\.\+\*\d]/ && $6 !~ /[\.\+\*\d]/' \
#>  maps-part2_Het3+_Hom2+_no_low_pc_EMS_noWTs_no_delns.txt
#| wc -l 
# 8,989,111 - gives 1209 libs				(last time: 9,001,867	- gives 1210 libs (1 more than predicted - the extra lib is coming from the HomMinCov==2 data)  )
# Copied the file to here: /tgac/workarea/group-tg/projects/hexaploid_capture/main_mutations_set/maps-part2_Het5+_Hom3+_no_low_pc_EMS_noWTs_no_deltns.txt
# Still need to remove hetzg muts  


# For calculating the # het muts at HetMinCov==2 minus WTs, no deltns etc for the %EMS vs coverage graph:
#cat maps-part2_-lX_-c10_-d2_-p10_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt \
#../PlateK/MAPS_bams_to_use_LIB15904-LIB15935/5DS/maps-part2_-l20_-c10_-d2_-p10_-s2.txt \
#../PlateK/MAPS_bams_to_use_LIB15904-LIB15935/6AL/maps-part2_-l20_-c10_-d2_-p10_-s2.txt \
#| grep -v ^Chrom \
#| grep -v "$testSamples54\|$plateH\|$plateI\|$plateN\|$testSamples54_WT\|$CadWTs" \
#| awk '$5 !~ /[\.\+\*\d]/ && $6 !~ /[\.\+\*\d]/' \
#| grep het | awk '$10 == 2' | awk '$11 ~ /GA|CT/' \
#| wc -l






# The calculations for the actual number of stop and missense mutations is best done with the file containing 
# one mutation per line i.e. VEP_-d3_-s3_[no]_new3B.vcf_plus_info.txt.
# Catted both files together:
#cat VEP_-d3_-s3_no_new3B.vcf_plus_info.txt VEP_-d3_-s3_new3B.vcf_plus_info.txt > VEP_-d3_-s3.vcf_plus_info_allSNPs.txt
### NB - the catted file has 2 header lines!!!!


# Total number of missense mutations identified:
#cat VEP_-d3_-s3.vcf_plus_info_allSNPs.txt | grep missense_variant | awk '($10 ~ /het/ && $12 >= 3) || ($10 ~ /hom/ && $12 >= 3) ' | wc -l
#cat VEP_-d3_-s3.vcf_plus_info_allSNPs.txt | grep missense_variant | awk '($10 ~ /het/ && $12 == 5)' | wc -l
#cat VEP_-d3_-s3.vcf_plus_info_allSNPs.txt | grep missense_variant | awk '($10 ~ /het/ && $12 == 4)' | wc -l
#cat VEP_-d3_-s3.vcf_plus_info_allSNPs.txt | grep missense_variant | awk '($10 ~ /het/ && $12 == 3)' | wc -l
# Total number of stop or splice site mutations:
#cat VEP_-d3_-s3.vcf_plus_info_allSNPs.txt | grep stop_gained | awk '($10 ~ /het/ && $12 >= 4) || ($10 ~ /hom/ && $12 >= 3) ' | wc -l
#cat VEP_-d3_-s3.vcf_plus_info_allSNPs.txt | grep stop_gained | awk '($10 ~ /het/ && $12 == 5)' | wc -l
#cat VEP_-d3_-s3.vcf_plus_info_allSNPs.txt | grep stop_gained | awk '($10 ~ /het/ && $12 == 4)' | wc -l
#cat VEP_-d3_-s3.vcf_plus_info_allSNPs.txt | grep stop_gained | awk '($10 ~ /het/ && $12 == 3)' | wc -l


# Get number of genes that have 1 or more stop or splice site changes ( I think these values can be found with wheat_6x_mutations.tab
# because I'm uniqifying the genes before counting):
#cat wheat_6x_mutations.tab | awk '($10 ~ /het/ && $12 >= 6) || ($10 ~ /hom/ && $12 >= 3)'| awk -F '\t' '$16 ~ /stop_gained/ || $16 ~ /splice_donor_variant/ || $16 ~ /splice_acceptor_variant/' | awk -F '\t' '{print $14}' | sort -u > wheat_6x_mutations.tab_-d6_stop_gained_gene_list
### 6.7.2015 - 33935 - had left out /splice_region_variant/ - when include this as well get 44,413 BUT Cristobal doesn't want this category the exon bases within the splice_region_variant will appear in the stop_gained category anyway.
# Now need to find the extra new genes at medium and low confidence that aren't in the high confidence set
#cat wheat_6x_mutations.tab | awk '$10 ~ /het/ && $12 == 5' | awk -F '\t' '$16 ~ /stop_gained/ || $16 ~ /splice_donor_variant/ || $16 ~ /splice_acceptor_variant/' | awk -F '\t' '{print $14}' | sort -u > wheat_6x_mutations.tab_-d5_stop_gained_gene_list
#cat wheat_6x_mutations.tab | awk '$10 ~ /het/ && ($12 == 4)' | awk -F '\t' '$16 ~ /stop_gained/ || $16 ~ /splice_donor_variant/ || $16 ~ /splice_acceptor_variant/' | awk -F '\t' '{print $14}' | sort -u > wheat_6x_mutations.tab_-d4_stop_gained_gene_list
#cat wheat_6x_mutations.tab | awk '$10 ~ /het/ && ($12 == 3)' | awk -F '\t' '$16 ~ /stop_gained/ || $16 ~ /splice_donor_variant/ || $16 ~ /splice_acceptor_variant/' | awk -F '\t' '{print $14}' | sort -u > wheat_6x_mutations.tab_-d3_stop_gained_gene_list
# Get counts:
# Number of high confidence SNPs with stop_gained:
#cat wheat_6x_mutations.tab_-d6_stop_gained_gene_list | wc -l
# Number of medium confidence SNPs with stop_gained:
#comm -13 wheat_6x_mutations.tab_-d6_stop_gained_gene_list wheat_6x_mutations.tab_-d5_stop_gained_gene_list | wc -l
# Number of low confidence SNPs with stop_gained:
#cat wheat_6x_mutations.tab_-d6_stop_gained_gene_list wheat_6x_mutations.tab_-d5_stop_gained_gene_list | sort -u > wheat_6x_mutations.tab_-d6_d5_stop_gained_gene_list
#comm -13 wheat_6x_mutations.tab_-d6_d5_stop_gained_gene_list wheat_6x_mutations.tab_-d4_stop_gained_gene_list | wc -l
#cat wheat_6x_mutations.tab_-d6_stop_gained_gene_list wheat_6x_mutations.tab_-d5_stop_gained_gene_list wheat_6x_mutations.tab_-d4_stop_gained_gene_list | sort -u > wheat_6x_mutations.tab_-d6_-d5_-d4_stop_gained_gene_list
#comm -13 wheat_6x_mutations.tab_-d6_-d5_-d4_stop_gained_gene_list wheat_6x_mutations.tab_-d3_stop_gained_gene_list | wc -l


# Get number of genes that have 1 or more missense/non-synonymous changes:
#cat wheat_6x_mutations.tab | grep missense_variant | awk '($10 ~ /het/ && $12 >= 6) || ($10 ~ /hom/ && $12 >= 3)' | awk -F '\t' '{print $14}' | sort -u > wheat_6x_mutations.tab_-d6_missense_variant_gene_list
# Now need to find the extra new genes at medium and low confidence that aren't in the high confidence set:
#cat wheat_6x_mutations.tab | grep missense_variant | awk '$10 ~ /het/ && $12 == 5' | awk -F '\t' '{print $14}' | sort -u > wheat_6x_mutations.tab_-d5_missense_variant_gene_list
#cat wheat_6x_mutations.tab | grep missense_variant | awk '$10 ~ /het/ && ($12 == 4 || $12 == 3)' | awk -F '\t' '{print $14}' | sort -u > wheat_6x_mutations.tab_-d4_3_missense_variant_gene_list
#cat wheat_6x_mutations.tab | grep missense_variant | awk '$10 ~ /het/ && ($12 == 4)' | awk -F '\t' '{print $14}' | sort -u > wheat_6x_mutations.tab_-d4_missense_variant_gene_list
#cat wheat_6x_mutations.tab | grep missense_variant | awk '$10 ~ /het/ && ($12 == 3)' | awk -F '\t' '{print $14}' | sort -u > wheat_6x_mutations.tab_-d3_missense_variant_gene_list
# Get counts:
# Number of high confidence SNPs with missense changes:
#cat wheat_6x_mutations.tab_-d6_missense_variant_gene_list | wc -l
# Number of medium confidence SNPs with missense changes:
#comm -13 wheat_6x_mutations.tab_-d6_missense_variant_gene_list wheat_6x_mutations.tab_-d5_missense_variant_gene_list | wc -l
# Number of low confidence SNPs with missense changes:
#cat  wheat_6x_mutations.tab_-d6_missense_variant_gene_list wheat_6x_mutations.tab_-d5_missense_variant_gene_list | sort -u > wheat_6x_mutations.tab_-d6_d5_missense_variant_gene_list
#comm -13 wheat_6x_mutations.tab_-d6_d5_missense_variant_gene_list wheat_6x_mutations.tab_-d4_3_missense_variant_gene_list | wc -l
# Low confidence at hetMinCov >= 4:
#cat wheat_6x_mutations.tab_-d6_missense_variant_gene_list wheat_6x_mutations.tab_-d5_missense_variant_gene_list wheat_6x_mutations.tab_-d4_stop_gained_gene_list| sort -u | wc -l

# Number of unique positions with a mutation (at HetMinCov >= 4):
#cat wheat_6x_mutations.tab | awk '($10 ~ /het/ && $12 >= 4) || ($10 ~ /hom/ && $12 >= 3)' | awk '{print $1 "\t" $5 "\t" $6}'| sort -u | wc -l
# Number of effects (at HetMinCov >= 4):
#cat wheat_6x_mutations.tab | awk '($10 ~ /het/ && $12 >= 4) || ($10 ~ /hom/ && $12 >= 3)' | wc -l
# Number of libraries:
#cat wheat_6x_mutations.tab | awk '{print $3}' | sort -u | wc -l

# Count the number of SNPs on subgenome A, B and D (at HetMinCov >= 4; no 3B):
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g_no_new3B.txt | awk '($8 ~ /het/ && $10 >= 4) || ($8 ~ /hom/ && $10 >= 3)' | grep IWGSC_CSS_.A._ | wc -l
#1326137
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g_no_new3B.txt | awk '($8 ~ /het/ && $10 >= 4) || ($8 ~ /hom/ && $10 >= 3)' | grep IWGSC_CSS_.D._ | wc -l
#1237728
#cat maps-part2_3x_4x_5x_more6x_lines_no_low_pc_EMS_noWTs_sort_g_no_new3B.txt | awk '($8 ~ /het/ && $10 >= 4) || ($8 ~ /hom/ && $10 >= 3)' |grep IWGSC_CSS_.B._ | wc -l
#1111708





















##############
# Old analysis
##############

# The above file is too big for Excel - only 1048576 lines are accepted into Excel, but
# this number of rows also causes problems with editing table e.g. deleting columns.
# So printing out blocks of 50,000 lines to separate files:
### 20.3.2015 - could split this up further - 25,000 
### 26.5.2015 - the file is now too big to do ecen this!!
#head -n 500000 VEP_-d3_-s3.vcf_plus_info_formatted.txt > VEP_-d3_-s3.vcf_plus_info_formatted_part1.txt
#sed -n '500001,1000000p; 1000001q' VEP_-d3_-s3.vcf_plus_info_formatted.txt > VEP_-d3_-s3.vcf_plus_info_formatted_part2.txt
#sed -n '1000001,1500000p; 1500001q' VEP_-d3_-s3.vcf_plus_info_formatted.txt > VEP_-d3_-s3.vcf_plus_info_formatted_part3.txt
### Do we need the split???
#sed -n '1000001,1500000p; 1500001q' VEP_-d3_-s3.vcf_plus_info_formatted.txt > VEP_-d3_-s3.vcf_plus_info_formatted_part3.txt
#sed -n '1000001,1500000p; 1500001q' VEP_-d3_-s3.vcf_plus_info_formatted.txt > VEP_-d3_-s3.vcf_plus_info_formatted_part3.txt
#sed -n '1000001,1500000p; 1500001q' VEP_-d3_-s3.vcf_plus_info_formatted.txt > VEP_-d3_-s3.vcf_plus_info_formatted_part3.txt
#sed -n '1000001,1500000p; 1500001q' VEP_-d3_-s3.vcf_plus_info_formatted.txt > VEP_-d3_-s3.vcf_plus_info_formatted_part3.txt
#sed -n '1000001,1500000p; 1500001q' VEP_-d3_-s3.vcf_plus_info_formatted.txt > VEP_-d3_-s3.vcf_plus_info_formatted_part3.txt
### 20.3.2015 - with more mutations calls, will have to print out more chunks!
# Import into Excel - press hyperlink in the Styles area to show hyperlink underlined in blue.










#Results - all SNPs:
#mapsFileAllSNPs=maps-part2_-l29_-c10_-d6_-p20_-s4_all_plate_weird_samples_rm.txt
#snpEffFile=snpEff_maps-part2_-l29_-c10_-d6_-p20_-s4_all_plate_weird_samples_rm.vcf	# snpEff_maps-part2_-l30_-c10_-d6_-p20_-s4.vcf

# Total number of mutant lines: 
#numbrMutantLines=`more $mapsFileAllSNPs | awk '{print $7}' | sort -u | wc -l`

#Total number of high confidence mutations:
#totalMutations=`more $mapsFileAllSNPs | wc -l` 
#echo "Total number of high confidence mutations: $totalMutations"

#Total number of heterozygous mutations:
#totalHet=`more $mapsFileAllSNPs | grep 'het' | wc -l` 
#echo "Total number of heterozygous mutations: $totalHet"

#Total number of homozygous mutations:
#totalHom=`more $mapsFileAllSNPs | grep 'hom' | wc -l` 
#echo "Total number of homozygous mutations: $totalHom"


#Total number of annotated genes = 108569
#source jre-7.11; java -Xmx4g -jar ~baileyp/ProgramFiles/snpEff_v4.0/snpEff.jar dump -v -txt Triticum_aestivum.IWGSP1.21 > Triticum_aestivum.IWGSP1.21.txt		

#Total number of effects = 128443		423218 (all samples to date)
#more $file | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | wc -l	
#Total number of predicted non-synonymous alleles				
#Total number of predicted stop codons gained & splice site alleles

#Total number of EFFECTS with non-synonymous changes = 34125
#more $file | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep missense_variant | wc -l
	
#Total number of EFFECTS with a splice site change or early stop codon = 2190 				
# acceptor = 317 		
#more snpEff_maps-part2_-l30_-c10_-d6_-p20_-s4.vcf | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep splice_acceptor_variant | wc -l
#more snpEff_maps-part2_-l30_-c10_-d6_-p20_-s4.vcf | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep splice_acceptor_variant | wc -l			
# donor = 292
#more snpEff_maps-part2_-l30_-c10_-d6_-p20_-s4.vcf | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep splice_donor_variant | wc -l
#stop = 1581
#more snpEff_maps-part2_-l30_-c10_-d6_-p20_-s4.vcf | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep stop_gained | wc -l
#% all genes [c] = 2 %				

#The above EFFECTS should be the same as for hans's calculation - check it is:
#@Paul: For those two numbers, I believe it's just the count of records in the combined VEP VCF file.  

#So, for non-synonmous:
#campus-018-171:current_combined_maps havasquezgross$ cut -f 8 Kronos.mapspart2.HetMinCov6HetMinPer15.ktc1-12_ktc15-81.vep.vcf | cut -f 5 -d "|" | sort | uniq -c | egrep "missense" | awk '{ sum += $1 } END { print sum }'
#340104

#Or for stop codons gained/splice site alleles:
#campus-018-171:current_combined_maps havasquezgross$ cut -f 8 Kronos.mapspart2.HetMinCov6HetMinPer15.ktc1-12_ktc15-81.vep.vcf | cut -f 5 -d "|" | sort | uniq -c | egrep "stop_gained|splice" | awk '{ sum += $1 } END { print sum }'
#53469


# Total number of GENES with with non-synonymous changes: 
#genesNonSyn=`more $snpEffFile | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep missense_variant | awk -F '|' '{print $6}' | sort -u | wc -l`
#echo "Total number of GENES with with non-synonymous changes: $genesNonSyn" 
#echo "	as a % of all genes:" `awk -vgenesNonSyn=$genesNonSyn 'BEGIN{printf "%.2f", (genesNonSyn/108569)*100}'` '%'



#Total number of GENES with a splice site change or early stop codon:			
# acceptor, donor and/or stop	
#acceptorDonorStop=`more $snpEffFile | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep 'splice_acceptor_variant\|splice_donor_variant\|stop_gained' | awk -F '|' '{print $6}' | sort -u | wc -l`	
# donor
#donor=`more $snpEffFile | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep splice_donor_variant | awk -F '|' '{print $6}' | sort -u | wc -l`
##stop = 1581
#stop=`more $snpEffFile | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep stop_gained | awk -F '|' '{print $6}' | sort -u | wc -l`
#% all genes [c] = 2 %				
#echo "Total number of GENES with a splice site change and/or early stop codon: $acceptorDonorStop"
#echo "	as a % of all genes:" `awk -vacceptorDonorStop=$acceptorDonorStop 'BEGIN{printf "%.2f", (acceptorDonorStop/108569)*100}'` '%'