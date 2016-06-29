#!/bin/bash
#########################################
# 35.ExonCapture_300_Samples_GetCounts.sh

# Paul bailey		24.1.2015
#########################################

# This script uses the hard coded values - check them before use!
#onTargSiz=95373308					# Original value for CSS v1 reference = 102152087; IWGSC_v2_ChrU.fa = 95373308 (minus any 3B contigs) 
#basesInGenomeRef=9500075743		# Number of bases in the CSS v1 reference = 10138701012; IWGSC_v2_ChrU.fa = 9500075743 (minus any 3B contigs)
# 7.6.2016 Updated these values for IWGSC_v2_ChrU.fa and a new on-target bedfile: 
basesInGenomeRef=10393511112		# for IWGSC_v2_ChrU.fa genome reference
onTargSizPlusPadding=116624907		# for "tophit" exon query bedfile (plus padding)
onTargSizMinusPadding=60369289		# for "tophit" exon query bedfile (minus padding)
 


printf '%b' "Lane\tSampleID\tDescription\tReadMates\tCadenzaNo\t \
readsR1\treadsR2\tadaptorContamR1_%\tadaptorContamR2_%\tgoodQualTrimR1_%FullLen\tgoodQualTrimR2_%FullLen\t \
BWA_FragLen_bp\tDupRate_%\tOpticalDups\tOpticalDups_%\tNon-dupReads\tNon-dupReads_%\t \
Non-dupReads_q20\tNon-dupReads_q20_%\t \
MatesInProperReadPairs\tProperReadPairs_%\t \
OnTargReadMates0bp_padding\tOnTargReadMates0bp_padding_%\t \
OnTargReadMates200bp_padding\tOnTargReadMates200bp_padding_%\t \
AvOn-targCovPerBase0bp_padding\tOn-targBases1xCov0bp_padding\tOn-targBases1xCov0bp_padding_%\tOn-targBases6xCov0bp_padding\tOn-targBases6xCov0bp_padding_%\tOn-targBases10xCov0bp_padding\tOn-targBases10xCov0bp_padding_%\t \
On-targBases15xCov0bp_padding\tOn-targBases15xCov0bp_padding_%\tOn-targBases20xCov0bp_padding\tOn-targBases20xCov0bp_padding_%\tOn-targBases30xCov0bp_padding\tOn-targBases30xCov0bp_padding_%\tOn-targBases100xCov0bp_padding\tOn-targtBases100xCov0bp_padding_%\t \
AvOn-targCovPerBase200bp_padding\tOn-targBases1xCov200bp_padding\tOn-targBases1xCov200bp_padding_%\tOn-targBases6xCov200bp_padding\tOn-targBases6xCov200bp_padding_%\tOn-targBases10xCov200bp_padding\tOn-targBases10xCov200bp_padding_%\t \
On-targBases15xCov200bp_padding\tOn-targBases15xCov200bp_padding_%\tOn-targBases20xCov200bp_padding\tOn-targBases20xCov200bp_padding_%\tOn-targBases30xCov200bp_padding\tOn-targBases30xCov200bp_padding_%\tOn-targBases100xCov200bp_padding\tOn-targtBases100xCov200bp_padding_%\t \
GenomeCovBases1xCov\tGenomeCovBases1xCov_%\tGenomeCovBases6xCov\tGenomeCovBases6xCov_%\tGenomeCovBases10xCov\tGenomeCovBases10xCov_%\tGenomeCovBases20xCov\tGenomeCovBases20xCov_%\n" \
> ExonCapture_alignment_results.txt
#On-targReadPairs_q>=20\tOn-targReadPairs_q>=20_%\t
#> ExonCapture_alignment_results_52_test_samples.txt
#> ExonCapture_alignment_results.txt
#>  35.ExonCapture_alignment_results.txt

# For looking at the qPCR tests:
#printf '%b' "Lane\tSampleID\tDescription \
#avOnTargCovPerBase	avOnTargCovperBase_Norm	avOnTargCovPerBase_IWGSP1_EnsemblPlnts22	avOnTargCovPerBase_JoseRef	\
#avOnTargCovPerBase_Rubsico_MDH	avOnTargCovPerBase_Rubsico_MDH_Norm	Rubsico_MDH_EF	\
#GenomeCovBases1xCov_%	\
#avOffTargCovPerBase	OnTarg_EF	avOffTargCovPerBase_IWGSP1_EnsemblPlnts22	OnTarg_IWGSP1_EnsemblPlnts22_EF	avOffTargCovPerBase_JoseRef	OnTarg_JoseRef_EF"
#> 35.ExonCapture_alignment_results.txt
### NB - this header doesn't print! Why not?





#plateToDo='PlateA\|PlateB\|PlateC\|PlateD\|PlateE\|PlateF\|PlateG'
plateToDo='PlateH\|PlateI\|PlateJ\|PlateK\|PlateL\|PlateM\|PlateN'
#plateToDo='54_test_samples'
#plateToDo='-v LIB5780\|LIB5782'
#tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_300_Samples/35.ExonCapture_PreparingTable.txt | grep $plateToDo |
#tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/54_test_samples/35.ExonCapture_PreparingTable_54_TestSamples.txt | grep $plateToDo |
#tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/ExonCapture_PlateH2N.txt | grep $plateToDo |
tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/VEP_analysis_aux_files/ExonCapture_PlateTABCDEFGHIJKLMN.tab | grep $plateToDo |
while read line; do

	lane=`echo $line | cut -d ' ' -f 1`	
	sampleID=`echo $line | cut -d ' ' -f 2` 
	description=`echo $line | cut -d ' ' -f 3` 
	numbrReads=`echo $line | cut -d ' ' -f 4`
	paths_to_libs=`echo $line | cut -d ' ' -f 5` 
	CadenzaNo=`echo $line | cut -d ' ' -f 6`	# | sed s/\n//`
	
	# Prepare basic library name:
	lib=`echo $sampleID | sed s/^...._// | sed s/_LDI....$// `
	
	# Total number of reads:
	readsR1=`cat $paths_to_libs/R1_scythe.stderr | grep 'contaminated:' | awk '{print $6}'`
	readsR2=`cat $paths_to_libs/R2_scythe.stderr | grep 'contaminated:' | awk '{print $6}'`
	
	# Adaptor contamination rates:
	adaptorContamR1=`cat $paths_to_libs/R1_scythe.stderr | grep 'contamination rate:' | awk '{print $3}'`
	adaptorContamR1Pc=`awk -vadaptorContamR1=$adaptorContamR1 'BEGIN{printf "%.2f", adaptorContamR1 * 100}'`
	adaptorContamR2=`cat $paths_to_libs/R2_scythe.stderr | grep 'contamination rate:' | awk '{print $3}'`
	adaptorContamR2Pc=`awk -vadaptorContamR2=$adaptorContamR2 'BEGIN{printf "%.2f", adaptorContamR2 * 100}'`
	
	# Number of full length reads after quality trimming:
	qualTrimR1=`cat $paths_to_libs/trimmed_R1_len.txt | grep '^101' | awk '{print $2}'`
	goodQualTrimR1Pc=`awk -vreadsR1=$readsR1 -vqualTrimR1=$qualTrimR1 'BEGIN{printf "%.2f", (qualTrimR1/readsR1) * 100}'`
	qualTrimR2=`cat $paths_to_libs/trimmed_R2_len.txt | grep '^101' | awk '{print $2}'` 
	goodQualTrimR2Pc=`awk -vreadsR2=$readsR2 -vqualTrimR2=$qualTrimR2 'BEGIN{printf "%.2f", (qualTrimR2/readsR2) * 100}'`
	
	### mapped read count (includes proper pairs and single mates):
	###tail -n 1  $lib/sampe_mapped_count.log | grep -P '^\d'

	# look at the fragment length sizes for each library:
	### 4.8.2015 - NB - it would probably be best to issue this command in the main script and save the result in a separate file
	### which can then be counted here:
	fragSize=`head -n 3000 $paths_to_libs/sampe_sort_markdup_rm_filter_-f2.log | grep 'inferred external isize from' | head -n 1| awk '{print $8}'`
	# For the UV's:
#	fragSize=`head -n 3000 $paths_to_libs/LIB*.e* | grep 'inferred external isize from' | head -n 1| awk '{print $8}'`
	
	# Get Picard's PERCENT_DUPLICATION column:
	dupRate=`cat $paths_to_libs/sampe_sort_markdup_metrics | grep '^Unknown' | awk '{print $9*100}'`
	
	
	# Get Picard's READ_PAIR_OPTICAL_DUPLICATES column:
	opticalDups=`cat $paths_to_libs/sampe_sort_markdup_metrics | grep '^Unknown' | awk '{print $8}'`
	opticalDupsPc=`awk -vnumbrReads=$numbrReads -vopticalDups=$opticalDups 'BEGIN{printf "%.2f", (opticalDups/numbrReads)*100}' `

		
	# Number of non-duplicated reads:
	nonDupReads=`tail -n 1  $paths_to_libs/sampe_sort_markdup_rm_count.log | grep -P '^\d'`
	nonDupReadsPc=`awk -vnumbrReads=$numbrReads -vnonDupReads=$nonDupReads 'BEGIN{printf "%.2f", (nonDupReads/numbrReads)*100}' `

	
	# Number of non-duplicated reads with mapping quality >= 20:
	nonDupReads_q20=`tail -n 1  $paths_to_libs/sampe_sort_markdup_rm_st_-F4_-q20_count.log | grep -P '^\d'`
	nonDupReads_q20Pc=`awk -vnumbrReads=$numbrReads -vnonDupReads_q20=$nonDupReads_q20 'BEGIN{printf "%.2f", (nonDupReads_q20/numbrReads)*100}' `
	

	#Number of mates in concordantly aligned read pairs:
	alnReadMates=`tail -n 1  $paths_to_libs/sampe_sort_markdup_rm_filter_-f2_count.log | grep -P '^\d' `
	alnReadPairsPc=`awk -vnumbrReads=$numbrReads -valnReadMates=$alnReadMates 'BEGIN{printf "%.2f", (alnReadMates/numbrReads)*100}' `

	
	# Number read mates on-target (no padding) - Samtools -L:
	onTargReadMates0=`tail -n 1  $paths_to_libs/sampe_sort_markdup_rm_st_-L_0padding_count.log | grep -P '^\d'`
	onTargReadMates0Pc=`awk -vnumbrReads=$numbrReads -vonTargReadMates0=$onTargReadMates0 'BEGIN{printf "%.2f", (onTargReadMates0/numbrReads)*100}' `
	

	# Number read mates on-target (plus padding) - Samtools -L:
	onTargReadMates200=`tail -n 1  $paths_to_libs/sampe_sort_markdup_rm_st_-L_200padding_count.log | grep -P '^\d'`
	onTargReadMates200Pc=`awk -vnumbrReads=$numbrReads -vonTargReadMates200=$onTargReadMates200 'BEGIN{printf "%.2f", (onTargReadMates200/numbrReads)*100}' `
	

	# Av_cov and on-target coverages (no padding):
	avCov_Ontarg0=`tail -n 1  $paths_to_libs/sampe_sort_markdup_rm_covBed_-d_av_cov_more_1_6_10_15_20_25_100x_covs_0padding.log`
	avCov_Ontarg0Pc=`echo $avCov_Ontarg0 | awk -vonTargSiz=$onTargSizMinusPadding '{print $1 \
	"\t" $2 "\t" ($2/onTargSiz)*100 "\t" $3 "\t" ($3/onTargSiz)*100 "\t" $4 "\t" ($4/onTargSiz)*100 "\t" $5 "\t" ($5/onTargSiz)*100 "\t" \
	$6 "\t" ($6/onTargSiz)*100 "\t" $7 "\t" ($7/onTargSiz)*100 "\t" $8 "\t" ($8/onTargSiz)*100 }' `

	
	# Av_cov and on-target coverages (plus padding) :
	avCov_Ontarg200=`tail -n 1  $paths_to_libs/sampe_sort_markdup_rm_covBed_-d_av_cov_more_1_6_10_15_20_25_100x_covs_200padding.log`
	avCov_Ontarg200Pc=`echo $avCov_Ontarg200 | awk -vonTargSiz=$onTargSizPlusPadding '{print $1 \
	"\t" $2 "\t" ($2/onTargSiz)*100 "\t" $3 "\t" ($3/onTargSiz)*100 "\t" $4 "\t" ($4/onTargSiz)*100 "\t" $5 "\t" ($5/onTargSiz)*100 "\t" \
	$6 "\t" ($6/onTargSiz)*100 "\t" $7 "\t" ($7/onTargSiz)*100 "\t" $8 "\t" ($8/onTargSiz)*100 }' `
	

	# genomecov coverage results:
	genomecov=`tail -n 1  $paths_to_libs/sampe_sort_markdup_rm_bt_genomecov_-d_more_1_6_10_20x_covs.log`
	genomecovPc=`echo $genomecov | awk -vbasesInGenomeRef=$basesInGenomeRef '{print $1 "\t" ($1/basesInGenomeRef)*100 \
	"\t" $2 "\t" ($2/basesInGenomeRef)*100 "\t" $3 "\t" ($3/basesInGenomeRef)*100 "\t" $4 "\t" ($4/basesInGenomeRef)*100 }'` 

		

	# Values for the qPCR tests: 
#	avCov_Ontarg=`awk -vnumbrReads=$numbrReads '{print $1}' $paths_to_libs/sampe_sort_markdup_rm_filter_-f2_covBed_orig_ref_more_1_6_10_15_20_25_100x_covs.log`
#	avCov_Ontarg_Norm=`awk -vavCov_Ontarg=$avCov_Ontarg -vnumbrReads=$numbrReads 'BEGIN{printf "%.2f", (avCov_Ontarg/numbrReads)*10000000}'`
	#	avCov_Ontarg_IWGSP1_EnsemblPlnts22=`awk '{print $1}' $paths_to_libs/sampe_sort_markdup_rm_filter_-f2_covBed_IWGSP1_EnsemblPlnts22_ref_more_1_6_10_15_20_25_100x_covs.log`
	#	avCov_Ontarg_JoseRef=`awk '{print $1}' $paths_to_libs/sampe_sort_markdup_rm_filter_-f2_covBed_Jose_ref_more_1_6_10_15_20_25_100x_covs.log`
#	avCov_Ontarg_Rubsico5ABDL_MDH=`awk '{print $1}' $paths_to_libs/sampe_sort_markdup_rm_filter_-f2_bedT_coverage_-d_for_qPCR.log`
#	avCov_Ontarg_Rubsico5ABDL_MDH_Norm=`awk -vavCov_Ontarg=$avCov_Ontarg_Rubsico5ABDL_MDH -vnumbrReads=$numbrReads 'BEGIN{printf "%.2f", (avCov_Ontarg/numbrReads)*10000000}'`
	
	#	genome1xcovPc=`echo $genomecov | awk -vbasesInGenomeRef=10138701012 '{print ($1/basesInGenomeRef)*100}'`
	
#	avCov_Offtarg=`awk '{print $1}' $paths_to_libs/sampe_sort_markdup_rm_filter_-f2_covBed_off-targ_orig_ref_more_1_6_10_15_20_25_100x_covs.log`
#	avCov_Offtarg_EF=`awk -vavCov_Ontarg=$avCov_Ontarg -vavCov_Offtarg=$avCov_Offtarg 'BEGIN{printf "%.2f", (avCov_Ontarg/avCov_Offtarg)}'`
#	avCov_Offtarg_Rubsico5ABDL_MDH_EF=`awk -vavCov_Ontarg=$avCov_Ontarg_Rubsico5ABDL_MDH -vavCov_Offtarg=$avCov_Offtarg 'BEGIN{printf "%.2f", (avCov_Ontarg/avCov_Offtarg)}'`
	#	avCov_Offtarg_IWGSP1_EnsemblPlnts22=`awk '{print $1}' $paths_to_libs/sampe_sort_markdup_rm_filter_-f2_covBed_off-targ_IWGSP1_EnsemblPlnts22_ref_more_1_6_10_15_20_25_100x_covs.log`
	#	avCov_Offtarg_IWGSP1_EnsemblPlnts22_EF=`awk -vavCov_Ontarg=$avCov_Ontarg_IWGSP1_EnsemblPlnts22 -vavCov_Offtarg=$avCov_Offtarg_IWGSP1_EnsemblPlnts22 'BEGIN{printf "%.2f", (avCov_Ontarg/avCov_Offtarg)}'`
	#	avCov_Offtarg_JoseRef=`awk '{print $1}' $paths_to_libs/sampe_sort_markdup_rm_filter_-f2_covBed_off-targ_Jose_ref_more_1_6_10_15_20_25_100x_covs.log`
	#	avCov_Offtarg_IWGSP1_EnsemblPlnts22_EF=`awk -vavCov_Ontarg=$avCov_Ontarg_JoseRef -vavCov_Offtarg=$avCov_Offtarg_JoseRef 'BEGIN{printf "%.2f", (avCov_Ontarg/avCov_Offtarg)}'`
	
	

	#	echo "$sampleID	$avCov_Ontarg	$avCov_Ontarg_Norm	$avCov_Ontarg_IWGSP1_EnsemblPlnts22	$avCov_Ontarg_JoseRef	\
	#	$avCov_Ontarg_Rubsico5ABDL_MDH	$avCov_Ontarg_Rubsico5ABDL_MDH_Norm	$avCov_Offtarg_Rubsico5ABDL_MDH_EF	\
	#	$genome1xcovPc	\
	#	$avCov_Offtarg	$avCov_Offtarg_EF	$avCov_Offtarg_IWGSP1_EnsemblPlnts22	$avCov_Offtarg_IWGSP1_EnsemblPlnts22_EF	$avCov_Offtarg_JoseRef	$avCov_Offtarg_IWGSP1_EnsemblPlnts22_EF"

# Print the data I need here:
echo "$lane	$sampleID	$description	$numbrReads	$CadenzaNo	\
$readsR1	$readsR2	$adaptorContamR1Pc	$adaptorContamR2Pc	$goodQualTrimR1Pc	$goodQualTrimR2Pc	\
$fragSize	$dupRate	$opticalDups	$opticalDupsPc	$nonDupReads	$nonDupReadsPc	\
$nonDupReads_q20	$nonDupReads_q20Pc	\
$alnReadMates	$alnReadPairsPc	\
$onTargReadMates0	$onTargReadMates0Pc	\
$onTargReadMates200	$onTargReadMates200Pc	\
$avCov_Ontarg0Pc	\
$avCov_Ontarg200Pc	\
$genomecovPc"
#$avCov_Ontarg	$avCov_Ontarg_Norm	$avCov_Offtarg	$avCov_Offtarg_EF	\
#$avCov_Ontarg_Rubsico5ABDL_MDH	$avCov_Ontarg_Rubsico5ABDL_MDH_Norm	$avCov_Offtarg_Rubsico5ABDL_MDH_EF"
# $onTargReadPairs_q20  $onTargReadPairs_q20Pc
done \
>> ExonCapture_alignment_results.txt
#>> ExonCapture_alignment_results_52_test_samples.txt
#>> 35.ExonCapture_alignment_results.txt
#>> 35.ExonCapture_alignment_results_qPCR_cols_only.txt


# Transfer to Excel and table should be ready to send.
# Now just adding the rows for the new data to a global 
#table so that the red marks are not destroyed each time 
# data is added
# Look for any unusual patterns and mark in red.

### To do:
### 1. Some of the frag lengths aren't coming in so need to increase head -n - done - still check

### 2. Change Reads to ReadPairs but will need to change other calculations


