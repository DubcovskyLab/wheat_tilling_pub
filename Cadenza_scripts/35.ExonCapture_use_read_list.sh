#!/bin/bash
#######################################################################
# 35.ExonCapture_300_Samples.sh			Modified from 34.ExonCapture.sh

# Paul Bailey	27.10.2014
#######################################################################

# Paths and dirs:
readsPath=/tgac/data/reads
#indexFile=/tgac/workarea/group-cg/baileyp/WheatLoLa/33.ExonCapture/indices/IWGSC_CSS_all_scaff_v1_bwa_index										  				# IWGSC_v1 - the original CSS assembly file
#indexFile=/tgac/workarea/group-cg/baileyp/IWGSCv1.0_to_IWGSCv2.0/IWGSCv2.0.fa_bwa_index														  					# IWGSC_v2 - my version containing new chr 3B+ seqs
indexFile=/tgac/workarea/group-cg/baileyp/IWGSC_v2_ChrU_Ref/IWGSC_v2_ChrU_bwa_index												  								# IWGSC_v2 - my version containing new chr 3B+ seqs + ChrU assembly
#indexFile=/tgac/workarea/group-cg/baileyp/EnsemblPlantsRelease25_added_CSS_chr3B_plus/EnsemblPlantsRelease25_IWGSC2_and_missing_CSS_in3BSEQ_bwa_index 				# IWGSC_v2 - Ensembl Plants version (Release 25) - contains the missing CSS chr 3B+ seqs
#indexFile=/tgac/workarea/group-tg/reference/wheat-genomes/tgac_cs_wgs_discovar/tgac_cs_wgs_discovar_pre-release_v1/discovar_5th_scaffolds_Nremapped.fa_bwa_index	# tgac_cs_wgs_discovar_pre-release_v1 assembly
#indexFile=/tgac/workarea/group-tg/reference/wheat-genomes/tgac_cs_wgs_discovar/Triticum_aestivum_CS42_TGACv1/Triticum_aestivum_CS42_TGACv1_all.fa_bwa_index			# Triticum_aestivum_CS42_TGACv1 plublic release v1 assembly - not used yet




# Building a BWA index on all wheat CSS contigs (RUN TIME: 7.75 hours):
# Reference fasta files:
# 1. Jon's final file = /tgac/references/external/projects/iwgsc/css/IWGSC_CSS_all_scaff_v1.fa
# 2. EnsemblPlantsRelease25_IWGSC2_and_missing_CSS_in3BSEQ.fa
# 3. IWGSC_v2_ChrU.fa
#bsub -J bwa_index -q Prod128 -R "rusage [mem=8000]" -oo bwa_index.log " \
#source testing_env; source bwa-0.7.7; bwa index -a bwtsw -p IWGSC_v2_ChrU_bwa_index \
#IWGSC_v2_ChrU.fa
#"

# Building a Bowtie2 index on all wheat CSS contigs (RUN TIME: ***************7.5*************** hours):
#prefix=IWGSC_CSS_all_scaff_v1			# NB - Jon's final file
#bsub -J ${prefix}_bt2_index -q Prod128 -R "rusage [mem=10000]" -oo ${prefix}_bt2_index.log " \
#source testing_env; source bowtie-2.2.1; bowtie2-build \
#/tgac/references/external/projects/iwgsc/css/${prefix}.fa \
#${prefix}_bt2_index \
#"





#####################################################################
# Mutation/SNP Calling Pipeline for concordantly aligned reads pairs: 
#####################################################################

# Explanation of steps:

# 1a.
# Adaptor and quality trimming (RUN TIME: up to 2 hours 20 minutes)
PRIOR=0.4
QUAL_THRESH=20
ADAPTERS=/tgac/workarea/group-cg/baileyp/Illumina_adaptors/adaptors.fasta


# 1b.
# BWA aln program - align reads against all wheat CSS chr arm contigs (RUN TIME: up to 80 minutes each sample; aln_sa_[12].sai=5.25GB/sample)
# Using aln command with each read file --> sampe to find the uniquely mapping read pairs
# Options:
# aln -n flag - number of mismatches you will allow between a read and a potential 
#				mapping location for that location to be considered a match.
#				The default is 4% of the read length.


# 2.
# BWA sampe command to get the sam output (RUN TIME: up to 6 hours; size of sampe.sam=39GB/sample ).
### Possible commands to count unique reads:
### -o 	- don't understand this option
### -n 1 - if set to 1, the XA flag will appear and will indicate a unique read pair (?)
### -N 1 - for discordant reads - it's using the same XA flag!
### Possible sam flags to count unique reads:	- these are related to single mates not read pairs 
### X0 Number of best hits
### X1 Number of suboptimal hits found by BWA 
### XT Type: Unique/Repeat/N/Mate-sw 	- Unique = unique alignment; Mate-sw = pair end rescue
### XA Alternative hits; format: (chr,pos,CIGAR,NM;)*

### XO should match X0:i:1	
### X1 should match X1:i:0	- NB: XO can be X0:i:1 and X1 e.g. X1:i:5 but XT is still XT:A:U - so can't use XT flag
### Decided to use X0:i:1 only as indicating a unique read


# Removing duplicates with Picard (Picard includes info on the levels of optical duplicates).
# 3.
# Converting to bam format, then sort by coord (RUN TIME: up to 3.75 hours; size of sampe_sort.bam=5-8GB per sample):	
		

# 4.
# Then removing duplicates (RUN TIME: up to 19.5 hours; 2.8 GB per sample)).
# OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 (default) - counts the number of duplicates within a similar area of the flow cell (ie that have replicated)


# 5.
# Count the number of non-duplicate reads produced above (RUN TIME: 3 minutes).
# NB - counting the output reveals odd numbers so Picard is not neccessarily removing read pairs, just ones that are duplicates


# 6.
# Getting properly paired reads using Samtools -f2 flag (RUN TIME: 1h --> 2 h; 2 GB mem).


# 7. Count the number of properly paired reads produced by samtools -f2 (RUN TIME: 8 minutes).


# 8.
# The sampe_sort_markdup_rm_filter_-f2.bam file should contain only properly paired reads.
# Reads with INT=2 may also map as pairs but each mate has the same orientation.
# Checked to see whether we have any of these reads: sam flags 67_131_115_179 and secondary (256).
# Result: there don't appear to be any read pairs both mapped in the same direction!
# NB - this step doesn't work in the compound command so running separately.


# 9.
# Now getting the on-target reads with Samtools (cmd from Dharanya) (RUN TIME: 35 minutes; 11GB mem).
# -L = outputs alignments overlapping the bed file.
# Note: the output only contains on-target reads. So if one mate isn't on-target,
# it is not included in the output.
paddingSize=0	# 0 2 100 150
#bedfile=/tgac/workarea/group-cg/baileyp/WheatLoLa/33.ExonCapture/all_wheat_singlectg.bedfile_plus_minus_${paddingSize}_1stbase_adjusted
#bedfile_IWGSP1_EnsemblPlnts22=/tgac/workarea/group-cg/baileyp/WheatLoLa/33.ExonCapture/indices/ensembl_plants_ref/Triticum_aestivum.IWGSP1.22.bedfile_plus_minus_0_1stbase_adjusted
#bedfile_JoseRef=/tgac/workarea/group-cg/baileyp/WheatLoLa/33.ExonCapture/indices/Jose_ref_July2014/wheat_genomeABD_update.UPDATES.bedfile_plus_minus_0_1stbase_adjusted


# 10.
# Count the number of on-target reads (RUN TIME: 14 minutes):
# The ontarget bam file may contain single reads or both read pairs on target but we need to know 
# the number of on targets as read pairs i.e. the number of read pairs where one or both mates are counted once.


# 11.
# To reduce the amount of data entering MAPS, removing reads with a mapping quality score of <20 (RUN TIME: 20 - 30 minutes).
# NB - also ran MAPS with reads with quality <20 removed from the bam file beforehand 
# and it makes no difference to the number of SNPs called - good - but using this file 
# doesn't seem to make the MAPS mpileup step quicker! It may start to with more samples.


# 12.
# Count the number of on-target reads with a quality score of >=20 (RUN TIME: 11 -24 minutes):
# The above ontarget bam file may contain single reads or both read pairs on target but we need to know 
# the number of on targets as read pairs i.e. the number of read pairs where one or both mates are counted once.


# 13.
# Using Bedtools coverageBed to estimate the coverage for the on target regions (cmd from Dharanya) (RUN TIME: 30 minutes for unique; 1 hour for concordant).
# 	NB - I think that the "coverageBed" command is the same as "Bedtools coverage".

# NB - need to have another bed file for bedtools where the coordinates are set +1 --> 0
# DS: Need to give your sorted, uniquely mapped, properly paired BAM file (whole one not the on-target one)
# and give your targets BED file which will give coverage at each base.
# The reason for using this BAM file is that you actually want to be able to see 
# the coverage of both the off-target and on-target reads when you apply the padding to the bedfile.

# Also, for the padding of the exons, it doesn't matter that the coordinates in the bedfile go beyond 
# the size of the reference contig at the 3' end - coverageBed and Samtools -L don't know 
# what these coordinates are - noting this here because this puzzled me for a while.
# (23.10.2014 - what about the sam header - doesn't that contain info on the length of each contig?)

# Reason for not using mpileup is where regions that have 0 coverage 
# samtools doesn't report them.

# -d is used to give depth at each base - this option gives this output:
# Cols:
# Overlapping	coords in bed	pos		coverage
# feature		file
#1AL_1027303     3986    4057    1       29
#1AL_1027303     3986    4057    2       30
#1AL_1027303     3986    4057    3       31
#1AL_1027303     3986    4057    4       33
#1AL_1027303     3986    4057    5       34


# 14. 
#  Dharanya's command to get the average coverage count


# 15.
# Get coverage at the full range of depths, >1x, >5x etc.


# 16.
# Line count on the bedfile - sampe_sort_markdup_rm_filter_-f2_covBed_orig_ref -
# to find the total number of on-target bases so the coverage values can be presented 
# as a % of the total.
# The bedfile for each sample should contain the same number of rows!
# Now counting on the fly instead of producing the output (RUN TIME: 8 hours for coverageBed and genomecov)


# 17.
# using Bedtools genomecov to look at whole genome coverage of the reads for each sample (RUN TIME: 7 hours; memory = 8 GB)
#bedT_GenomeFile=/tgac/workarea/group-cg/baileyp/WheatLoLa/33.ExonCapture/indices/IWGSC_CSS_all_scaff_v1_bedT.genome
#bedT_GenomeFileName=`basename $bedT_GenomeFile`
# 7.6.2016 - need to update this genome file to inlcude new 3B and ChrU (It is sorted w.r.t. contig size:
bedT_GenomeFile=/tgac/workarea/group-cg/baileyp/IWGSC_v2_ChrU_Ref/IWGSC_v2_ChrU_lengths_sort_-n_reorder_contigId_tab_length
bedT_GenomeFileName=`basename $bedT_GenomeFile`


# 18. Parse genomecov output for coverage values (RUN TIME: 7 hours):


# 19. 
# Count number of bases in the genome file RUN TIME: instant): 
# Use this value to find % of bases with >=1x, >=6x, >10=x coverage (Removed >19 
# and >29x coverge counts - taking too long)
# Number of bases in the current reference should be 10138701012
numbrBasesInIWGSC_v1=10,393,511,112		# Total number of bases in IWGSC_v2_ChrU.fa; number of bases in IWGSC_CSS_all_scaff_v1.fa = 10138701012


# 20.
# Investigating the coverage, specifically for 3 genes: (RUN TIME: 15 mins, mem=4GB)
# RUBISCO on 5AL, 5BL and 5DL and MDH
# Also calculating the off-target coverage using these genome off-target bed files (RUN TIME: ? mins, mem=17GB - seemed to need this amount of max mem!) :
off_targBedFile=/tgac/workarea/group-cg/baileyp/WheatLoLa/33.ExonCapture/indices/IWGSC_CSS_all_scaff_v1_off-target_regions_singlectg.bed
off_targBedFile_IWGSP1_EnsemblPlnts22=/tgac/workarea/group-cg/baileyp/WheatLoLa/33.ExonCapture/indices/IWGSC_CSS_all_scaff_v1_off-target_regions_IWGSP1.22.bed
off_targBedFile_JoseRef=/tgac/workarea/group-cg/baileyp/WheatLoLa/33.ExonCapture/indices/IWGSC_CSS_all_scaff_v1_off-target_regions_Jose_ref.bed


# 21. 
# Remove these files at the end:
#aln_sa_[12].sai
#sampe.sam
#sampe_sort.bam
#sampe_sort_markdup_rm.bam
#sampe_sort_markdup_rm_filter_-f2_ontarget_orig_ref_-q20.bam
#sampe_sort_markdup_rm_filter_-f2_covBed_orig_ref - no longer producing this file
#sampe_sort_markdup_rm_filter_-f2_bedT_genomecov_-d - no longer producing this file 

#####################################################################


# Pipeline starts here:
#:<<***string_marker_for_commenting_out_code***

# Use this array to redo specific samples - need to specify 
# the full path to each sample directory:
#samplesToRedo=( \
#'/tgac/data/reads/wheatsLoLaexomeCadenza/141114_SN7001150_0303_AC59LVACXX/Sample_1147_LIB10701_LDI8903' \
#$readsPath/wheatsLoLaexomeCadenza/141118_D00507_0043_BC5RJUACXX/Sample_1149_LIB10703_LDI8905 \
#$readsPath/wheatsLoLaexomeCadenza/141219_D00507_0056_AC5Y4KACXX/Sample_1238_LIB10705_LDI8907 \
#$readsPath/wheatsLoLaexomeCadenza/141219_D00507_0056_AC5Y4KACXX/Sample_1238_LIB10706_LDI8908 \
#$readsPath/wheatsLoLaexomeCadenza/141016_D00507_0038_BHALDRADXX/Sample_1116_LIB10020_LDI8275 \
#$readsPath/wheatsLoLaexomeCadenza/141028_D00534_0027_AC5L3NANXX/Sample_1128_LIB10021_LDI8276 \
#$readsPath/wheatsLoLaexomeCadenza/141031_D00507_0041_AC49M1ACXX/Sample_1134_LIB10495_LDI8690 \
#$readsPath/wheatsLoLaexomeCadenza/141029_SN790_0391_Ac4960acxx/Sample_1132_LIB10313_LDI8531 \
#$readsPath/wheatsLoLaexomeCadenza/141029_SN790_0391_Ac4960acxx/Sample_1132_LIB10314_LDI8532 \
#$readsPath/wheatsLoLaexomeCadenza/141029_SN790_0392_BC499TACXX/Sample_1133_LIB10317_LDI8535 \
#$readsPath/wheatsLoLaexomeCadenza/141029_SN790_0392_BC499TACXX/Sample_1133_LIB10318_LDI8536 \
#$readsPath/wheatsLoLaexomeCadenza/141203_SN319_0388_AC5PPPACXX/Sample_1162_LIB11008_LDI9085 \
#$readsPath/wheatsLoLaexomeCadenza/141203_SN319_0388_AC5PPPACXX/Sample_1162_LIB11009_LDI9086 \
#$readsPath/wheatsLoLaexomeCadenza/141209_SN7001150_0306_AC5RJ1ACXX/Sample_1171_LIB11494_LDI9433 \
#$readsPath/wheatsLoLaexomeCadenza/141209_SN7001150_0306_AC5RJ1ACXX/Sample_1171_LIB11495_LDI9434)

#plateToDo=35.ExonCapture_PlateG		# string to search in field 5 of file - to redo individual libraries, the string is e.g. 'LIB10667 \|LIB10668 \|LIB10669 \|LIB10671 \|LIB10703 \|LIB10672 \|LIB10673 \|LIB10674 \|LIB10675 \|LIB10676 \|LIB10677 \|LIB10679 \|LIB10704 \|LIB10680 \|LIB10681 \|LIB10682 \|LIB10683 \|LIB10688 \|LIB10689 \|LIB10690 \|LIB10691 \|LIB10692 \|LIB10693 \|LIB10695 \|LIB10702 \|LIB10696 \|LIB10698 \|LIB10699 \|LIB10705 \|LIB10628 \|LIB10630 \|LIB10687 \|LIB10706' 
#plateToDo=1143_LIB104
#plateToDo='LIB10604\|LIB10605\|LIB10606\|LIB10608\|LIB10609\|LIB10610\|LIB10611\|LIB10616\|LIB10617\|LIB10618\|LIB10619\|LIB10620\|LIB10621\|LIB10622\|LIB10623\|LIB10624\|LIB10625\|LIB10626\|LIB10627\|LIB10632\|LIB10633\|LIB10634\|LIB10635\|LIB10636\|LIB10638\|LIB10701'
plateToDo='PlateB'
#plateToDo='LIB10607\|LIB10614\|LIB10615\|LIB10612\|LIB10637\|LIB10656\|LIB10666\|LIB10679\|LIB10704\|LIB10683'
#plateToDo='LIB11398\|LIB11399'	#PlateG: 'LIB11398\|LIB11399' #PlateF: 'LIB10912\|LIB10913'	# PlateE: 'LIB10289\|LIB10290' #PlateD'LIB10403\|LIB10404' #PlateC: 'LIB11215\|LIB11216' # PlateB: 'LIB9924\|LIB9925' #PlateA: 'LIB10604\|LIB10605'


tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_300_Samples/35.ExonCapture_PreparingTable.txt | grep $plateToDo |
while read line; do

	pathsToLibs=`echo $line | cut -d ' ' -f 5`
	CadenzaId=`echo $line | cut -d ' ' -f 6`
	read1FileName=`echo $line | cut -d ' ' -f 7`
	pathToReads=`echo $line | cut -d ' ' -f 8`
	
	echo $pathsToLibs
	
	sampleDir=`basename $pathsToLibs`

	# Only make directory for the first run: 
	if [[ ! -d $sampleDir ]]; then mkdir $sampleDir; fi
		
	cd $sampleDir
	echo $sampleDir
		
	# Small modification for alpha design samples:
	#cd ${sampleDir}_v1
	#echo ${sampleDir}_v1
				
	#echo `ls $readsPath/$dir/${sampleDir}/*R1.fastq.gz`
	#fileR1=`ls $readsPath/$dir/${sampleDir}/*R1.fastq.gz`
	echo `ls ${pathToReads}${read1FileName}`
	fileR1=`ls ${pathToReads}${read1FileName}`

	#echo `ls $readsPath/$dir/${sampleDir}/*R2.fastq.gz`
	#fileR2=`ls $readsPath/$dir/${sampleDir}/*R2.fastq.gz`
	echo `ls ${pathToReads}${read1FileName} | sed s/_R1.fastq.gz/_R2.fastq.gz/`
	fileR2=`ls ${pathToReads}${read1FileName} | sed s/_R1.fastq.gz/_R2.fastq.gz/`

	# For bsub/qsub -J/-N cmd:
	jobId=`echo $sampleDir | awk -F '_' '{print $3}'`
	echo $jobId

	# To run via qsub:
	# 1. comment out the bsub line and add in echo "cd $PWD; to the line below:
	#		echo "cd $PWD; \
	# 2. remove the final quote then add in the qsub line at the bottom.
	# 3. cd to /uv2000-scratch/ and everything in the script should work the same.
		 		
	# If Tim produces exe's for uv2k2, then I would need to switch in the sourcing
	# or full path to the exe's here!
	source_testing_env='source /tgac/software/production/bin/testing_env'			# For uv2k2: source /tgac/software/production/bin/testing_env;		for hpc: source testing_env 
	source_production_env='source /tgac/software/production/bin/production_env'		# For uv2k2: source /tgac/software/production/bin/production_env;	for hpc: source production_env
	# sourceBWA='' 
	# sourceSamtools='' 
	# sourcePicard='' 
		  
	# Notes:
	# Wanted to time the pipeline. Problem with using "time" cmd is that it will only give the time for the first step!
	# Alternative to using compound commands would be to use bsub -K followed by wait - see notes.
	# 		Q: would the shell script hang? - think so - then I would need to bsub a wrapper script.  
	# Another alternative to running all steps consectutively is to put all the commands in a shell script, then bsub the whole script.
	#		This avoids issues with escaping various characters e.g. $ \` and " - can only specify rusage and -n once which is what I'm doing with the compound command anyway. 
	# Yet another way (Gonzalez): BWA ... && echo success > logfile - if program is successful, success will be echoed to logfile and script won't progress to next step until finished.
	# Yet another way is: Matt C. said to echo the command by piping it into a bsub (seems a bit like qsub) : echo "cmd" | bsub
	# Which way is best?
	# NB - 5.1.2015 - Could ||elize the coverage counting steps but then then I might loose my place in the queue, although they are very small jobs - testing it out...
	# 25.2.2015 - Jon told me about the bsub -w flag which looks like the best option - waits for the specified job name or jobId to finish before starting


#:<<***string_marker_for_commenting_out_code***
	cpu=1			# Set to 4 to fit into a UV2 nodes that have multiples of 4 cpus - I think I used 8 - see other 35. scripts; only need to set to 1 cpu for Tyson's script 
	mem_mb=15000	# Required for bsub; for Picard MarkDuplicates.jar, require 59000; for processing the mm reads (masking=18000 and Tyson approach=16000); mapping and sorting steps need 32 GB
	mem_gb=15		# Required for  Picard AND qsub; for Picard MarkDuplicates.jar, require 59; 18 for processing the mm reads (masking and Tyson approaches)
			
	# NB - to run via qsub, the bsub line can't be inside any "" because it itself has a double quote at the end; to run via bsub, just copy this line back to below 'echo "cd $PWD '.
	#bsub -J $jobId -n $cpu -q Prod128 -R "rusage [mem=${mem_mb}]" -oo sampe_sort_markdup_rm_filter_-f2.log "   OR   bsub -J $jobId -n $cpu -q Prod128 -R "rusage [mem=${mem_mb}]" -oo mm_mapping.log "
	# Normal settings for bsub: -n=15 and -R=60000! - 21.3.2015 - reduced resources to -n 9 rather than 15 and mem=53000 (to run on one host only: -R "span[hosts=1]" ) 
	#echo "cd $PWD
	#export TMPDIR=uv_tmp
	
	export readsPath \
	indexFile \
	PRIOR \
	QUAL_THRESH \
	ADAPTERS \
	paddingSize \
	bedfile \
	bedT_GenomeFile \
	bedT_GenomeFileName \
	numbrBasesInIWGSC_v1 \
	off_targBedFile \
	off_targBedFile_IWGSP1_EnsemblPlnts22 \
	off_targBedFile_JoseRef \
	source_testing_env \
	source_production_env \
	cpu \
	mem_mb \
	mem_gb \
	maxFileHandles \
	fileR1 \
	fileR2 \
	CadenzaId
	
	sbatch -J $jobId -c 1 --mem=14000 -p tgac-long -t 1-0:00   -o count_ontarget_reads.log   -e count_ontarget_reads.log_err   39.ExonCapture_coverage_slurm.sh
	### 1.2.2016 - BEFORE SUBMITTING TO GITHUB, MODIFY THIS SCRIPT WITH THE extra CODE FROM 39.ExonCapture_slurm.sh 
	#bsub -J $jobId -n $cpu -q normal -R "rusage [mem=${mem_mb}]" -oo bwa_aln_before_dup_rm.log "39.ExonCapture_slurm.sh"
	#bsub -J $jobId -n $cpu -q Prod256 -R "rusage [mem=${mem_mb}]" -oo count_ontarget_reads.log "
	
#	set +o posix
#$source_production_env; source sickle-1.2.2013
#	time sickle pe -t sanger -q $QUAL_THRESH \
#	-f <(/usr/users/ga002/baileyp/ProgramFiles/seqqs-master/seqqs -e -p raw_R1 $fileR1 | /usr/users/ga002/baileyp/ProgramFiles/scythe-master/scythe -a $ADAPTERS -p $PRIOR - 2> R1_scythe.stderr) \
#	-r <(/usr/users/ga002/baileyp/ProgramFiles/seqqs-master/seqqs -e -p raw_R2 $fileR2 | /usr/users/ga002/baileyp/ProgramFiles/scythe-master/scythe -a $ADAPTERS -p $PRIOR - 2> R2_scythe.stderr) \
#	-o >(/usr/users/ga002/baileyp/ProgramFiles/seqqs-master/seqqs -e -p trimmed_R1 - | gzip > R1_trimmed.fq.gz) \
#	-p >(/usr/users/ga002/baileyp/ProgramFiles/seqqs-master/seqqs -e -p trimmed_R2 - | gzip > R2_trimmed.fq.gz) \
#	-s >(/usr/users/ga002/baileyp/ProgramFiles/seqqs-master/seqqs -e -p singles - | gzip > singles_trimmed.fq.gz) 2> sickle.stderr
#$source_testing_env; source bwa-0.7.7
#	gzip -dc R1_trimmed.fq.gz | \
#	bwa aln -t $cpu $indexFile \
#	/dev/fd/0 \
#	> aln_sa_1.sai
#	gzip -dc R2_trimmed.fq.gz | \
#	bwa aln -t $cpu $indexFile \
#	/dev/fd/0 \
#	> aln_sa_2.sai
#bwa sampe -n 10 -N 0 $indexFile \
#	aln_sa_1.sai aln_sa_2.sai \
#	R1_trimmed.fq.gz R2_trimmed.fq.gz \
#	> sampe.sam
#$source_production_env; source samtools-0.1.18
#	samtools view -bS sampe.sam \
#	| samtools sort -m 10000000000   /dev/fd/0   sampe_sort
#$source_production_env; source jre-7.11; source picardtools-1.84
#	mkdir tmp
#	java -jar -Xmx${mem_gb}g -XX:ParallelGCThreads=$cpu -Djava.io.tmpdir=tmp /tgac/software/production/picardtools/1.84/x86_64/bin/MarkDuplicates.jar \
#	INPUT=sampe_sort.bam \
#	OUTPUT=sampe_sort_markdup_rm.bam \
#	METRICS_FILE=sampe_sort_markdup_metrics \
#	REMOVE_DUPLICATES=true \
#	AS=TRUE VALIDATION_STRINGENCY=LENIENT \
#	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 \
#	TMP_DIR=tmp
#samtools view -c   sampe_sort_markdup_rm.bam > sampe_sort_markdup_rm_count.log
#samtools index sampe_sort_markdup_rm.bam 
#samtools view -f2  sampe_sort_markdup_rm.bam   -bo sampe_sort_markdup_rm_filter_-f2.bam
#samtools view -c   sampe_sort_markdup_rm_filter_-f2.bam > sampe_sort_markdup_rm_filter_-f2_count.log
#samtools view -L $bedfile   sampe_sort_markdup_rm_filter_-f2.bam    -bo sampe_sort_markdup_rm_filter_-f2_ontarget_orig_ref.bam
#samtools view      sampe_sort_markdup_rm_filter_-f2_ontarget_orig_ref.bam \
#	| awk '{print \$1}' | sort -u | wc -l \
#	> sampe_sort_markdup_rm_filter_-f2_ontarget_orig_ref_count_read_pairs.log
#samtools view -q20 sampe_sort_markdup_rm_filter_-f2_ontarget_orig_ref.bam   -bo sampe_sort_markdup_rm_filter_-f2_ontarget_orig_ref_-q20.bam sampe_sort_markdup_rm_filter_-f2_covBed_orig_ref
#samtools view sampe_sort_markdup_rm_filter_-f2_ontarget_orig_ref_-q20.bam \
#	| awk '{print \$1}' | sort -u | wc -l \
#	> sampe_sort_markdup_rm_filter_-f2_ontarget_orig_ref_-q20_count.log


# 22.1.2016 - finalising the on-target coverage stats (RUN TIME: up to 4.5h; 10 GB mem)
#bedfile=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/sequence_space_coverage/blast/exons_0-based.bed
#bedfile=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/sequence_space_coverage/blast/exons_0-based_200bp_padding.bed
#bedfile=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/sequence_space_coverage/blast/bedfile_to_get_max_cov.bed
#bedfile=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/sequence_space_coverage/blast/exons_0-based_200bp_padding_sort_-k11_-k22n_bt_merge_comm_-23uniq.bed


#fraction=0.2	# also try a 1bp overlap by missing out this flag; NB -wa doesn't have to be specified to get the reads to be outputted
#$source_production_env; source bedtools-2.17.0
#bedtools intersect -abam *sampe_sort_markdup_rm.bam \
#-b \$bedfile \
#-f 0.000000001 \
#-u \
#> sampe_sort_markdup_rm_on-target_bedt_-f1bp.bam
#$source_production_env; source samtools-0.1.18
#samtools view -c   sampe_sort_markdup_rm_on-target_bedt_-f1bp.bam > sampe_sort_markdup_rm_on-target_bedt_-f1bp_count.log
#samtools view -L \$bedfile   *sampe_sort_markdup_rm.bam    -bo sampe_sort_markdup_rm_samtools_-L_ontarget.bam
#samtools view -c   sampe_sort_markdup_rm_samtools_-L_ontarget.bam > sampe_sort_markdup_rm_samtools_-L_ontarget_count.log


# Assessing coverage of contigs with no mutations with this bedfile: exons_0-based_200bp_padding_sort_-k11_-k22n_bt_merge_comm_-23uniq.bed (only present in in the exon capture - no mutations)
#coverageBed -abam *sampe_sort_markdup_rm.bam -b \$bedfile -d \
#	| ~baileyp/bin/35.getCoverageStats_On-targetRegions.pl \
#	> sampe_sort_markdup_rm_covBed_more_1_6_10_15_20_25_100x_covs_comm_-23uniq.log


#bedfile=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/sequence_space_coverage/blast/exons_0-based_200bp_padding_sort_-k11_-k22n.bed (I think I should have used the bt_merge file but it shoudln't make any difference)
#bedfile=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/sequence_space_coverage/blast/exons_0-based_200bp_padding_sort_-k11_-k22n_bt_merge.bed
#coverageBed -abam *sampe_sort_markdup_rm.bam -b \$bedfile -d \
#	| ~baileyp/bin/35.getCoverageStats_On-targetRegions.pl \
#	> sampe_sort_markdup_rm_covBed_more_1_6_10_15_20_25_100x_covs.log
#" # End of bsub

#coverageBed -abam sampe_sort_markdup_rm_filter_-f2.bam -b ${bedfile_IWGSP1_EnsemblPlnts22} -d \
#	| ~baileyp/bin/35.getCoverageStats_On-targetRegions.pl \
#	> sampe_sort_markdup_Summary Stats for The Exon Capturerm_filter_-f2_covBed_IWGSP1_EnsemblPlnts22_ref_more_1_6_10_15_20_25_100x_covs.log
#coverageBed -abam sampe_sort_markdup_rm_filter_-f2.bam -b $bedfile_JoseRef -d \
#	| ~baileyp/bin/35.getCoverageStats_On-targetRegions.pl \
#	> sampe_sort_markdup_rm_filter_-f2_covBed_Jose_ref_more_1_6_10_15_20_25_100x_covs.log
#bedtools genomecov -d -ibam sampe_sort_markdup_rm_filter_-f2.bam -g $bedT_GenomeFile \
#	| ~baileyp/bin/35.getCoverageStats_WholeGenome.pl $numbrBasesInIWGSC_v1 \
#	> sampe_sort_markdup_rm_filter_-f2_bedT_genomecov_-d_more_1_6_10_20x_covs.log
#bedtools coverage -d -abam sampe_sort_markdup_rm_filter_-f2_ontarget_orig_ref.bam \
#	-b /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_300_Samples/genes_for_qPCR_zero_based.txt \
#	| ~baileyp/bin/35.getCoverageStats_On-targetRegions.pl \
#	> sampe_sort_markdup_rm_filter_-f2_bedT_coverage_-d_for_qPCR.log			
#coverageBed -abam sampe_sort_markdup_rm_filter_-f2.bam -b ${off_targBedFile} -d \
#	| ~baileyp/bin/35.getCoverageStats_On-targetRegions.pl \
#	> sampe_sort_markdup_rm_filter_-f2_covBed_off-targ_orig_ref_more_1_6_10_15_20_25_100x_covs.log
#coverageBed -abam sampe_sort_markdup_rm_filter_-f2.bam -b ${off_targBedFile_IWGSP1_EnsemblPlnts22} -d \
#	| ~baileyp/bin/35.getCoverageStats_On-targetRegions.pl \
#	> sampe_sort_markdup_rm_filter_-f2_covBed_off-targ_IWGSP1_EnsemblPlnts22_ref_more_1_6_10_15_20_25_100x_covs.log			
#coverageBed -abam sampe_sort_markdup_rm_filter_-f2.bam -b ${off_targBedFile_JoseRef} -d \
#	| ~baileyp/bin/35.getCoverageStats_On-targetRegions.pl \
#	> sampe_sort_markdup_rm_filter_-f2_covBed_off-targ_Jose_ref_more_1_6_10_15_20_25_100x_covs.log
#rm aln_sa_[12].sai sampe.sam
# sampe_sort.bam  sampe_sort_markdup_rm_filter_-f2.bam  sampe_sort_markdup_rm_filter_-f2_ontarget_orig_ref.bam
#rm -fR tmp
#rm singles_trimmed.fq.gz

# Additional code to obtain unmapped reads from the sampe_sort_markdup_rm.bam file (15 samples only). 
#	$source_production_env; source samtools-0.1.18
	# Extract unmapped read whose mate is mapped:
#	samtools view -f 4 -F 264  -bo sampe_sort_markdup_rm_tmp1.bam  sampe_sort_markdup_rm.bam
	# Extract mapped read whose mate is unmapped:
#	samtools view -f 8 -F 260  -bo sampe_sort_markdup_rm_tmp2.bam  sampe_sort_markdup_rm.bam
	# Extract paired reads that are both unmapped:
#	samtools view -f 12 -F 256  -bo sampe_sort_markdup_rm_tmp3.bam  sampe_sort_markdup_rm.bam
	# Merge the 3 tmp files:
#	samtools merge -u  sampe_sort_markdup_rm_unmap_rds_merg.bam  sampe_sort_markdup_rm_tmp1.bam  sampe_sort_markdup_rm_tmp2.bam  sampe_sort_markdup_rm_tmp3.bam
#	samtools sort -m 10000000000 -n  sampe_sort_markdup_rm_unmap_rds_merg.bam  sampe_sort_markdup_rm_unmap_rds_merg_sort-n
#$source_production_env; source bedtools-2.17.0
#	bedtools bamtofastq -i sampe_sort_markdup_rm_unmap_rds_merg_sort-n.bam \
#   -fq ../unmapped_reads/${jobId}_unmapped_R1.fastq \
#   -fq2 ../unmapped_reads/${jobId}_unmapped_R2.fastq
    #rm sampe_sort_markdup_rm_tmp[123].bam  sampe_sort_markdup_rm_unmap_rds_merg.bam  sampe_sort_markdup_rm_unmap_rds_merg_sort-n.bam
    


	# Analyse all the reads that map ambiguously and prepare to remap and run MAPS using a masked reference.
	# Run masking script, contig-oriented approach, in non-masking mode to get the mm read Id's:
#$source_production_env; source samtools-0.1.18
#	samtools view *sampe_sort_markdup_rm.bam \
#	| ~baileyp/bin/39.WheatReferenceMasking_BJsWay.pl \
#	-r mm_reads_list \
#	-m no \
#	> 39.WheatReferenceMasking.log 2>&1

	# Sort unique the read list so each read is extracted exactly once from the fastq files:
#	sort -u mm_reads_list > mm_reads_list_sort_-u

	# Get the fastq records for each read:
#set +o posix
	# Trimmed: R1_trimmed.fq.gz OR untrimmed: $fileR1
#	~baileyp/ProgramFiles/seqtk/seqtk-master/seqtk \
#   subseq <(zcat R1_trimmed.fq.gz) mm_reads_list_sort_-u > mm_reads_list_R1.fastq
#	~baileyp/ProgramFiles/seqtk/seqtk-master/seqtk \
#	subseq <(zcat R2_trimmed.fq.gz) mm_reads_list_sort_-u > mm_reads_list_R2.fastq
	
	# The master masked reference is prepared once for 10 samples with 39.K-masking.sh, then used in the alignment steps below.

#	refVersion=10bams	# use '\' to access the contents
#	mm_indexFile=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_masked_ref/39.ExonCapture_masked_ref_10_CadWT_samples_v1/masked_contigs_bwa_index
	# Align to BWA and count (RUN TIME: up to 4h 15 mins; up to 18 GB mem)
#$source_testing_env; source bwa-0.7.7
#   bwa aln -t $cpu \$mm_indexFile \
#	mm_reads_list_R1.fastq \
#	> mm_aln_sa_1_\$refVersion.sai
#	bwa aln -t $cpu \$mm_indexFile \
#	mm_reads_list_R2.fastq \
#	> mm_aln_sa_2_\$refVersion.sai
#	bwa sampe -n 10 -N 0 \$mm_indexFile \
#	mm_aln_sa_1_\$refVersion.sai mm_aln_sa_2_\$refVersion.sai \
#	mm_reads_list_R1.fastq mm_reads_list_R2.fastq \
#	> mm_sampe_\$refVersion.sam
#$source_production_env; source samtools-0.1.18
#	samtools view -bS mm_sampe_\$refVersion.sam \
#	| samtools sort -m 10000000000   /dev/fd/0   mm_sampe_sort_\$refVersion

	# Stats new and original bam files:
	# New bam file: 
	# Number of reads:
#	samtools view -c  mm_sampe_sort_\$refVersion.bam > mm_sampe_sort_counts_\$refVersion.log
	# Number of reads with multiple hits:
#	samtools view mm_sampe_sort_\$refVersion.bam | grep XA | wc -l > mm_sampe_sort_XA_counts_\$refVersion.log
	# Number of reads with mapping score = 60:
#	samtools view  mm_sampe_sort_\$refVersion.bam | awk '\$5 ~ /60/' | wc -l > mm_sampe_sort_mapqual60_counts_\$refVersion.log
	# Number of reads with mapping score >= 20:
#	samtools view  mm_sampe_sort_\$refVersion.bam | awk '\$5 >= 20' | wc -l > mm_sampe_sort_mapqual20_counts_\$refVersion.log
	# Original bam file:
	# Number of reads:
#	samtools view -c  *sampe_sort_markdup_rm.bam > sampe_sort_counts.log
	# Number of reads with multiple hits:
#	samtools view  *sampe_sort_markdup_rm.bam | grep XA | wc -l > sampe_sort_XA_counts.log 
	# Number of reads with mapping score = 60:
#	samtools view  *sampe_sort_markdup_rm.bam | grep XA | awk '\$5 ~ /60/' | wc -l > sampe_sort_XA_mapqual60_counts.log
	# Number of reads with mapping score >= 20:
#	samtools view  *sampe_sort_markdup_rm.bam | grep XA | awk '\$5 >= 20' | wc -l > sampe_sort_XA_mapqual20_counts.log
	
	# Remove intermediate files:
#	rm   mm_reads_list   mm_reads_list_sort_-u   mm_aln_sa_[12]_\$refVersion.sai   mm_sampe_\$refVersion.sam   mm_reads_list_R[12].fastq

	### Now run MAPS on the mm_sampe_sort.bam files - see script 35.runMAPS_InChrArmChunks_mm.sh
	
	
	# Tyson's approach to rescue multiply mapped reads (RUN TIME: 45 minutes; 15 GB mem):
	# (Note the use of both the -A and -I flags)
#$source_production_env; source samtools-0.1.18
#	bamFileName=`ls *_sampe_sort_markdup_rm.bam | sed s/\.bam//`
#	set +o posix
#	source /tgac/software/testing/bin/python-2.7.5
#	python ~baileyp/bin/multi-map-corrector-V1.6.py \
#	-i <(samtools view \${bamFileName}.bam) \
#	-l /tgac/workarea/group-cg/baileyp/IWGSC_v2_ChrU_Ref/IWGSC_v2_ChrU_lengths_sort_-n \
#	-A \
#	-I \
#	> ty_mm_reads_v1.6.sam

	# Get header from original Cadenza*_sampe_sort_markdup_rm.bam
#	samtools view -H \${bamFileName}.bam > \${bamFileName}_header_only.sam

#	cat \${bamFileName}_header_only.sam ty_mm_reads_v1.6.sam > ty_mm_reads_v1.6_plus_headr.sam
	
	# Sort new sam by coord:
#	samtools view -bS ty_mm_reads_v1.6_plus_headr.sam \
#   | samtools sort -m 10000000000   /dev/fd/0   ty_mm_reads_v1.6_plus_headr_sort
	
    # Count number of reads: 
#    samtools view -c  ty_mm_reads_v1.6_plus_headr_sort.bam > ty_mm_reads_v1.6_plus_headr_sort_counts.log
    
#	rm ty_mm_reads_v1.6.sam   *_sampe_sort_markdup_rm_header_only.sam   ty_mm_reads_v1.6_plus_headr.sam
	
	### Can now run MAPS on the ty_mm_reads_sort_plus_headr.bam files - see script 35.runMAPS_InChrArmChunks_mm.sh and Tyson's compression version, 35.runMAPS_TY_Cadenza.sh
#	"
	#| qsub -N $jobId -q Prod -l select=1:ncpus=$cpu:mem=${mem_gb}gb -j eo
	
	# 9.7.2015 - Chris said that I can use:
	# qsub -N $jobId -q Prod -l select=1:ncpus=$cpu:mem=${mem_gb}gb -o /dev/null -e /dev/null
																     # NBNB - trying to name the qsub log file here - but it might not work!!!!!
																     # -oe sampe_sort_markdup_rm.log - didn't give an error but might not work
																     # for the moment using -j flag like i have done in qsub scripts 
#***string_marker_for_commenting_out_code***
														
	# Coming out of current directory to the one above, ready for next sample:
	cd ../
done
#***string_marker_for_commenting_out_code***




# 30.11.2014 - Reminder - counting genome bases in reference:
#bsub "awk '{sum +=\$2} END {print \"Number of bases in the reference:\" sum}' /tgac/workarea/group-cg/baileyp/WheatLoLa/33.ExonCapture/indices/IWGSC_CSS_all_scaff_v1_bedT.genome" 
# Gives 10138701012
# Tried to count each sampe_sort_markdup_rm_filter_-f2_bedT_genomecov_-d file but it gives a different number of lines per sample: 
#wc -l sampe_sort_markdup_rm_filter_-f2_bedT_genomecov_-d > sampe_sort_markdup_rm_filter_-f2_bedT_genomecov_-d_base_count.log
# It should give the same number of lines as the number of bases in the genomecov file - why doesn't it?
# Maybe it is because not all contigs are covered by ANY reads in the bam file so the contigs are not printed out - 11.1.2014 - this seems the only explanation




# Step 8 above to look for sam flags 67_131_115_179_256 failed and I don't know why.
# So checking here in a separate loop (RUN TIME: xxx minutes):
:<<***string_marker_for_commenting_out_code***	# Works except for the backticks apparently
for dir in ${outerReadsDir[@]}; do

	echo `ls -d $readsPath/$dir`

	for samplePath in `ls -d $readsPath/$dir/${samplePrefix}*/`; do
	
		sampleDir=`basename $samplePath`

		cd $sampleDir
		echo $sampleDir
		ls sampe_sort_markdup_rm_filter_-f2.bam

		bsub -J flag_check -q Prod128 -R "rusage[mem=3000]" -oo sampe_sort_markdup_rm_filter_-f2_check_for_67_131_115_179.log " \
		source samtools-0.1.18; \
		samtools view sampe_sort_markdup_rm_filter_-f2.bam | awk '\$2 ~ /^67|131|115|179|256/ {print \$0}' \
		"
		
		# Coming out of current directory to the one above, ready for next sample:
		cd ../
	done
done
***string_marker_for_commenting_out_code***
# Result: there don't appear to be any read pairs both mapped in the same direction!










#################################################################
# Mutation/SNP Calling Pipeline for UNIQUELY aligned reads pairs:
#################################################################

# Explanation of steps:

# 1.
# Finding unique and concordantly aligned read pairs from this file e.g. sampe_sort_markdup_rm_filter_-f2.bam (RUN TIME: up to 3 hours)
# i.e. all reads that are properly paired, on- and off-target. can't use the on-target file because not all 
# reads are in there as read pairs. 

# Then performing SNP calling on the UNIQUELY aligned reads.

# First sort by name so the read pairs are next to each other. (RUN TIME: up to 3 hours)
# Then extract uniquely aligned reads using a Perl oneliner (RUN TIME: up to xxxxxxxxxx hours)


#################################################################


# Pipeline starts here:
:<<***string_marker_for_commenting_out_code***		# Works except for the backticks apparently
for dir in ${outerReadsDir[@]}; do

	echo `ls -d $readsPath/$dir`

	for samplePath in `ls -d $readsPath/$dir/${samplePrefix}*/`; do

		sampleDir=`basename $samplePath`

		cd $sampleDir
		echo $sampleDir
		ls sampe_sort_markdup_rm_filter_-f2.bam

		
		# If Tim produces exe's for uv2k2, then I would need to switch in the sourcing
		# or full path to the exe's here!
		source_testing_env='source /tgac/software/production/bin/testing_env'			# For uv2k2: source /tgac/software/production/bin/testing_env;		for hpc: source testing_env 
		source_production_env='source /tgac/software/production/bin/production_env'		# For uv2k2: source /tgac/software/production/bin/production_env;	for hpc: source production_env
		# sourceBWA='' 
		# sourceSamtools='' 
		# sourcePicard='' 


#		bsub -J srt_by_nam -q Prod128 -R "rusage[mem=25000]" -oo sampe_sort_markdup_rm_filter_-f2_srt_by_nam.log " \
#
#$source_production_env;	source samtools-0.1.18; \
#		samtools sort -n -m 10000000000  sampe_sort_markdup_rm_filter_-f2.bam  sampe_sort_markdup_rm_filter_-f2_srt_by_nam \
#		"

		### UPTOHERE 13.11.2014
		### NB - think about zipping up the above bam file and piping it in - ditto for the ontarget bam file (?) - might have issues in MAPS
		


		# Coming out of current directory to the one above, ready for next sample:
		cd ../
	
	done
done
***string_marker_for_commenting_out_code***



# Now get all the unique and concordant read pairs (RUN TIME: 2 hours)
# then piping output straight into samtools to convert to bam:
### NB - HAVE COMBINED THIS STEP WITH THE STEP ABOVE
# NB - notice that the sampe_sort_markdup_rm_filter_-f2_srt_by_nam.bam file
# will be identical for both X0 and X0_X1 data but if you want to run them 
# together the outputs will role into the same file - not good. So creating 2
# identical files here. 
#flags=X0_X1	# X0  X0_X1 - also need to alter conditional
#for sampleDir in `ls -d Sample_1071_LIB92*/`; do

#	cd $sampleDir
#	echo $sampleDir
#	ls sampe_sort_markdup_rm_filter_-f2.bam

#	bsub -J get_${flags} -q Prod128 -R "rusage[mem=25000]" -oo sampe_sort_markdup_rm_filter_-f2_srt_by_nam_${flags}_sort_coord_recheck.log " \
#	source samtools-0.1.18; \
#	samtools sort -n -m 10000000000  sampe_sort_markdup_rm_filter_-f2.bam  sampe_sort_markdup_rm_filter_-f2_srt_by_nam_${flags}_recheck;
#	samtools view -h sampe_sort_markdup_rm_filter_-f2_srt_by_nam_${flags}_recheck.bam \
#	| perl -e '
#	\$lastLine = \"\";
#	\$LastReadName = \"\";
#	#$uniqReadPairCountr = 0;
#	#$numbrReadsInFile = 0;

#	while(\$line = <>)	{

#		if(\$line =~ m/^@/)	{# The header lines
#			print \$line;
#			next; 	
#		}

#		@fields = split /\t/, \$line;
### ***** YOU NEED TO ALTER THESE CONDITIONALS DEPENDING ON WHETHER UNIQUE READS ONLY ARE REQUIRED  *****
#		if( (\$fields[0] eq \$LastReadName) && (\$lastLine =~ m /X0:i:1\t/) && (\$line =~ m /X0:i:1\t/) )	{
		###if( (\$fields[0] eq \$LastReadName) && (\$lastLine =~ m/X0:i:1\t/) && (\$lastLine =~ m/X1:i:0\t/) && (\$line =~ m/X0:i:1\t/) && (\$line =~ m/X1:i:0\t/)  )	{

#			print \$lastLine, \$line;
			#$uniqReadPairCountr++;
#		}		
#		\$lastLine = \$line;
#		\$LastReadName = \$fields[0];
		#$numbrReadsInFile++;
#	}
#	' \
#	| samtools view -bS /dev/fd/0 \
#	| samtools sort -m 10000000000  /dev/fd/0  sampe_sort_markdup_rm_filter_-f2_srt_by_nam_${flags}_sort_coord_recheck \	CALL it REsort
#	"
	# Note - I am also piping the sam output plus header lines directly
	# back into bam format with samtools
	
	# Coming out of current directory to the one above, ready for next sample:
#	cd ../
#done
# Running perl this way via bsub, need to backslash these chars: "  $ 	any more?

# Then return to the samtools -f2 step above to count the properly paired reads
# then count the number of uniq+concordant alns
### NBNB - I should have sorted them by coordinate here as well I guess like so:
#	| samtools view -bS /dev/fd/0 \
#	| samtools sort -m 10000000000  /dev/fd/0  sampe_sort_L001_markdup_rm_filter_-f2_srt_by_nam_${flags}_sort_coord 
# Then extract the on-target reads - does indexing the bam help with this?
# Then count the on-target reads
# continue with pipeline and call SNPs on these files