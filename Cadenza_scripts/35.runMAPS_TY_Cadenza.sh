#!/bin/bash
#######################################################################
# 35.runMAPS_TY.sh		Modified from 35.runMAPS_InChrArmChunks_H-VG.sh		
#
# Paul Bailey	8.12.2015
#######################################################################


referenceFile=/tgac/workarea/group-cg/baileyp/IWGSC_v2_ChrU_Ref/IWGSC_v2_ChrU.fa					# The full and final reference for the project (IWGSC CSS minus CSS 3B + new 3BSeq + CSS 3B not in 3BSeq + ChrU)
#referenceFile=/tgac/workarea/group-cg/baileyp/references/wheat_genomes/my_iwgsc_v2/my_iwgsc_v2_split/3B_v443_contigs_only.fa	# had to redo MAPS just for the 3B v443 contigs

#MAPS-part1:
#minLib=23	# Now taken from the MAPSruns* files
minCov=10
hetOneMinPer=20
maxCov=10000			# default=2000

#MAPS-part2:
hetMinPer=10			# Running with -p=20 (original value I used), 15, 25 and 40
#homMinCov=3


cpu=1			# Used 28 cpu for Tyson's script 				
mem_mb=2000		# 86000 MB seems OK for all samples
mem_gb=2		# 86 GB

# Running MAPS on ty mm bam files (Total Run Time: 24-42h, max mem 85 GB):

# These MAPS runs have to be set off by hand:
### 18.2.2016 - actually I have added these now to the control file, useful for repeating some of the MAPS steps if the bam file links already exist.
# PlateH + 12 Plate20 samples:
#libsToDo='LIB15657\|LIB15660\|LIB15661\|LIB15662\|LIB15665\|LIB15666\|LIB15667\|LIB15670\|LIB15671\|LIB15672\|LIB15675\|LIB15676\|LIB15677\|LIB15678\|LIB15679\|LIB15680' # PlateH
#libsToDo='LIB17351\|LIB17352\|LIB17353\|LIB17354\|LIB17355\|LIB17356\|LIB17357\|LIB17358\|LIB17359\|LIB17360\|LIB17361\|LIB17362'	# Plate20 (these are in PlateN dir)

# PlateL + 5 Plate20 samples:
#libsToDo='LIB16232\|LIB16233\|LIB16234\|LIB16235\|LIB16236\|LIB16237\|LIB16238\|LIB16239\|LIB16240\|LIB16241\|LIB16242\|LIB16243\|LIB16244\|LIB16245\|LIB16246\|LIB16247'	# PlateL
#libsToDo='LIB15419\|LIB15424\|LIB15441\|LIB15443\|LIB15445'	# Plate20 (these are in PlateN dir)

# Poor captures:
# Plates A --> G:
# List of captures (On-targReadPairs < 5m) - adding these to a separate run of the 19 5m -10m samples - T = 27:
#libsToDo='LIB10607\|LIB10612\|LIB10614\|LIB10615\|LIB10637\|LIB10672\|LIB10231\|LIB10235\|LIB10241\|LIB11002\|LIB11484\|LIB11486\|LIB10943\|LIB10944\|LIB10945\|LIB10701' # Plate[ACEFG]
#libsToDo='LIB9944\|LIB9946\|LIB9947\|LIB9956\|LIB9957\|LIB9958\|LIB9959\|LIB10454\|LIB10460' # Plate[BD]

#List of captures (On-targReadPairs ~5m - ~10m):
#libsToDo='LIB10612\|LIB10614\|LIB10615\|LIB10637\|LIB10672\|LIB10231\|LIB10235\|LIB10241\|LIB11002\|LIB11484\|LIB11486\|LIB10943\|LIB10944\|LIB10945\|LIB10701' # Plate[ACEFG] 
#libsToDo='LIB9956\|LIB9957\|LIB9958\|LIB9959\|LIB10454\|LIB10460' # Plate[BD]

#List of captures (On-targReadPairs < ~10m - ~15m):
#libsToDo='LIB10674\|LIB10675\|LIB10676\|LIB10677\|LIB10679\|LIB10689\|LIB10695\|LIB10628\|LIB10687\|LIB10994\|LIB11457\|LIB11201\|LIB11227\|LIB11251\|LIB11255\|LIB11256\|LIB11268\|LIB11205\|LIB11213\|LIB11237\|LIB11242\|LIB11245' # Plate[ACEFG] 
#libsToDo='LIB9940\|LIB9941\|LIB9942\|LIB9943\|LIB10012\|LIB10013\|LIB10014\|LIB10015' # Plate[BD]

#List of captures (On-targReadPairs < ~15m - ~20m):
#libsToDo='LIB10664\|LIB10665\|LIB10666\|LIB10667\|LIB10668\|LIB10669\|LIB10671\|LIB10644\|LIB10645\|LIB10646\|LIB10647\|LIB10648\|LIB10649\|LIB10650\|LIB10651\|LIB10673\|LIB10943\|LIB10944\|LIB10945\|LIB10946\|LIB10701\|LIB10702\|LIB10705\|LIB10706' # Plates [ACEFG]:

# Poor captures - plates H to N:
#bamFileLoctn=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[JKLM]	# Plates [JKLM]:
#bamFileLoctn=/tgac/scratch/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[HN]	# Plates [HN]:
#List of captures (On-targReadPairs < 5 million) - adding these to a separate run of the 19 5-10m samples - T = 27:
#libsToDo='LIB15925\|LIB15926\|LIB15927\|LIB16227\|LIB16099\|LIB16148\|LIB16150\|LIB16158\|LIB15891\|LIB15895\|LIB16220\|LIB16569'	# Plate JKLM
#libsToDo='LIB15597\|LIB15658\|LIB15437\|LIB15640\|LIB15659\|LIB15664\|LIB15406\|LIB15407\|LIB15408\|LIB15428\|LIB15430\|LIB15435\|LIB15478\|LIB15417\|LIB15448'	# Plate H+N

#List of captures (On-targReadPairs  ~5m - ~10m):  
#libsToDo='LIB16099\|LIB16148\|LIB16150\|LIB16158\|LIB15891\|LIB15895\|LIB16220\|LIB16569' # PlateJKLM
#libsToDo='LIB15640\|LIB15659\|LIB15664\|LIB15406\|LIB15407\|LIB15408\|LIB15428\|LIB15430\|LIB15478\|LIB15417\|LIB15448' # PlateH+N





#fileList=/tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/54_test_samples/35.ExonCapture_PreparingTable_54_TestSamples.txt
fileList=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/VEP_analysis_aux_files/ExonCapture_PlateTABCDEFGHIJKLMN.txt

bamFileName=ty_mm_reads_v1.6_plus_headr_sort.bam

MAPSruns=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/MAPSruns_minLibInfo_ty_mm.txt
runsToIgnore='-v #NBNB'
#runsToDo='ty_mm_less5m_On-targReadPairs\|ty_mm_5m-10m_On-targReadPairs\|ty_mm_10m-15m_On-targReadPairs\|ty_mm_15m-20m_On-targReadPairs\|ty_mm_less_5m_On-targReadPairs\|ty_mm_5m-10m_On-targReadPairs\|ty_mm_MAPS_bams_to_use_LIB15657-LIB15680-12_Plate20_samples\|ty_mm_MAPS_bams_to_use_LIB16232-LIB16247_5_Plate20_samples'
#runsToDo='ty_mm_MAPS_bams_to_use_LIB16072-LIB16103\|ty_mm_MAPS_bams_to_use_LIB16104-LIB16135\|ty_mm_MAPS_bams_to_use_LIB16136-LIB16167\|ty_mm_MAPS_bams_to_use_LIB15880-LIB15951\|ty_mm_MAPS_bams_to_use_LIB15936-LIB15974\|ty_mm_MAPS_bams_to_use_LIB15904-LIB15935\|ty_mm_MAPS_bams_to_use_LIB16168-LIB16207\|ty_mm_MAPS_bams_to_use_LIB16208-LIB16231\|ty_mm_MAPS_bams_to_use_LIB16530-LIB16561\|ty_mm_MAPS_bams_to_use_LIB16562-LIB16593\|ty_mm_MAPS_bams_to_use_LIB16594-LIB16624\|ty_mm_MAPS_bams_to_use_LIB15385-LIB15439\|ty_mm_MAPS_bams_to_use_LIB15449-LIB15400\|54_test_samples/ty_mm_MAPS_bams_to_use_LIB5772-LIB8427\|54_test_samples/ty_mm_MAPS_bams_to_use_LIB8428-LIB9232\|ty_mm_less5m_On-targReadPairs\|ty_mm_5m-10m_On-targReadPairs\|ty_mm_10m-15m_On-targReadPairs\|ty_mm_15m-20m_On-targReadPairs\|ty_mm_less_5m_On-targReadPairs\|ty_mm_5m-10m_On-targReadPairs\|ty_mm_MAPS_bams_to_use_LIB15657-LIB15680-12_Plate20_samples\|ty_mm_MAPS_bams_to_use_LIB16232-LIB16247_5_Plate20_samples'
#runsToDo='ty_mm_MAPS_bams_to_use_LIB16208-LIB16231\|ty_mm_MAPS_bams_to_use_LIB16530-LIB16561\|ty_mm_MAPS_bams_to_use_LIB16562-LIB16593\|ty_mm_MAPS_bams_to_use_LIB16594-LIB16624\|ty_mm_MAPS_bams_to_use_LIB15385-LIB15439\|ty_mm_MAPS_bams_to_use_LIB15449-LIB15400'
#runsToDo='54_test_samples/ty_mm_MAPS_bams_to_use_LIB5772-LIB8427\|54_test_samples/ty_mm_MAPS_bams_to_use_LIB8428-LIB9232'
runsToDo='ty_mm_MAPS_bams_to_use_LIB10912-LIB10938\|ty_mm_MAPS_bams_to_use_LIB10947-LIB10970\|ty_mm_MAPS_bams_to_use_LIB10975-LIB11009\|ty_mm_MAPS_bams_to_use_LIB11394-LIB11491\|ty_mm_MAPS_bams_to_use_LIB11458-LIB11426\|ty_mm_MAPS_bams_to_use_LIB11426-LIB11493'	# PlateE

:<<***string_marker_for_commenting_out_code***
#tail -n +2 $MAPSruns | grep $runsToIgnore |
tail -n +2 $MAPSruns | grep $runsToDo |
while read line; do

	fullPathToMAPSrun=`echo $line | cut -d ' ' -f 1`
	minLib=`echo $line | cut -d ' ' -f 2`
	libsToDo=`echo $line | cut -d ' ' -f 3`

	bamFileLoctn=`dirname $fullPathToMAPSrun`
	# For creating symlinks by hand with the loop below:
	#bamFileLoctn=/tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[ACEFG]
	#bamFileLoctn=/tgac/scratch/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[BD]
	#bamFileLoctn=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[JKLM]
	#bamFileLoctn=/tgac/scratch/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[HN]
		
	# Only make the MAPS run directory if it doesn't exist: 
#	if [[ ! -d $fullPathToMAPSrun ]]; then mkdir $fullPathToMAPSrun; fi

	echo
	echo
	echo
	ls -d $fullPathToMAPSrun 
	ls -d $bamFileLoctn
	echo $libsToDo
	
	cd $fullPathToMAPSrun
	#cd $fullPathToMAPSrun/3B_v443_contigs
	#cd ${fullPathToMAPSrun}_includingLIB16222
	pwd
	
	# Make symlink to full path of ty mm bam file:
#	tail -n +2 $fileList | grep $libsToDo |
#	while read line; do

#		pathsToLibs=`echo $line | cut -d ' ' -f 5`

#		sampleDir=`basename $pathsToLibs`
		#echo $sampleDir
		# NB - newer seqing runs and newer libraries have longer ids
		### 29.6.0215 - need to check this again for the old and new libraries -  seem to need single quotes round the sed cmd
		### 21.7.2015 - could also test this line again for PlatesH-->N - improved the search syntax
#		samplePrefix=`echo $sampleDir | sed 's/Sample_[0-9]\{3,\}_//' | sed 's/_LDI[0-9]\{4,\}//'`
#		echo $samplePrefix
#		#echo $bamFileName
#		ls $bamFileLoctn/${sampleDir}/$bamFileName
#		if [[ ! -a ${samplePrefix}_sample.bam ]]; then ln -s $bamFileLoctn/${sampleDir}/$bamFileName ${samplePrefix}_sample.bam; fi 
#	done

:<<***string_marker_for_commenting_out_code***
	#echo "cd $PWD;
	bsub -J maps_ml$minLib -q Prod128 -n $cpu -R "rusage[mem=${mem_mb}]" -oo maps_run_-l${minLib}.log "

	set +o posix
	source /tgac/software/testing/bin/python-2.7.5
	
	#Use less than max threads for this command, gzip is called in the script for compressing temp files
	# (Run time of this step: 50 minutes, 54GB max mem) - 18.2.2016 - one of these time valeus must be wrong because the whole pipeline takes 24h or more!
#	time python /usr/users/ga002/baileyp/bin/beta-run-mpileup-TH.py \
#	-t 25 \
#	-d 8000 \
#	-r $referenceFile \
#	-o - \
#	-s /tgac/software/production/samtools/0.1.18/x86_64/bin/samtools \
#	--bamname _sample.bam \
#	| tee >(wc -l > run_mpileup.out.len) \
#	| /usr/users/ga002/baileyp/ProgramFiles/pigz-2.3.3/pigz --best -p 28 > run_mpileup.out.gz

	# (Run time of this step: 3h 20 minutes, 69MB max mem)
#	time python /usr/users/ga002/baileyp/bin/beta-mpileup-parser-TH.py \
#	-t 25 \
#	-f run_mpileup.out.gz \
#	-l run_mpileup.out.len \
#	| tee >(wc -l > parsed_run_mpileup.out.len) \
#	| /usr/users/ga002/baileyp/ProgramFiles/pigz-2.3.3/pigz --best -p 28 > parsed_run_mpileup.out.gz 

	# RUN TIME: 1.5 h, max mem 4 GB):
#	time python /usr/users/ga002/baileyp/bin/beta-maps1-v2-TH.py \
#	-t 28 \
#	-f parsed_run_mpileup.out.gz \
#	-m m \
#	-L parsed_run_mpileup.out.len \
#	-l $minLib \
#	-c $minCov \
#	-i $hetOneMinPer \
#	--maxCov $maxCov \
#	-o maps-part1_l${minLib}_v${minCov}_i${hetOneMinPer}.txt

	# RUN TIME: 40 minutes for -s 3 and 4): 
	hetMinCovs=(1 2 3 4 5 6 7 8 9 10 11 12)
	homMinCov=2		# Access this variable like so: \$homMinCov!
	for hetMinCov in \${hetMinCovs[@]}; do

		echo \"hetMinCov=\$hetMinCov	hetMinPer=$hetMinPer	homMinCov=\$homMinCov\"
	
		# bsubing these to run in ||el (RUN TIME: 3.5 hours):
		bsub -J maps2_c\${hetMinCov} -q Prod128 -n 1 -R "rusage[mem=1000]" -oo maps-part2_-l${minLib}_-c${minCov}_-d\${hetMinCov}_-p${hetMinPer}_-s\$homMinCov.log \"

		python ~baileyp/ProgramFiles/maps_fortgac/maps-part2-v2.py \
		-f maps-part1_l${minLib}_v${minCov}_i${hetOneMinPer}.txt \
		-o maps-part2_-l${minLib}_-c${minCov}_-d\${hetMinCov}_-p${hetMinPer}_-s\$homMinCov.txt \
		-m m \
		-l ${minLib} \
		-c $minCov \
		-d \$hetMinCov \
		-p $hetMinPer \
		-s \$homMinCov
		\"
	done
	homMinCov=3	# Access this variable like so: \$homMinCov!

	for hetMinCov in \${hetMinCovs[@]}; do

		echo \"hetMinCov=\$hetMinCov	hetMinPer=$hetMinPer	homMinCov=\$homMinCov\"

		bsub -J maps2_c\${hetMinCov} -q Prod128 -n 1 -R "rusage[mem=1000]" -oo maps-part2_-l${minLib}_-c${minCov}_-d\${hetMinCov}_-p${hetMinPer}_-s\$homMinCov.log \"

		python ~baileyp/ProgramFiles/maps_fortgac/maps-part2-v2.py \
		-f maps-part1_l${minLib}_v${minCov}_i${hetOneMinPer}.txt \
		-o maps-part2_-l${minLib}_-c${minCov}_-d\${hetMinCov}_-p${hetMinPer}_-s\$homMinCov.txt \
		-m m \
		-l ${minLib} \
		-c $minCov \
		-d \$hetMinCov \
		-p $hetMinPer \
		-s \$homMinCov
		\"
	done
	#9.2.2016 - Could remove mpileup files here like so:
	#if [[ -a run_mpileup.out.gz ]]; then rm run_mpileup.out.gz; fi
	#if [[ -a parsed_run_mpileup.out.gz ]]; then rm parsed_run_mpileup.out.gz; fi
	" 
	# End of main bsub
***string_marker_for_commenting_out_code***

	# NB - have to run this next step separately, after the MAPS runs are complete (outside of the above bsub!):
#:<<***string_marker_for_commenting_out_code***
	# Get contig sets from the mm bamfile for each MAPS run (RUN TIME: 1 hour per sample, up to 135 GB mem (for kronos: 10 minutes, 26 GB mem)):
	# Generate unique list of all individuals with mutations called (NB - need to set the minLib at top of script!):
	coverageLevel=het4+hom3+
	
	# Rather than put output into main directory, output files to the het/hom maps files used. 
	# Only make the MAPS run directory if it doesn't exist: 
	if [[ ! -d $coverageLevel ]]; then mkdir $coverageLevel; fi
	
	
	awk '{print $7}' maps-part2_-l${minLib}_-c10_-d4_-p10_-s3.txt | tail -n +2 | sort -u > $coverageLevel/individual_list.txt
	#awk '{print $7}' maps-part2_-l${minLib}_-c10_-d3_-p10_-s2.txt | tail -n +2 | sort -u > individual_list.txt

	cpu=1
	mem_mb=125000 	# 26000 for Kronos; 62 for Cadenza (NB - there were quite a lot of samples using memory > 80 GB up to ~ 130 GB; 135 GB was OK for all but one sample LIB9205)
	mem_gb=125
	# Collect all lines belonging to each individual from the MAPS output:
	#cat individual_list_reruns.txt | \	# For re-runs of individual libs if they fail
	cat $coverageLevel/individual_list.txt | \
	while read line; do \
		echo $line
		#cat maps-part2_-l${minLib}_-c10_-d3_-p10_-s2.txt   3B_v443_contigs/maps-part2_-l${minLib}_-c10_-d3_-p10_-s2.txt \
		cat maps-part2_-l${minLib}_-c10_-d4_-p10_-s3.txt   3B_v443_contigs/maps-part2_-l${minLib}_-c10_-d4_-p10_-s3.txt \
		| grep $line > $coverageLevel/${line}.muts

		echo "bamfile: ${line}_sample.bam"
		
		#bsub -J $line -q normal -n $cpu -R "rusage[mem=${mem_mb}]" -oo ${line}_getContigSets.bsub_l
		
		echo "cd $PWD
		source /tgac/software/production/bin/production_env
		source samtools-0.1.18
		samtools view ${line}_sample.bam \
		| ~baileyp/bin/39.getContigSets.pl -fa $referenceFile \
		-mi $coverageLevel/${line}.muts \
		-mo $coverageLevel/${line}.muts_plus_contig_sets.txt \
		-s $coverageLevel/${line}.muts_contig_sets \
		> $coverageLevel/${line}.muts_contig_sets.log 2>&1 "| qsub -N $line -q Prod -l select=1:ncpus=$cpu:mem=${mem_gb}gb -j eo
	done
#***string_marker_for_commenting_out_code***

	# The contig_sets file (minus coordinates) can be uniqified to produce a file of non-redundant contig sets:
	#cat contig_sets | sort -u > contig_sets_sort_-u
	#cat contig_sets_sort_-u | wc -l
done
***string_marker_for_commenting_out_code***




# Merge the MAPS output for each sample and get the SNPs that are unique to the ty_mm list:
:<<***string_marker_for_commenting_out_code***
cpu=1				
mem_mb=3000	# 85 GB	
mem_gb=3
runsToDo='ty_mm_MAPS_bams_to_use_LIB16208-LIB16231\|54_test_samples/ty_mm_MAPS_bams_to_use_LIB5772-LIB8427\|54_test_samples/ty_mm_MAPS_bams_to_use_LIB8428-LIB9232'
#tail -n +2 $MAPSruns | grep $runsToDo | 
# To run most runs, missing out a comments line and specific runs:
#tail -n +2 $MAPSruns | grep -v '^"#\|ty_mm_MAPS_bams_to_use_LIB16208-LIB16231\|54_test_samples/ty_mm_MAPS_bams_to_use_LIB5772-LIB8427\|54_test_samples/ty_mm_MAPS_bams_to_use_LIB8428-LIB9232' |
tail -n +2 $MAPSruns | grep $runsToDo |
while read line; do

	fullPathToMAPSrun=`echo $line | cut -d ' ' -f 1`
	minLib=`echo $line | cut -d ' ' -f 2`
	libsToDo=`echo $line | cut -d ' ' -f 3`

	bamFileLoctn=`dirname $fullPathToMAPSrun`

	echo
	echo
	echo
	ls -d $fullPathToMAPSrun 
	ls -d $bamFileLoctn
	echo "minLib=$minLib"
	echo $libsToDo
	
	cd $fullPathToMAPSrun		### 3.4.2016 - now need to cd into a subdir for a particulat het/hom value
	
	MAPS_OrigRunDirName=`basename $fullPathToMAPSrun | sed 's/ty_mm_//'`
	MAPSrun=`basename $fullPathToMAPSrun | sed 's/ty_mm_MAPS_bams_to_use_//'`
	echo $MAPS_OrigRunDirName 
	echo $MAPSrun

	# RUN TIME: 3 h 11 minutes, 85 GB mem (except 54_test_samples ty_mm_MAPS_bams_to_use_LIB8428-LIB9232 (110 GB mem):
	bsub -J $MAPSrun -q Prod128 -R "rusage[mem=${mem_mb}]" -oo ${MAPSrun}_SNP_list_mm_uniq_MAPS_lines_withPerl.log "

	cat *.muts_plus_contig_sets.txt > all_mutations_plus_contig_sets.txt

	echo 'Number of SNPs in the original run:'				# NB - Replaced ${minLib} with * because some original MAPS runs have a different minLib than for the ty_mm MAPS runs!
	tail -n +2 ../$MAPS_OrigRunDirName/*/maps-part2_-l*_-c10_-d3_-p10_-s2.txt| awk '{print \$1 \"\t\" \$2}'| grep -v '==>' | sort -u | wc -l
	tail -n +2 ../$MAPS_OrigRunDirName/*/maps-part2_-l*_-c10_-d3_-p10_-s2.txt| awk '{print \$1 \"\t\" \$2}'| grep -v '==>' | sort -u > ${MAPSrun}_SNP_list
	echo 'Number of SNPs in mm run:'
	cat all_mutations_plus_contig_sets.txt | awk '{print \$1 \"\t\" \$2}' | sort -u | wc -l 
	cat all_mutations_plus_contig_sets.txt | awk '{print \$1 \"\t\" \$2}' | sort -u > ${MAPSrun}_SNP_list_mm
	echo 'Number of SNPs common to both files:'
	comm -12 ${MAPSrun}_SNP_list ${MAPSrun}_SNP_list_mm | wc -l
	echo 'Number of SNPs unique to the mm list:'
	comm -13 ${MAPSrun}_SNP_list ${MAPSrun}_SNP_list_mm | wc -l
	comm -13 ${MAPSrun}_SNP_list ${MAPSrun}_SNP_list_mm > ${MAPSrun}_SNP_list_mm_uniq
	
	# Extract MAPS lines matching above file (RUN TIME: 3h 11 minutes, 85 GB mem; instantly using Perl):
	# NB - it is much quicker to use a Perl script than grep -f for list with > 10,000's of lines! (RUN TIME instant!)
	39.filterUniq_ty_mm_MAPS_Lines.pl -f all_mutations_plus_contig_sets.txt -l ${MAPSrun}_SNP_list_mm_uniq -o ${MAPSrun}_SNP_list_mm_uniq_MAPS_lines_withPerl
	#grep -f ${MAPSrun}_SNP_list_mm_uniq   all_mutations_plus_contig_sets.txt > ${MAPSrun}_SNP_list_mm_uniq_MAPS_lines
	echo 'Number of SNPs unique to the mm list again (after extracting from the MAPS output - just to check):'
	cat ${MAPSrun}_SNP_list_mm_uniq_MAPS_lines_withPerl | wc -l
	"
done
***string_marker_for_commenting_out_code***
### NB - in the last step above to get the ${MAPSrun}_SNP_list_mm_uniq_MAPS_lines file, there may be some lines with nested matches.
### This just means that one or two SNPs that are in the original set will get through to the final ty_mm set. 



# Now prepare large file with all MAP output from these files: *SNP_list_mm_uniq_MAPS_lines (RUN TIME: *******8********** minutes):
# Only preparing this file for the stats and VEP analysis: maps-part2_-l${minLib}_-c10_-d4_-p10_-s3.txt (i.e. the *SNP_list_mm_uniq_MAPS_lines files derived from it) 

:<<***string_marker_for_commenting_out_code***
# Count to check the number of files in each set:
ls /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/54_test_samples/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl | wc -l	# Expect 2 files
more /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/54_test_samples/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log | grep Successful | wc -l	# Expect 2 files
tail -n 12 /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/54_test_samples/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log
ls /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[ACEFG]/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl | wc -l	# Expect 14 files
more /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[ACEFG]/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log | grep Successful | wc -l	# Expect 14 files
tail -n 12 /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[ACEFG]/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log
ls /tgac/scratch/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[BD]/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl | wc -l	# Expect 6 files
more /tgac/scratch/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[BD]/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log | grep Successful | wc -l	# Expect 6 files
tail -n 12 /tgac/scratch/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[BD]/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log
ls /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[JKLM]/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl | wc -l	# Expect 12 files
more /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[JKLM]/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log | grep Successful | wc -l	# Expect 12 files
tail -n 12 /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[JKLM]/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log
ls /tgac/scratch/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[HIN]/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl | wc -l	# Expect 8 files
more /tgac/scratch/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[HIN]/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log | grep Successful | wc -l	# Expect 8 files
tail -n 12 /tgac/scratch/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[HIN]/ty_mm_MAPS_bams_to_use_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log
ls /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/PlateA-G_poor_captures_MAPS/ty_mm_*/*SNP_list_mm_uniq_MAPS_lines_withPerl | wc -l	# Expect 4 files
more /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/PlateA-G_poor_captures_MAPS/ty_mm_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log | grep Successful | wc -l	# Expect 4 files
tail -n 12 /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/PlateA-G_poor_captures_MAPS/ty_mm_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log
ls /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/PlateH2N_poor_captures_MAPS/ty_mm_*/*SNP_list_mm_uniq_MAPS_lines_withPerl | wc -l	# Expect 2 files
more /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/PlateH2N_poor_captures_MAPS/ty_mm_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log | grep Successful | wc -l	# Expect 2 files
tail -n 12 /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/PlateH2N_poor_captures_MAPS/ty_mm_*/*SNP_list_mm_uniq_MAPS_lines_withPerl.log
# Expecting 48 files
		
# Prepare final MAPS output for plates A --> N
# NB - also needed to count before removal of SNPs that are in the original set. The easiest thing to do was to use this this code 
# but use a variable for the infile and outfile:
infile=all_mutations_plus_contig_sets.txt					# *SNP_list_mm_uniq_MAPS_lines_withPerl
outfile='all_mutations_plus_contig_sets_-d3+_-s2+'		# TABCDEFGHIJKLMN_poorA-G_poorH-N_-d3+_-s2+
bsub -J mps_cat -q Prod128 -n 1 -oo ${outfile}_-d3+_-s2+_cat.log "
# Test captures (the libraries indicated are present in other parts of the data so are removed here):
cat /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/54_test_samples/ty_mm_MAPS_bams_to_use_*/$infile \
| grep -v 'LIB577[2345679]\|LIB5780' \
> ${outfile}.tab
# Plates A to G:
cat /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[ACEFG]/ty_mm_MAPS_bams_to_use_*/$infile \
>> ${outfile}.tab
cat /tgac/scratch/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[BD]/ty_mm_MAPS_bams_to_use_*/$infile \
>> ${outfile}.tab
# Plates H --> N:			
cat /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[JKLM]/ty_mm_MAPS_bams_to_use_*/$infile \
>> ${outfile}.tab
cat /tgac/scratch/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[HIN]/ty_mm_MAPS_bams_to_use_*/$infile \
>> ${outfile}.tab	
# Poor captures for plates A to G:
# 0-5k: prefiltered to keep only LIB9944\|LIB9946\|LIB9947:
cat /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/PlateA-G_poor_captures_MAPS/ty_mm_less5m_On-targReadPairs/$infile \
| grep -v 'LIB10612\|LIB10614\|LIB10615\|LIB10637\|LIB10672\|LIB10674\|LIB10675\|LIB9956\|LIB9957\|LIB9958\|LIB9959\|LIB10454\|LIB10460\|LIB10231\|LIB10235\|LIB10241\|LIB11002\|LIB11484\|LIB11486\|LIB10943\|LIB10944\|LIB10945\|LIB10946\|LIB10701' \
>> ${outfile}.tab
# 5 - 10k: removing 3 extra samples present elsewhere in the poor captures:
cat /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/PlateA-G_poor_captures_MAPS/ty_mm_5m-10m_On-targReadPairs/$infile \
| grep -v 'LIB10674\|LIB10675\|LIB10946' \
>> ${outfile}.tab
# 10 - 15k: removing 1 sample present elsewhere in the poor captures:
cat /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/PlateA-G_poor_captures_MAPS/ty_mm_10m-15m_On-targReadPairs/$infile \
| grep -v 'LIB10673' \
>> ${outfile}.tab
# 15 - 20k
cat /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/PlateA-G_poor_captures_MAPS/ty_mm_15m-20m_On-targReadPairs/$infile \
>> ${outfile}.tab
# Poor captures for plates H to N:				
# 0-5k: prefiltered to keep only LIB9944\|LIB9946\|LIB9947:
cat /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/PlateH2N_poor_captures_MAPS/ty_mm_less_5m_On-targReadPairs/$infile \
| grep -v 'LIB15640\|LIB15659\|LIB15664\|LIB15406\|LIB15407\|LIB15408\|LIB15428\|LIB15430\|LIB15478\|LIB15417\|LIB15448\|LIB16099\|LIB16148\|LIB16150\|LIB16158\|LIB15891\|LIB15895\|LIB16220\|LIB16569' \
>> ${outfile}.tab
# 5 - 10k: 
cat /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/PlateH2N_poor_captures_MAPS/ty_mm_5m-10m_On-targReadPairs/$infile \
>>${outfile}.tab
"
***string_marker_for_commenting_out_code***





# Now check that I have the expected number of libraries:
#awk '{print $7}' TABCDEFGHIJKLMN_poorA-G_poorH-N_more-d4_more-s3.txt | sort -u | wc -l
# Total number samples = 1244 - correct!





# Proceed with the SNP annotation pipeline with 35.getSNP_ty_mm_Annotation.sh script