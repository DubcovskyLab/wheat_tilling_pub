#!/usr/bin/env perl
#########################
# 39.getContigSets.pl

# Paul Bailey	9.12.2015
######################### 
use strict;
use warnings;
use Getopt::Long;


my ($fasta, $mmReads, $contigSets, $MAPS_In, $MAPS_Out, $usage);

GetOptions(
       'fa|fasta:s'    => \$fasta,                                            
       'r|reads:s' => \$mmReads,
       's|contigSets:s' => \$contigSets,
       'mi|mapsIn:s' => \$MAPS_In,
       'mo|mapsOut:s' => \$MAPS_Out,
       #'h|help:s'    => die $usage
       );

$usage="
--------------------------------------------------------------------
usage: perl script.pl -fa <fastafile> -mi <mapsInFile> -mo <mapsOutFile> -s <contigSets>  

required options:
-fa or --fasta:		fasta file for original reference
-r  or --mmreads:		filename for the multiply mapping reads
-s  or --contigSets:		list of contig sets (contig Id of the MAPS output is printed first)
-mi or --mapsInFile:		MAPS infile
-mo or --mapsOutFile:		MAPS outfile 			
-------------------------------------------------------------------- 
";





# Load fasta into a hash %sequence (required for preparing coordinates in the next loop through bam lines):
my $sequence;	# ref to hash
dateAndTimeNow('Start of loading reference into hash');
$sequence = fastaSeqsToHash($fasta, $usage);
#print "\nTesting a seq in the hash: $$sequence{IWGSC_CSS_1AL_scaff_1168724}";
dateAndTimeNow('End of loading reference into hash');




 
my %contigs;							# Data structure holding all read info with the primary hit as the top key
my $n = 0;								# Read Counter; becomes a unique id for each read
my $uniqMapReads = 0;					# Counter for the number of uniquely mapping reads
my $multiArmMapQualScore_more20 = 0;	# Number of reads with mapping quality >= 20
my $multiArmMapQualScore_less20 = 0;	# Number of reads with mapping quality < 20 - using for MAPS
my $unmappedReads = 0;					# Number of unmapped reads


# Data input from a merged bam file:
dateAndTimeNow('Start loop through bam');
while(my $line = <>)     {

	$n++;
	#if($n >= 100000)	{
	#	print "Processed up to this line:\n$line";
	#	last;
	#}
	
	my @fields = split /\t/, $line;

	# Collect all info for this read:
	my $XA_Field;		# NB the XA Flag can appear in either the 12th or 14th column of the sam output; if in 12th column, the CIGAR string reveals a higher freq of indels
	if($line =~ m/X0:i:1/ && $line =~ m/X1:i:0/)						{  $uniqMapReads++; next  }											# Read maps uniquely
	elsif($fields[11] && $fields[11] =~ m/XA:Z:/ && $fields[4] >= 20)	{  $XA_Field = $fields[11]; $multiArmMapQualScore_more20++;  }								# All multi mapped reads have mapping quality of 255 - we want these reads
	elsif($fields[13] && $fields[13] =~ m/XA:Z:/ && $fields[4] >= 20)	{  $XA_Field = $fields[13]; $multiArmMapQualScore_more20++;  }								# All multi mapped reads have mapping quality of 255 - we want these reads
	elsif($fields[11] && $fields[11] =~ m/XA:Z:/ && $fields[4] < 20)	{  $multiArmMapQualScore_less20++; next  }	# We don't want these reads but they shouldn't exist
	elsif($fields[13] && $fields[13] =~ m/XA:Z:/ && $fields[4] < 20)	{  $multiArmMapQualScore_less20++; next  }	# We don't want these reads but they shouldn't exist
	else	{  $unmappedReads++; next  }# Read must be unmapped

	
	my $primaryContig = $fields[2];
	#print "\n\nReadId: $n; \$primaryContig: $primaryContig"; 
	
	my $primaryContigHitCoord = $fields[3];
	# Use the samflag to interpret whether the top hit coordinate is on + or - strand:
	if($fields[1] & 0x10)	{  $primaryContigHitCoord = "-$primaryContigHitCoord"  } # sequence is revcom'ed and so on the - strand of the reference
	else	{  $primaryContigHitCoord = "+$primaryContigHitCoord"  } # Just makes coordinate format identical to that in the alternative hits 
	
	my $alnLen = length $fields[9];

	# Parse the alternative contigs in the XA:Z: flag.
	# Example of format:
	#XA:Z:IWGSC_CSS_5BL_scaff_10861556,-5698,101M,1;IWGSC_CSS_4BL_scaff_7033493,+3505,101M,1;IWGSC_CSS_7AL_scaff_4554281,+213,101M,1;
	$XA_Field =~ s/XA:Z://;
	my @contigsInfo = split /;/, $XA_Field;	# no need to use pop because null fields are stripped
	
	# Format the primary hit like the alternative hits then push onto alt hits to process together:
	my $primaryContigInfo = "$primaryContig,$primaryContigHitCoord,-,-";
	push @contigsInfo, $primaryContigInfo;
	#print "\n\@contigsInfo: @contigsInfo";
	  
	# Prepare a temporary array with this contigsInfo:
	my @contigIdArray;
	foreach my $contigInfo (@contigsInfo)	{
		my ($contigId, $pos, $CIGAR) = split /,/, $contigInfo;
		push @contigIdArray, $contigId;
	}
	my @sortedContigIdArray = sort { $a cmp $b } @contigIdArray;
	#print "\nContigs sorted: @sortedContigIdArray";

	# Populate main hash with all info from bam file:
	foreach my $contigInfo (@contigsInfo)	{
		
		my ($contigId, $pos, $CIGAR) = split /,/, $contigInfo;
		#my $contigLength = 0;
		my $contigLength = length($$sequence{$contigId});
		
		# Define the start + end coordinates:
		#my $start=0; my $end = 0;
        my($start, $end, $strand) = getCoordinates($pos, $alnLen, $contigLength);
		my $mapPos = "$start,$end,$strand";
		
		push @{ $contigs{ $contigId }{'reads'}{'non-uniq'}{$n}{'all-map-pos'} }, $mapPos;# Array of mapped regions for this read on current contig
		### NB - by counting size of this array, we know if there are multiple mapping to this contig
		# For the largest contig, all these regions for the current/same read all need masking (except one)
			
		if($contigId eq $primaryContig)	{# Also add mapping pos coordinates to 'primary-contig-pos' key
			push @{ $contigs{ $contigId }{'reads'}{'non-uniq'}{$n}{'primary-contig-pos'} }, $mapPos;
		}

		# Delete current contig from alt contig list (I think this deletes all entries of the same contig):
		my @copySortedContigIdArray = @sortedContigIdArray;
		my $countr = 0;
		foreach my $contig (@copySortedContigIdArray)	{
			### 17.8.2015 - NB - this is going to remove all hits for this contig - check if this matters
			### Doesn't it mean that no self contigs are stored? Yes, but you still have all the hit regions for the 
			### same contig being stored above in the "all-map-pos" key   
			if($contigId eq $contig)	{
				#my $deleted = delete $copySortedContigIdArray[$countr]; # Didin't work!
				my $deleted = splice @copySortedContigIdArray, $countr, 1;
				last;    
			}
			$countr++;
		}

		$contigs{ $contigId }{'reads'}{'non-uniq'}{$n}{'alt-contigs'} = [ @copySortedContigIdArray ];
			
		if($contigId eq $primaryContig)	{

			push @{ $contigs{ $contigId }{'reads'}{'non-uniq'}{$n}{'primary-map-pos'} }, $mapPos;
		}
			
		# Sanity check to see the areas I'm masking have high sequence similarity! 
		# NB - can remove this check once happy with the code
		# Get the substring, reverse complement if on -ve strand, then print strings to compare regions of similarity to make sure they agree
		#print "ReadId: $n;  \$contigId: $contigId:\n";
		#print "ReadSeq: $fields[9]\n";
		
#		my $similarRegion;
#		if($strand eq '-')	{
														# NB - OFFSET is not quite POS ($localStart here) - OFFSET starts at zero (i.e. no offset) so a coordinate starts at POS - 1 (i.e. == number of chars to leave off string at beginning)  
#			$similarRegion = substr $$sequence{$contigId}, ($start - 1), $alnLen;
			#print "befor_rc ", $similarRegion, "\n";
#			my $rev = reverse $similarRegion;
#			$rev =~ tr/GATC/CTAG/;
#			$similarRegion = $rev;
#		}
#		else	{
#			$similarRegion = substr $$sequence{$contigId}, ($start - 1), $alnLen;	
#		}
		#print "Strand:", $strand, " ", $similarRegion, "\n";
	}

	#print  "\nReadId: $n; printing contigs for this read only: ";
	#foreach my $contigId ( keys %contigs )	{
	#	if(exists $contigs{$contigId}{'reads'}{'non-uniq'}{$n})	{
	#		print "\nTop level contig: $contigId; alt-contigs (not self):"; 
	#		print @{$contigs{$contigId}{'reads'}{'non-uniq'}{$n}{'alt-contigs'}};
	#	}
	#}
}
#######close MM_READS;
dateAndTimeNow('End of loop through bam');
print "\n# of reads in bam file\t$n";
print "\n# of uniquely mapped reads\t$uniqMapReads";
print "\n# of multi mapped reads (map qual >= 20)\t$multiArmMapQualScore_more20";
print "\n# of multi mapped reads (map qual < 20 - used for masking)\t$multiArmMapQualScore_less20";
print "\n# of unmapped reads\t$unmappedReads";
print "\n\n\n";



open MAPS_IN, $MAPS_In or die "cannot open MAPS in file\n$usage";
open MAPS_OUT, ">$MAPS_Out" or die "cannot open MAPS out file\n$usage";
open CONTIG_SETS, ">$contigSets" or die "cannot open MAPS out file\n$usage";



my $SNPCountr = 0;
my $contigsCountr = 0;
foreach my $line (<MAPS_IN>)	{
	
	chomp $line;

	$SNPCountr++;

	my($contigId, $SNP_Pos) = split /\t/, $line;

	if(exists $contigs{$contigId} )	{# Only deal with contigs that are in the bam file(s)

		$contigsCountr++;
		print "\n********** $contigId **********";
		
		my %altContigHash;
	  	READLOOP: foreach my $n (keys %{ $contigs{$contigId}{'reads'}{'non-uniq'} })	{ # Foreach read that maps to the current contig
	  	
	  		foreach my $mapPos ( @{ $contigs{ $contigId }{'reads'}{'non-uniq'}{$n}{'all-map-pos'} } )	{# Foreach position read maps to in current contig
	  		
	  			# Get the coordinates:
	  			my($start, $end, $strand) = split /,/, $mapPos;
	  			 
	  			if($SNP_Pos >= $start && $SNP_Pos <= $end)	{
		  			
		  			#print "\nRead No.\$n: $n; mapPos:  $mapPos";
	  			
		  			# Get the alt contigs for this read and store in hash.
		  			foreach my $altHit (@{$contigs{$contigId}{'reads'}{'non-uniq'}{$n}{'alt-contigs'}})	{
			  			
		  				$altContigHash{$altHit} = '';
  					}
  					#last READLOOP; # Commented out to collect alt contigs from all reads
				}
			}
		}
		# Print the line and the alt contigs for each read:
		if (%altContigHash)	{# Only print if hash has any elements but i think it always will have
			#print "\n$contigId - Read Id: $n", $line, "\t"; 
			print "\n", $line, "\t";
			#print MAPS_OUT $line, "\t";
			my $contigSet = $line . "\t"; 
			print CONTIG_SETS $contigId; 
			foreach my $key (sort keys %altContigHash)	{
				print $key, ",";
				#print MAPS_OUT $key, ",";	
				$contigSet = $contigSet . $key . ",";						
				print CONTIG_SETS ",", $key; # Still need to uniqify the contig_sets file outside this script
			}
			$contigSet =~ s/,$//; # Remove the last comma
			print MAPS_OUT "$contigSet\n";
			#print MAPS_OUT "\n";
			print CONTIG_SETS "\n";
		}
	}
}
print "\nTotal number of SNPs in MAPS infile: $SNPCountr";
print "\nNumber of SNPs on contigs with multiple read mappings: $contigsCountr\n";





#####################
sub fastaSeqsToHash	{
#####################

my ($fastaFile, $usage) = @_;

	#open fasta, store headers and sequence
	my %sequence;
	my @query_info;
	my $id;

	open (FASTAFILE, "<", $fastaFile) or die "cannot open fasta file\n$usage";

	# load fasta into a hash %sequence
	my $contigCountr = 0;
	while ( my $line = <FASTAFILE> ){
    
		if($line =~ m /^>/)	{
			$contigCountr++;
		}
	
    	chomp $line; 

    	if ($line=~ m/\>/)	{

			@query_info = split (/\s/, $line);

			$id = shift(@query_info);

			$id =~ s/\>//;
    	}

    	elsif($sequence{$id}){

			#$sequence{$id} = $sequence{$id} . $line;
			$sequence{$id} .= $line;	# This concatenate is much much much much quicker, the perfect solution for the chr 3B pseudomolecule!
    	}
    	else	{
			$sequence{$id} = $line;
    	}
    
    	#if($contigCountr == 1000)	{
		#	last;
		#}
	}
	close FASTAFILE;
	return \%sequence;
}# End of fastaSeqsToHash subR





#####################
sub getCoordinates	{
### NBNB - getting the coordinates can be greatly simplified now I understand that the the '-' indicates strand and 'pos' can be used
#####################

	my($pos, $alnLen, $contigLength)  = @_;

	my $start = 0;
	my $end = 0;
	my $strand;
	if ($pos < 0)	{ 

		# This works; the -ve strand coordinate is in the same location as on the +ve strand.
		# The "-" char just informs that read matches on the -ve strand of ref, pointing upstream
		$end = 0 - $pos + $alnLen; # NB minus minus $pos is +ve
		# When adding $aln_length to $pos, if there are deletions in the reference, $end could end up
		# being longer than the contig length. If it is, adjusting $end to be the length of the contig
		if($end > $contigLength)	{
				$end = $contigLength;
		}
		$start = 0 - $pos; 
		$strand = '-';
	}
	else	{

		$start = $pos;
		$end = $pos + $alnLen;
		# See just above for explanation of this conditional:
		if($end > $contigLength)	{
			$end = $contigLength;
		}
		$strand =  '+';
	}
	return $start, $end, $strand;
}# End of getCoordinates subR





#################
sub dateAndTimeNow	{
#################

	my ($startEnd) = @_;
	my($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime;
	my $month=$mon + 1;				# $mon is zero-based, as is $wday
	my $realYear = (1900 + $year);	# 4 digit year is 1900 + $year
	print "$startEnd time: $mday.$month.$realYear $hour:$min:$sec\n";
}
