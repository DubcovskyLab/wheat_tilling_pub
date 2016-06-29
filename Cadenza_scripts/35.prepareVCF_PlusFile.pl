#!/usr/bin/perl
use strict;
use warnings;

my $table = shift;
my $IWGSC1_Coords= shift;

my @homoeologs123Sets = (1, 2, 3);
my %homoeologsOn123Hash;
foreach my $set	(@homoeologs123Sets)	{

	open CONTIGS_ALL_HOMOEOLOGS_FOR_SNPS, "/tgac/workarea/group-cg/baileyp/WheatLoLa/27.ExonCapture_exon_fishing/genes_with_${set}_homoeologs_for_SNPs" or die "genes_with_${set}_homoeologs_for_SNPs: $!"; 

	while(my $lineForHash = <CONTIGS_ALL_HOMOEOLOGS_FOR_SNPS>)	{
	
		my @fields = split /\t/, $lineForHash;
		# There are more than one gene on some contigs so will need to store them all:
		my $contigID_Start_Stop = join "\t", $fields[0], $fields[1], $fields[2];

		if(!exists $homoeologsOn123Hash{$set}{$fields[0]}{$contigID_Start_Stop})	{
			$homoeologsOn123Hash{$set}{$fields[0]}{$contigID_Start_Stop} = $contigID_Start_Stop;
			#print "\n\$homoeologsOn123Hash{$fields[0]}{$contigID_Start_Stop}: $homoeologsOn123Hash{$fields[0]}{$contigID_Start_Stop}";
		}
	}
}

# Getting Cadenza line Ids from the exon capture table:
#open TABLE, "/tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_300_Samples/35.ExonCapture_PreparingTable.txt" or die "35.ExonCapture_PreparingTable.txt : $!";
# Looks like I have to merge all the tables I have for this task!  
open TABLE, $table or die "$table : $!";

my %tableHash;
while(my $lineForHash = <TABLE>)	{
	
	my @fields = split /\t/, $lineForHash;
	chomp @fields;
	# Reformat the library name:
	$fields[1] =~ s/^...._LIB/LIB/;
	$fields[1] =~ s/_LDI\d{4,}//;	# Plates H to N have a longer LDI number - 5 digits!
	$lineForHash = join "\t", @fields;
	#print "\n", $fields[1];
	#if(exists $tableHash{$fields[1]}) 	{
		#print "\nLibrary found twice: $fields[1]   $lineForHash\n";
		### 26.1.2015 -  At the moment 4 libraires appear twice in the list:
		### LIB9980, LIB9981, LIB9982, LIB9983 - these libs were redone
		### I checked that they had the same Cadenza IDs numbers.
		### So we can ignore thee and not exit
		#exit;	
	#}
	$tableHash{$fields[1]} = $lineForHash;
}


# Also getting the new chromosome position in IWGSC_v2 pseudochromosomes:
# Line format example:
# scaffold:IWGSC2:IWGSC_CSS_1AL_scaff_1000404:4355:4355:1,chromosome:IWGSC2:1A:231871519:231871519:1	<-- need to pass this line to 35.prepareVCF_PlusFile.pl 

open IWGSC1_COORDS, $IWGSC1_Coords or die "$IWGSC1_Coords : $!";

my %iwgscCoordsHash;
while(my $lineForHash = <IWGSC1_COORDS>)	{
	
	my @fields = split /:/, $lineForHash;
	chomp @fields;
	
	my $iwgsc_v1ContigId = $fields[2];
	my $iwgsc_v1ContigCoord = $fields[3];
	my $iwgsc_v2ContigCoord = $fields[8];
	my $iwgsc_v2ChrArm = $fields[7];
	#print "\n$iwgsc_v1ContigId	$iwgsc_v1ContigCoord	$iwgsc_v2ContigCoord $iwgsc_v2ChrArm";	
	
	# Populate hash:
	$iwgscCoordsHash{$iwgsc_v1ContigId}{$iwgsc_v1ContigCoord} = $iwgsc_v2ContigCoord;
	#print "\n", $iwgscCoordsHash{$iwgsc_v1ContigId}{$iwgsc_v1ContigCoord};  
}
#exit;


# Main code:
my $line;
while($line = <>)	{# MAPS output

	my @fields = split /\t/, $line;
	my $library = $fields[6];
	### 26.1.2015 - Having an issue with 10 lines having an empty $library.
	### No line then gets printed - No idea wheats going on
	### Will print a line with a message for the moment.
	### 9.2.2015 - the file being looped through is a cat of all the MAPS output files.
	### There could be extra empty lines at the end of each file!   
	if(!$library)	{
		print "No information available for this SNP!\t\t\t\t\t\t\t\t\t\t\t\t\n";
		next;
	}
	my $cadenzaLineId;
	my $setMembership = ' ';
	
	# Get position on the IWGSC_v2 wheat pseudomolecules:
	my $IWGSC_v1_Id = $fields[0];
	my $IWGSC_v1_Coord = $fields[1]; 
	my @fields1 = split '_', $fields[0];
	my $chrArm = $fields1[2];
	my $IWGSC_v2_Coord = '.'; 
	
	# Fetch the position for the chromosome pseudomolecule, if there is one:
	if($IWGSC_v1_Id eq '3B' || $IWGSC_v1_Id =~ m/v443_/)	{
		$IWGSC_v2_Coord = $IWGSC_v1_Coord;
		$IWGSC_v1_Coord = '.';
		$chrArm = '3B'
	} 
	elsif( exists $iwgscCoordsHash{$IWGSC_v1_Id}{$IWGSC_v1_Coord} ) 	{
		$IWGSC_v2_Coord = $iwgscCoordsHash{$IWGSC_v1_Id}{$IWGSC_v1_Coord};
	}  
	
#	if( $fields[2] =~ m/[\.\+\*\d]/ || $fields[5] =~ /[\.\+\*\d]/ )	{	
#		next;
#	}	
	if(exists $tableHash{$library} )	{
			
		my @hashFields = split /\t/, $tableHash{$library};
			
		if($hashFields[5])	{
			$cadenzaLineId = $hashFields[5];
		}
		else	{
			$cadenzaLineId = '.';
		}
	}
	else	{
		$cadenzaLineId = '.';
	}
	#Testing whether SNP is on a gene with 1,2 or 3 homoeologs:
	foreach my $set	(@homoeologs123Sets)	{
	
		if(exists $homoeologsOn123Hash{$set}{$fields[0]})	{
		#print "\nHowdee";
		
			foreach my $contigsSites ( keys %{$homoeologsOn123Hash{$set}{$fields[0]} } )	{
				my @localFields = split /\t/, $homoeologsOn123Hash{$set}{$fields[0]}{$contigsSites};
				#print "\nHowdee: $homoeologsOn123Hash{$set}{$fields[0]}{$contigsSites}";
				if($fields[1] >= $localFields[1] && $fields[1] <= $localFields[2] )	{
					# NB - the coordinates from the annotation file are are just start/stop of the mRNA

					$setMembership = $setMembership . $set . ' ';		# enables to look for contflicting results between these homoeolog sets
				}
			}
		}	
	}
	# Also determining confidence of the mutation based on hetMinCov:
	my $confidence = 0;
#	if($fields[9] >= 6)		{ $confidence = 'high'}
#	elsif($fields[9] == 5 && $fields[7] eq 'het')	{ $confidence = 'medium'}
#	elsif($fields[9] == 4 && $fields[7] eq 'het')	{ $confidence = 'low'}
#	elsif($fields[9] == 3 && $fields[7] eq 'het')	{ $confidence = 'low'}
#	elsif($fields[9] == 5 && $fields[7] eq 'hom')	{ $confidence = 'high'}
#	elsif($fields[9] == 4 && $fields[7] eq 'hom')	{ $confidence = 'high'}
#	elsif($fields[9] == 3 && $fields[7] eq 'hom')	{ $confidence = 'high'}

	if($fields[9] >= 7)		{ $confidence = '99%'}
	elsif($fields[9] == 6 && $fields[7] eq 'het')	{ $confidence = '98%'}	
	elsif($fields[9] == 5 && $fields[7] eq 'het')	{ $confidence = '97%'}
	elsif($fields[9] == 4 && $fields[7] eq 'het')	{ $confidence = '84%'}
	elsif($fields[9] == 3 && $fields[7] eq 'het')	{ $confidence = '64%'}
	
	elsif($fields[9] == 6 && $fields[7] eq 'hom')	{ $confidence = '99%'}	
	elsif($fields[9] == 5 && $fields[7] eq 'hom')	{ $confidence = '99%'}
	elsif($fields[9] == 4 && $fields[7] eq 'hom')	{ $confidence = '98%'}
	elsif($fields[9] == 3 && $fields[7] eq 'hom')	{ $confidence = '97%'}
	elsif($fields[9] == 2 && $fields[7] eq 'hom')	{ $confidence = '73%'}


	# Formatting the vcf_plus line and inserting the cadenza line Id and number of homoeologs present for each gene triad:
	# Also adding in column for chrArm and position in the chromosome pseudomolecule:
	print "$fields[0]\t$chrArm\t$library\t$cadenzaLineId\t$IWGSC_v1_Coord\t$IWGSC_v2_Coord\t$fields[2]\t$fields[4]\t$fields[5]\t$fields[7]\t$fields[8]\t$fields[9]\t$confidence\t";
	
	# This $setMembership variable is giving conflicting results - in which case, dont print anything: 
	if( ($setMembership =~ m/1/ && $setMembership	=~ m/2/) || ($setMembership =~ m/1/ && $setMembership =~ m/3/) || ($setMembership =~ m/2/ && $setMembership	=~ m/3/) )	{
		$setMembership = ' ';
	}
	if($setMembership eq ' ')	{ # i.e. if there is no homoeolog number defined 
		print "$setMembership\t\n";
	}
	else	{
		my @fields = split ' ', $setMembership; # NB - split splits on white space after skipping leading whitespace, therefore I want to print $field[0] 
		print "$fields[0]\t\n";
	}
}