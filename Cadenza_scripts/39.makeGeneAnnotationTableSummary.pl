#!/usr/bin/perl
##################################################
my $script='39.makeGeneAnnotationTableSummary.pl';

# Modified from 41.getFastaLengths.pl

# Paul Bailey	23.2.2016
##################################################

use strict;
use warnings;
use Getopt::Long;

my ($infile, $table, $table1, $table2, $table3, $table4, $table5, $table6, $table7, $table8, $table9, $table10, $table11, $table12, $table13, $table14, $table15, $table16, $table17, $table18, $table19, $table20, $table21, $outfile);


my $usage="

program description: takes in a list of wheat trancript Ids (from a table) and adds 
                     information associated with each gene - in this case - number 
                     of mutant variant types per gene plus other info, including from 
                     BioMart. 

usage: perl $script -f <infile> -t <table> - outfile.tab

required options:
-f or --infile:		list of wheat trancript Ids - format: length	geneId	gc	geneId
-t or --table:		gene_list_numbr_muts - format: count	transcriptId
-t1 or --table1:	gene_list_numbr_stop_gained - format: count	transcriptId
-t2 or --table2:	gene_list_numbr_splice_donor_variant - format: count	transcriptId
-t3 or --table3:	gene_list_numbr_splice_acceptor_variant - format: count	transcriptId
-t4 or --table4:	gene_list_numbr_missense_variant - format: count	transcriptId
-t5 or --table5:	gene_list_numbr_synonymous_variant - format: count	transcriptId
-t6 or --table6:	gene_list_numbr_downstream_gene_variant - format: count	transcriptId
-t7 or --table7:	gene_list_numbr_upstream_gene_variant - format: count	transcriptId
-t8 or --table8:	gene_list_numbr_5_prime_UTR_variant - format: count	transcriptId
-t9 or --table9:	gene_list_numbr_3_prime_UTR_varian - format: count	transcriptId
-t10 or --table10:	gene_list_numbr_initiator_codon_variant - format: count	transcriptId
-t11 or --table11:	BioMart data - tab format: e.g. GeneStableID	TranscriptStableI	InterProDescription	PfamID	Chromosome/ScaffoldName	GeneStart(bp)	GeneEnd(bp)	%GC_Content
-t12 or --table12:	number of exons per gene (from Ricardo) - tab format: gene	avg_cov	bases_in_exons	exon_count
-t13 or --table13:	average read coverage for each gene - tab  format: coverage	transcriptId (From MAPS output)
-t14 or --table14:	percent EMS per gene - tab format: percent_EMS	transcriptId
-t15 or --table15: 	iwgsc_v1ContigId_gene_list - format: count	contigId	transcriptId
-t16 or --table16:	list of genes with deletions (from Ricardo): tab format: transcriptId	LIBRARY	CadenzaId	contigId
-t17 or --table17:	average read coverage for each gene - tab  format: coverage	transcriptId (From 28 bam files)
-t18 or --table18:	gene_list_numbr_frameshift - format: count	transcriptId
-t19 or --table19:	gene_list_numbr_muts_NO_up_downstream - format: count	transcriptId
-t20 or --table20:	gene_list_numbr_sift_less0.05 - format: count	transcriptId
-t21 or --table21:	gene_list_numbr_sift_less0.01 - format: count	transcriptId
-o or --outfile:	list of wheat trancript Ids with variant counts for each gene - format: length	geneId	gc	geneId	numbr_muts	numbr_stop_gained etc
";
# NB - 25.4.2016 - for next time, it's best to prepare these files with columns separated by tabs - would help with parsing
# NB - 25.4.2016 - I think the script could be made generic by dealing loading x number of table columns - all must start off in tab-delimited format


GetOptions(
       'f|infile:s'    => \$infile,
       't|table:s'    => \$table,
       't1|table1:s'    => \$table1,
       't2|table2:s'    => \$table2,
       't3|table3:s'    => \$table3, 
       't4|table4:s'    => \$table4,
       't5|table5:s'    => \$table5,
       't6|table6:s'    => \$table6,
       't7|table7:s'    => \$table7,
       't8|table8:s'    => \$table8,
       't9|table9:s'    => \$table9,
       't10|table10:s'    => \$table10,
       't11|table11:s'    => \$table11, 
       't12|table12:s'    => \$table12,
       't13|table13:s'    => \$table13,
       't14|table14:s'    => \$table14,
	   't15|table15:s'    => \$table15,
	   't16|table16:s'    => \$table16,
	   't17|table17:s'    => \$table17,
	   't18|table18:s'    => \$table18,
	   't19|table19:s'    => \$table19,
	   't20|table20:s'    => \$table20,
	   't21|table21:s'    => \$table21,
        'o|outfile:s'    => \$outfile,                             
       #'h|help:s'    => die $usage
       );

# Create Hash(es):
my %tableHash = createHash($table);
my %tableHash1 = createHash($table1);
my %tableHash2 = createHash($table2);
my %tableHash3 = createHash($table3);
my %tableHash4 = createHash($table4);
my %tableHash5 = createHash($table5);
my %tableHash6 = createHash($table6);
my %tableHash7 = createHash($table7);
my %tableHash8 = createHash($table8);
my %tableHash9 = createHash($table9);
my %tableHash10 = createHash($table10);
my %tableHash18 = createHash($table18);
my %tableHash19 = createHash($table19);
my %tableHash20 = createHash($table20);
my %tableHash21 = createHash($table21);

my %tableHash13 = createHashFromTabFormat($table13);
my %tableHash14 = createHashFromTabFormat($table14);
my %tableHash17 = createHashFromTabFormat($table17);



# Need to create BioMart hash separately:
open TABLE, $table11 or die "$table11: $!\n\n$usage";
my %tableHash11;
while(my $line = <TABLE>)	{

	chomp $line;
	#print $line, "\n";
	my($geneId, $transcriptId, $InterProDesc, $pFAM, $IWGSC_ContigId, $start, $end, $GC_Content) = split /\t/, $line;
	#print $id, "\t", $count, "\n";
	if(exists $tableHash11{$transcriptId})	{
		#print "$transcriptId - this contig already exits ...", "\n";
		$tableHash11{$transcriptId}{pFAM}{$pFAM} = $pFAM;
		$tableHash11{$transcriptId}{interProDesc} = join '; ', $tableHash11{$transcriptId}{interProDesc}, $InterProDesc;
	}
	else	{
		$tableHash11{$transcriptId}{interProDesc} = $InterProDesc;  
		$tableHash11{$transcriptId}{IWGSC_ContigId} = $IWGSC_ContigId;
		$tableHash11{$transcriptId}{pFAM}{$pFAM} = $pFAM;
		$tableHash11{$transcriptId}{start} = $start;
		$tableHash11{$transcriptId}{end} = $end;
		$tableHash11{$transcriptId}{GC_Content} = $GC_Content;
		#print $id, "\t", $tableHash{$id}, "\n"
	}
	#foreach my $key ( keys %{$tableHash11{$transcriptId}{pFAM} } )	{ print "$key \n" }
	#print "pFAM Ids:", keys %{$tableHash11{$transcriptId}{pFAM} }, "\n"
}


# Need to create gene coverage exon number hash separately:
open TABLE, $table12 or die "$table12: $!\n\n$usage";	
my %tableHash12;
while(my $line = <TABLE>)	{

	chomp $line;
	#print $line, "\n";
	my($transcriptId, $avg_cov, $bases_in_exons, $numbrExons) = split /\t/, $line;
	#print $id, "\t", $count, "\n";
	if(exists $tableHash12{$transcriptId})	{ print "$transcriptId - this id already exits ... exiting\n"; exit }
	else	{
		$tableHash12{$transcriptId}{avg_cov} = $avg_cov;
		$tableHash12{$transcriptId}{number_exons} = $numbrExons;
		$tableHash12{$transcriptId}{bases_in_exons} = $bases_in_exons;
		#print $id, "\t", $tableHash{$id}, "\n"
	}
	#foreach my $key ( keys %{$tableHash11{$transcriptId}{pFAM} } )	{ print "$key \n" }
	#print "pFAM Ids:", keys %{$tableHash11{$transcriptId}{pFAM} }, "\n"
}


# Need to create iwgsc_v1 ContigId geneId list hash separately:
open TABLE, $table15 or die "$table15: $!\n\n$usage";
my %tableHash15;
while(my $line = <TABLE>)	{

	chomp $line;
	#print $line, "\n";
	my($contigId, $transcriptId) = split /\t/, $line;
	#print $id, "\t", $count, "\n";
	if(exists $tableHash15{$transcriptId})	{ print "$transcriptId - this id already exits ... exiting\n"; exit }
	else	{
		$tableHash15{$transcriptId} = $contigId;
	}
}


# Create hash of deletions from Ricardo separately:
open TABLE, $table16 or die "$table16: $!\n\n$usage";	
my %tableHash16;
while(my $line = <TABLE>)	{

	chomp $line;
	#print $line, "\n";
	my($transcriptId, $LIBRARY, $CadenzaId, $contigId) = split /\t/, $line;
	#print $id, "\t", $count, "\n";
	if(exists $tableHash16{$transcriptId})	{
		#print "$transcriptId - this contig already exits ...", "\n";
		$tableHash16{$transcriptId} = join '; ', $tableHash16{$transcriptId}, $CadenzaId;
	}
	else	{
		$tableHash16{$transcriptId} = $CadenzaId;
	}
}



# Main code:
open INFILE, $infile or die "$infile: $!\n\n$usage";
open OUTFILE, ">$outfile" or die "$outfile: $!\n\n$usage";

# Header line:
print OUTFILE "transcriptId\tlength_of_CDS\t%gc\t#bases_in_exons\t#exons\tavg_cov(FromMAPS_Output)\tavg_cov_CDS(From28BamFiles)\tIWGSC_v1ContigId\t#muts\t%_EMS\t#muts_MINUS_up_down_stream_variants\t#stop_gained\t#splice_donor_variant\t#splice_acceptor_variant\t#missense_variant\t#missense_variant_sift<0.05\t#missense_variant_sift<0.01\t#synonymous_variant\t#downstream_gene_variant\t#upstream_gene_variant\t#5_prime_UTR_variant\t#3_prime_UTR_variant\t#initiator_codon_variant\t#OthrVariants\tLargeDeletions\tSmallDeletions(main_dataset_het5+_hom3+_only)\tInterProDesc\tpFAM_Ids\tIWGSC_gContigId\tgstart\tgend\t%gc(BioMart)\n";


while(my $line = <INFILE>)	{
	
	chomp $line;
	#print $line, "\n";
	my($length, $id, $gc, $id2) = split /\t/, $line;
	#print $id, "\t", $length, "\t", $gc, "\n";

	if(exists $tableHash{$id})	{ print OUTFILE $id, "\t", $length, "\t", $gc}
	else	{ print OUTFILE $id, "\t", $length, "\t", $gc}		# NB - 25.4.2016 - What's the purpose of printing exactly the same?
	
	if(exists $tableHash12{$id})	{ print OUTFILE "\t", $tableHash12{$id}{bases_in_exons}, "\t", $tableHash12{$id}{number_exons} }
	else	{ print OUTFILE "\t", "-", "\t", "-" }

	printDataFromHashes($id, $length, $gc, \%tableHash13); # coverage (from MAPS)
	
	printDataFromHashes($id, $length, $gc, \%tableHash17); # coverage (from bam files)	
	
	if(exists $tableHash15{$id})	{ print OUTFILE "\t", $tableHash15{$id} } # IWGSC_v1 contig ids
	else	{ print OUTFILE "\t", "-"}
	
	if(exists $tableHash{$id})	{ print OUTFILE "\t", $tableHash{$id}} #  # muts
	else	{ print OUTFILE "\t", "0"}
	
	printDataFromHashes($id, $length, $gc, \%tableHash14); # % EMS
	
	printDataFromHashes($id, $length, $gc, \%tableHash19); # number mutations minus up and downstream effects
	
	printDataFromHashes($id, $length, $gc, \%tableHash1);	# NB 25.4.2016 - why did I send $length and $gc to the printing subR?
	printDataFromHashes($id, $length, $gc, \%tableHash2);
	printDataFromHashes($id, $length, $gc, \%tableHash3);
	printDataFromHashes($id, $length, $gc, \%tableHash4);	# missense variants
	
	printDataFromHashes($id, $length, $gc, \%tableHash20);	# sift score <0.05
	printDataFromHashes($id, $length, $gc, \%tableHash21);	# sift score <0.01
	
	
	printDataFromHashes($id, $length, $gc, \%tableHash5);
	printDataFromHashes($id, $length, $gc, \%tableHash6);
	printDataFromHashes($id, $length, $gc, \%tableHash7);
	printDataFromHashes($id, $length, $gc, \%tableHash8);
	printDataFromHashes($id, $length, $gc, \%tableHash9);
	printDataFromHashes($id, $length, $gc, \%tableHash10);

	if(exists $tableHash16{$id})	{ print OUTFILE "\t", $tableHash16{$id} }
	else	{ print OUTFILE "\t", "-"}
	
	printDataFromHashes($id, $length, $gc, \%tableHash18); #small deletions


	if(exists $tableHash11{$id})	{
		print OUTFILE "\t", $tableHash11{$id}{interProDesc}, "\t", join " ", keys %{$tableHash11{$id}{pFAM} }, "\t", $tableHash11{$id}{IWGSC_ContigId}, "\t", $tableHash11{$id}{start}, "\t", $tableHash11{$id}{end}, "\t", $tableHash11{$id}{GC_Content};
	}
	else	{ print OUTFILE "\t", "-", "\t", "-", "\t", "-", "\t", "-", "\t", "-", "\t", "-" }
	
	print OUTFILE "\n"
}


# Subroutines:
#################
sub createHash	{
#################

	my($table) = @_;
	open TABLE, $table or die "$table: $!\n\n$usage";	
	my %tableHash;
	while(my $line = <TABLE>)	{

		chomp $line;
		#print $line, "\n";
		my($blank, $count, $id) = split /\s+/, $line;
		#print $id, "\t", $count, "\n";
		if(exists $tableHash{$id})	{ print "$table\t$id - this id already exits ... exiting\n"; exit }
		else	{
			$tableHash{$id} = $count;
			#print $id, "\t", $tableHash{$id}, "\n"
		}
	}
	return %tableHash
}


#############################
sub createHashFromTabFormat	{
#############################

	# Loading gene coverage column separately - not quite the same format as for e.g. gene_list_numbr_muts
	my($table) = @_;
	open TABLE, $table or die "$table: $!\n\n$usage";
	my %tableHash;
	while(my $line = <TABLE>)	{

		chomp $line;
		#print $line, "\n";
		my($value, $id) = split /\t/, $line;
		#print $id, "\t", $value, "\n";
		if(exists $tableHash{$id})	{ print "$table\t$id - this id already exits ... exiting\n"; exit }
		else	{
			$tableHash{$id} = $value;
		}
	}
	return %tableHash
}



#########################
sub printDataFromHashes	{
#########################
	
	my($id, $length, $gc, $tableHash) = @_;
	
	if(exists $$tableHash{$id})	{ print OUTFILE "\t", $$tableHash{$id} }
	else	{ print OUTFILE "\t", "0" }	
}