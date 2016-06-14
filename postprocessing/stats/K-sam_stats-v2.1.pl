#!/usr/bin/perl -w
# Author        : Ksenia V Krasileva
 
# new in v2.1
# option for controlling minimum depth coverage to report
# new in v1.1
# -s/--sam input option changed to -b/--bam
# takes in -o/--output option and prints results in the output file instead of sdout
 
use strict;
use warnings;
 
#use Bio::Seq;
#use Bio::SeqIO;
#use Bio::SearchIO;
use Getopt::Long;
 
my $usage="
Usage: perl K-sam-stats-v2.1.pl 
-b <mapping.sorted.bam> 
-f <reference.fasta> 
-c <minimum reference coverage to report [default=1]> 
-o <output.txt>";
 
my ($filein, $fileout, $fasta, $mincov);
 
GetOptions(
       'b|bam:s'        => \$filein,
       'o|output:s'        => \$fileout,
       'f|fasta:s'        => \$fasta,
       'c|mincov:s'        => \$mincov,
);  
 
die $usage unless(defined($filein) and defined($fasta));
 
$mincov=1 unless defined($mincov);
 
if ( defined($fileout) ){
    open (FILEOUT1, ">", $fileout);
    }
else{
    open (FILEOUT1, ">", $filein."-summary.txt");
    }
     
#load fasta reference and store it in a hash
 
open (FASTAFILE, "<", $fasta);
 
my %sequence;
my $header;
 
while ( my $line = <FASTAFILE> ){
     
chomp $line; 
 
    if ($line=~ m/\>/)   { ($header) = split (/\s/, $line); $header =~ s/>//;}
 
    elsif($sequence{$header}) { $sequence{$header} = $sequence{$header} . $line; }
 
    else{ $sequence{$header} = $line; }
 
}
 
close FASTAFILE;
 
#calculate total length of the reference
 
my $ref_total = 0;
 
foreach my $key (keys %sequence){ $ref_total += length($sequence{$key}); }
 
#calculating basic statistics on bam file with `samtools view -c`
 
my $total_aln = `samtools view -c $filein`;
my $mapped_aln = `samtools view -c -F4 $filein`;
my $secondary_aln = `samtools view -c -f 256 $filein`;
my $unmapped = `samtools view -c -f4 $filein`;
my $mapped_pairs = `samtools view -c -F256 -f2 $filein`; #considering primary alignments only
 
chomp($total_aln, $mapped_aln, $secondary_aln, $unmapped, $mapped_pairs);
 
my $total_reads = $unmapped + $mapped_aln - $secondary_aln;
my $mapped_reads = $mapped_aln - $secondary_aln;
 
#mapping quality
 
my $mq3 = `samtools view -c -q 3 $filein`;
my $mq10 = `samtools view -c -q 10 $filein`;
my $mq20 = `samtools view -c -q 20 $filein`;
my $mq30 = `samtools view -c -q 30 $filein`;
 
chomp ($mq3, $mq10, $mq20, $mq30);
 
#calculating number of clipped reads
 
open (PS, "samtools view -F4 $filein | ");
 
my $mapped_sclipped = 0;
 
while (my $line=<PS>){
 
    chomp $line;
         
    if ($line =~ m/^@/){    next;    }
     
    else{
     
    my ($qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual) = split (/\t/, $line);
 
        if ($cigar =~ m/S/) {$mapped_sclipped++;}
         
        }
    }
 
close PS;
 
#calculating the coverage
 
open (GCOV, "bedtools genomecov -d -ibam $filein | ");
 
my $ref_length_1 = 0;
 
while (my $line=<GCOV>){
 
    chomp $line;
     
    my ($ref, $pos, $cov) = split (/\t/, $line);
 
    if ($cov >= $mincov){ $ref_length_1++};
             
}
 
close GCOV;
     
print FILEOUT1"
 
       File analyzed: $filein
 
       READ STATISTICS\n
       Total number of reads: ", $total_reads,"
       Number of unmapped reads: $unmapped
       Number of mapped reads: ", $mapped_reads, " ( ", $mapped_reads/$total_reads, " )      
       Number of reads mapped in pairs: ", $mapped_pairs, " ( ", $mapped_pairs/$total_reads, " )
        
       ALIGNMENT STATISTICS\n
       Number of alignments: $mapped_aln
       Number of secondary alignments: $secondary_aln
       Number of soft-clipped alignments: $mapped_sclipped
 
       Mapping quality =>3: $mq3
       Mapping quality =>10: $mq10
       Mapping quality =>20: $mq20
       Mapping quality =>30: $mq30
        
       REFERENCE STATISTICS\n
       Total length of the reference: $ref_total
       Length of reference covered >$mincov read: ", $ref_length_1, " ( ", $ref_length_1/$ref_total, " )  
 
"; 
 
close FILEOUT1
