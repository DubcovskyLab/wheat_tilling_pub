#!/bin/bash
## By: Hans Vasquez-Gross

IN_FILES=*.tsv

for FILE in $IN_FILES
do

    prefix=$(basename $FILE .tsv)
    awkcmd="awk -F '\t' '{print \$0\"\t\"\$7;}' $FILE | sed 's/_bin_[0-9]*-[0-9]*\t/\t/' > $prefix.fixed.tsv" 
    echo $awkcmd
    eval $awkcmd 

    sedcmd="sed -i.bak 's/Lib/Gene/2' $prefix.fixed.tsv"
    echo $sedcmd
    eval $sedcmd

    vcfcmd="python generate_VCF_file_from_MAPs.py $prefix.fixed.tsv $prefix.unsorted.vcf $prefix.unsorted.indel.vcf"
    echo $vcfcmd
    eval $vcfcmd

    vcfsortcmd="vcf-sort $prefix.unsorted.vcf > $prefix.vcf"
    echo $vcfsortcmd
    eval $vcfsortcmd

done

echo "Working on indel files now"

IN_INDEL_FILES=*.unsorted.indel.vcf
for FILE in $IN_INDEL_FILES
do

    prefix=$(basename $FILE .unsorted.indel.vcf)

    vcfsortcmd="vcf-sort $FILE > $prefix.indel.vcf"
    echo $vcfsortcmd
    eval $vcfsortcmd

done

rm *.bak
rm *.unsorted.*
rm *.fixed.tsv
