Wheat TILLING Project
===============

#Bioinformatics Pipeline Steps:
1. Trim sequencing data with sickle/scythe/seqqs
2. Run mapping step with bwa for genome or novoalign for capture design
3. Run PicardTool's MarkDuplicates to remove PCR and Optical Duplicates
4. Run MAPS pipeline
5. Convert MAPS output to VCF file for different confidence intervals
6. Run python script to split up the VCF files correctly
7. Run updated VCF files through python script to generate SQL insert statements to update BLAST; update SQL database with new sql file


