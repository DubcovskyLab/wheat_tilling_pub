Wheat TILLING Project
===============

The tetraploid tilling project homepage can be found [here](http://dubcovskylab.ucdavis.edu/wheat-tilling) and the hexaploid project homepage can be found [here](http://wheat-tilling.com). These data will be available to the public after publication, but early access can be requested by sending an email to <jdubcovsky@ucdavis.edu>, <havasquezgross@ucdavis.edu>, and <Cristobal.Uauy@jic.ac.uk>.


##Bioinformatics Pipeline Steps:
1. Trim sequencing data with sickle/scythe/seqqs
2. Run mapping step with bwa for genome or novoalign for capture design
3. Run PicardTool's MarkDuplicates to remove PCR and Optical Duplicates
4. Run MAPS pipeline
5. Convert MAPS output to VCF file for different confidence intervals
6. Run python script to split up the VCF files correctly
7. Run updated VCF files through python script to generate SQL insert statements to update BLAST; update SQL database with new sql file


For archive purposes, we have added our version of the MAPS pipeline used for these analyses.  If you would like to get the current version of MAPS pipeline, please visit the Comai Lab's project page [here](http://comailab.genomecenter.ucdavis.edu/index.php/MAPS)

