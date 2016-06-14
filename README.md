Wheat TILLING Project
===============


The tetraploid tilling project homepage can be found [here](http://dubcovskylab.ucdavis.edu/wheat-tilling) and the hexaploid project homepage can be found [here](http://wheat-tilling.com). These data will be available to the public after publication, but early access can be requested by sending an email to <jdubcovsky@ucdavis.edu>, <havasquezgross@ucdavis.edu>, and <Cristobal.Uauy@jic.ac.uk>.


###Bioinformatics Pipeline Steps:
1. Trim sequencing data with sickle/scythe/seqqs
2. Run mapping step with bwa for genome or novoalign for capture design
3. Run PicardTool's MarkDuplicates to remove PCR and Optical Duplicates
4. Run MAPS pipeline
5. Convert MAPS output to VCF file for different confidence intervals
6. Run python script to split up the VCF files correctly
7. Run updated VCF files through python script to generate SQL insert statements to update BLAST; update SQL database with new sql file


For archive purposes, we have added our version of the MAPS pipeline used for these analyses.  If you would like to get the current version of MAPS pipeline, please visit the Comai Lab's project page [here](http://comailab.genomecenter.ucdavis.edu/index.php/MAPS).

###Folder Descriptions:
  * preprocessing: contains scripts to perform read trimming with sickle/scythe/seqqs, mapping trimmed reads, and then duplicate removal with PicardTool's MarkDuplicates to produce BAM files
  * maps_pipeline: contains programs for processing resulting BAM files with the MAPs pipeline
  * postprocessing: contains programs for identifying residual heterogeneity, mutations present in multiple individuals, downstream file conversion to VCF, and stats from the results
  * jbrowse_config: contains custom configuration for displaying VCF mutation information on JBrowse



###List of authors:
  * K.V. Krasileva<sup>1,2,3</sup>
  * H. Vasquez-Gross<sup>1</sup> 
  * T. Howell<sup>1</sup>
  * P. Bailey<sup>3</sup> 
  * F. Paraiso<sup>1</sup>
  * L. Clissold<sup>3</sup> 
  * J. Simmonds<sup>4</sup> 
  * R. Ramirez-Gonzalez<sup>3,4</sup> 
  * X. Wang<sup>1</sup>
  * P. Borrill<sup>4</sup> 
  * C. Fosker<sup>3</sup> 
  * S. Ayling<sup>3</sup> 
  * A. Phillips<sup>5</sup> 
  * C. Uauy<sup>4,7</sup>
  * J. Dubcovsky<sup>1,6,7</sup>


1. Dept. of Plant Sciences, University of California, Davis, CA 95616, USA 
2. The Sainsbury Laboratory, Norwich Research Park, Norwich, NR4 7UH United Kingdom 
3. The Genome Analysis Centre, Norwich Research Park, Norwich, NR4 7UH United Kingdom
4. John Innes Centre, Norwich Research Park, Norwich, NR4 7UH, United Kingdom
5. Rothamsted Research, Harpenden, Hertfordshire, AL5 2JQ, United Kingdom
6. Howard Hughes Medical Institute (HHMI), Chevy Chase, MD 20815, USA
7. The last two authors contributed equally to the direction of this project and are co-corresponding authors.

