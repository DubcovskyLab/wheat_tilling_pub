Wheat TILLING Project
===============

These programs are the versions used for the MAPS pipeline.  Towards the end of the project, our MAPS analysis was IO bound. To ameliorate this, we modifed the scripts to allow for compression.

###Example pipeline:

```bash

python beta-run-mpileup.py -t 30 -r $REF_DIR/IWGSC_CSS_AB-TGAC_UCW_v1.fa -o Kronos_mpileup.txt -s /usr/bin/samtools --bamname rmdup.bam

python beta-mpileup-parser.py -t 30 -f Kronos_mpileup.txt

python beta-maps1-v2.py -f parsed_Kronos_mpileup.txt -t 30 -l 20 -c 10000 -o Kronos.mapspart1.defaults.txt

python maps-part2-v2.py -f $FILE -l $libcutoff -s $k -d $j -p $i -o $basename.mapspart2.HetMinCov${j}HetMinPer${i}HomMinCov${k}.tsv -m m

```

