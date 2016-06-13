Wheat TILLING Project
===============

Since our analyses on our cluster were IO bound, we modified the MAPS scripts to allow piping so we can use compression.

###Examples of how to use the compression scripts:

```bash

python beta-run-mpileup-compression.py -t 18 -r $REF_DIR/IWGSC_CSS_AB-TGAC_UCW_v1.fa -o - -s /usr/bin/samtools --bamname .bam | tee >(wc -l > Kronos.mpileup.len) | pigz --best -p 28 > Kronos_mpileup.txt.gz

python beta-mpileup-parser-compression.py -t 18 -f Kronos_mpileup.txt.gz -l Kronos.mpileup.len | tee >(wc -l > parsed_Kronos.mpileup.len) | pigz --best -p 18 > parsed_Kronos_mpileup.txt.gz"

python beta-maps1-v2-compression.py -f parsed_Kronos_mpileup.txt.gz -t 18 -L parsed_Kronos.mpileup.len -l 20 -C 10000 -o Kronos.mapspart1.defaults.txt

```

