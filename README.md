# OICR_Task
This repository contains programming task assigned by OICR; it parses through a Variant calling format file and outputs the stats in JSON

In order to run it you will have to provide the sample name in the command line eg:- 'python OICR_task.py HG00096' and 'python OICR_task.py HG00097'.

Make sure to have the hg96-hg97.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf file in the same directory as the python script.

I used python library PyVCF to parse through the vcf file for this test. 

It can be installed using the command 'python -m pip install PyVCF' on the command line.  
