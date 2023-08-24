import os
path1 = "amr_detection/arg_contigs_tax"
path2 = "amr_detection/arg_contigs_tax/kreports"
try:
   os.makedirs(path1)
   os.makedirs(path2)
except FileExistsError:
   pass
SAMPLES, = glob_wildcards("amr_detection/arg_contigs/{sample}.arg_contigs.fa")
rule all:
   input:expand("amr_detection/arg_contigs_tax/kreports/{sample}.kreport.txt",sample=SAMPLES)
rule sketch:
   input:"amr_detection/arg_contigs/{sample}.arg_contigs.fa"
   output:"amr_detection/arg_contigs_tax/{sample}.sig.gz"
   shell:
       "sourmash sketch dna {input} -o {output} --name {wildcards.sample}"
rule gather:
   input:
      sig="amr_detection/arg_contigs_tax/{sample}.sig.gz",
      db="/group/ctbrowngrp/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k31.zip"
   output:"amr_detection/arg_contigs_tax/{sample}.gather.csv"
   params:
      threshold=1
   shell:
      "sourmash gather {input.sig} {input.db} -o {output} --threshold-bp={params.threshold}"
   
rule taxonomy:
   input:
      gather="amr_detection/arg_contigs_tax/{sample}.gather.csv",
      taxonomy="/group/ctbrowngrp/sourmash-db/gtdb-rs207/gtdb-rs207.taxonomy.csv"
   output:"amr_detection/arg_contigs_tax/kreports/{sample}.kreport.txt"
   params:"amr_detection/arg_contigs_tax/kreports/{sample}"
   shell:
      "sourmash tax metagenome -g {input.gather} -t {input.taxonomy}  -F kreport -o {params}"