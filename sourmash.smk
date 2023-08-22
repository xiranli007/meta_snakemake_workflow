SAMPLES, = glob_wildcards("host_rm_reads/splited_reads/{sample}.host_rm_R1.fq")
rule all:
   input:expand("smash_files/kreports/{sample}.kreport.txt",sample=SAMPLES)
rule sketch:
   input:
       r1="host_rm_reads/splited_reads/{sample}.host_rm_R1.fq",
       r2="host_rm_reads/splited_reads/{sample}.host_rm_R2.fq"
   output: "smash_files/{sample}.sig.gz"
   shell:
       "sourmash sketch dna {input.r1} {input.r2} -o {output} --name {wildcards.sample}"
rule gather:
   input:
      sig="smash_files/{sample}.sig.gz",
      db="/group/ctbrowngrp/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k31.zip"
   output:"smash_files/{sample}.gather.csv"
   params:
      threshold=1
   shell:
      " sourmash gather {input.sig} {input.db} -o {output} --threshold-bp={params.threshold}"
   
rule taxonomy:
   input:
      gather="smash_files/{sample}.gather.csv",
      taxonomy="/group/ctbrowngrp/sourmash-db/gtdb-rs207/gtdb-rs207.taxonomy.csv"
   output:"smash_files/kreports/{sample}.kreport.txt"
   shell:
      "sourmash tax metagenome -g {input.gather} -t {input.taxonomy}  -F kreport -o {output}"