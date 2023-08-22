import os
path = "assembly_per_sample/contigs/"
path2 = "assembly_per_sample/long_contigs/"
path3 = "waafle/"
path4="waafle/results/"
if not os.path.isdir(path):
    os.mkdir(path)
if not os.path.isdir(path2):
    os.mkdir(path2)
(SAMPLES,) = glob_wildcards("host_rm_reads/splited_reads/{sample}.host_rm_R1.fq")
rule all:
    input:
        expand("waafle/results/{sample}-1000.lgt.tsv",sample=SAMPLES),
        expand("waafle/results/{sample}-1000.no_lgt.tsv", sample=SAMPLES),
        expand("waafle/results/{sample}-1000.unclassified.tsv",sample=SAMPLES)
rule gather_fasta:
    input: "assembly_per_sample/{sample}/{sample}.contigs.fa"
    output:"assembly_per_sample/contigs/{sample}.contigs.fa"
    shell:"""
        mv {input} {output}
    """
rule extract_long_contigs:
    input: "assembly_per_sample/contigs/{sample}.contigs.fa"
    output:"assembly_per_sample/long_contigs/{sample}-1000.contigs.fa"
    shell:"""
        bioawk -c fastx '{{if(length($seq) > 1000) {{print ">"$name; print $seq }}}}' {input} > {output}
    """
rule waafle_search:
    message:" Step 1 in the WAAFLE pipeline (MGE detection), for more info, refer to https://github.com/biobakery/biobakery/wiki/WAAFLE"
    input:"assembly_per_sample/long_contigs/{sample}-1000.contigs.fa"
    output:"waafle/{sample}.blastout"
    threads: 24
    shell:"""
        waafle_search {input} /home/xiran007/davis_water/waafle/waafle_db/waafledb --threads {threads} --out {output}
    """
rule waafle_genecaller:
    message: "Step 2 in the WAAFLE pipeline"
    input: "waafle/{sample}.blastout"
    output:"waafle/{sample}.gff"
    shell: """
        waafle_genecaller {input} --gff {output}
    """
rule waafle_orgscorer:
    message:"Step 3 in the WAAFLE pipeline"
    input:
       contigs = "assembly_per_sample/long_contigs/{sample}-1000.contigs.fa",
       blastout = "waafle/{sample}.blastout",
       gff = "waafle/{sample}.gff" ,
       tax = "/home/xiran007/davis_water/waafle/waafle_db/waafledb_taxonomy.tsv"
    output:
       "waafle/results/{sample}-1000.lgt.tsv",
       "waafle/results/{sample}-1000.no_lgt.tsv",
       "waafle/results/{sample}-1000.unclassified.tsv"
    params:
       outdir = "waafle/results"
    shell: """
       waafle_orgscorer {input.contigs} {input.blastout} {input.gff} {input.tax} --outdir {params.outdir}
    """

    # if you want to change the name of the waafle_orgscorer use --basename