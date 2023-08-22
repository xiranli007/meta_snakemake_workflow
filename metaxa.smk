SAMPLES, = glob_wildcards("host_rm_reads/splited_reads/{sample}.host_rm_R2.fq")
rule all:
    input: expand("amr_detection/metaxa/{sample}_metaxa_out", sample =SAMPLES)
rule metaxa:
    input: 
        R1 = "host_rm_reads/splited_reads/{sample}.host_rm_R1.fq",
        R2 = "host_rm_reads/splited_reads/{sample}.host_rm_R2.fq",
    output:"amr_detection/metaxa/{sample}_metaxa_out"
    shell: """
         metaxa2 -o {output} -1 {input.R1} -2 {input.R2} --cpu 24 -f fastq --mode metagenome --plus T -t b
    """