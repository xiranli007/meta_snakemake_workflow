SAMPLES, = glob_wildcards("host_rm_reads/splited_reads/{sample}.host_rm_R2.fq")
rule all:
    input: 
        expand("amr_detection/metaxa/{sample}_metaxa_out.archaea.fasta",sample=SAMPLES),
        expand("amr_detection/metaxa/{sample}_metaxa_out.bacteria.fasta",sample=SAMPLES),
        expand("amr_detection/metaxa/{sample}_metaxa_out.chloroplast.fasta",sample=SAMPLES),
        expand("amr_detection/metaxa/{sample}_metaxa_out.eukaryota.fasta",sample=SAMPLES),
        expand("amr_detection/metaxa/{sample}_metaxa_out.extraction.fasta",sample=SAMPLES),
        expand("amr_detection/metaxa/{sample}_metaxa_out.mitochondria.fasta",sample=SAMPLES),
        expand("amr_detection/metaxa/{sample}_metaxa_out.taxonomy.txt",sample=SAMPLES),
        expand("amr_detection/metaxa/{sample}_metaxa_out.uncertain.fasta",sample=SAMPLES),
        expand("amr_detection/metaxa/{sample}_metaxa_out.summary.txt", sample=SAMPLES)
rule metaxa:
    input: 
        R1 = "host_rm_reads/splited_reads/{sample}.host_rm_R1.fq",
        R2 = "host_rm_reads/splited_reads/{sample}.host_rm_R2.fq",
    output:
        "amr_detection/metaxa/{sample}_metaxa_out.archaea.fasta",
        "amr_detection/metaxa/{sample}_metaxa_out.bacteria.fasta",
        "amr_detection/metaxa/{sample}_metaxa_out.chloroplast.fasta",
        "amr_detection/metaxa/{sample}_metaxa_out.eukaryota.fasta",
        "amr_detection/metaxa/{sample}_metaxa_out.extraction.fasta",
        "amr_detection/metaxa/{sample}_metaxa_out.mitochondria.fasta",
        "amr_detection/metaxa/{sample}_metaxa_out.taxonomy.txt",
        "amr_detection/metaxa/{sample}_metaxa_out.uncertain.fasta",
        "amr_detection/metaxa/{sample}_metaxa_out.summary.txt"
    params:
        prefix = "amr_detection/metaxa/{sample}_metaxa_out"
    shell: """
         metaxa2 -o {params.prefix} -1 {input.R1} -2 {input.R2} --cpu 24 -f fastq --mode metagenome --plus T -t b --graphical F
    """