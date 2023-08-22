(SAMPLES,) = glob_wildcards("host_rm_reads/splited_reads/{sample}.host_rm_R1.fq")
rule all:
    input: 
        expand("assembly_per_sample/{sample}",sample= SAMPLES), 
rule megahit:
    '''
        Assembling fastq files using megahit.
        All files created by megahit are stored in a temporary folder,
        and only the fasta file is kept for later analysis.
    '''
    input: 
        R1 = "host_rm_reads/splited_reads/{sample}.host_rm_R1.fq",
        R2 = "host_rm_reads/splited_reads/{sample}.host_rm_R2.fq",
    output: directory ("assembly_per_sample/{sample}")
    threads: 24
    shell: """
        megahit -1 {input.R1} -2 {input.R2} -t {threads} -m 0.9 -o {output} --out-prefix {wildcards.sample}
    """