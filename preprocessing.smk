(SAMPLES,) = glob_wildcards("renamed_raw_data/{sample,[^/]+}_R1.fastq.gz")
# [^/]+ to contrain glob_wildcards to go to subfolders
rule all:
    input: 
        expand("host_rm_reads/splited_reads/{sample}.host_rm_R1.fq",sample=SAMPLES),
        expand("host_rm_reads/splited_reads/{sample}.host_rm_R2.fq",sample=SAMPLES)

rule quality_trim:
    message:"Trims given paired-end reads with given parameters"
    input:
        R1 = "renamed_raw_data/{sample}_R1.fastq.gz",
        R2 = "renamed_raw_data/{sample}_R2.fastq.gz",
        adapters = "/home/xiran007/apps/Trimmomatic-0.39/adapters/TruSeq3-PE.fa", 
    output: 
        forwards_paired = "trimmed/{sample}_R1.paired.fastq.gz",
        forwards_unpaired = "trimmed/{sample}_R1.unpaired.fastq.gz",
        backwards_paired = "trimmed/{sample}_R2.paired.fastq.gz",
        backwards_unpaired = "trimmed/{sample}_R2.unpaired.fastq.gz",
    shell:
        """
        trimmomatic PE -phred33 {input.R1} {input.R2} {output.forwards_paired} {output.forwards_unpaired} {output.backwards_paired} {output.backwards_unpaired} ILLUMINACLIP:{input.adapters}:2:30:10:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36   
        """
rule remove_host:
    message: "map the trimmed reads to host genome for further host removal"
    input:
        R1 = "trimmed/{sample}_R1.paired.fastq.gz",
        R2 = "trimmed/{sample}_R2.paired.fastq.gz",
        ref = "host/catfish_genome.fna"
    output: 
        sam = "host_rm_reads/sam/{sample}.x.catfish_genome.sam",
        bam = "host_rm_reads/bam/{sample}.x.catfish_genome.bam",
    shell: """
        minimap2 -ax sr {input.ref} {input.R1} {input.R2} > {output.sam}
        samtools view -bS {output.sam} > {output.bam}
    """
rule sam_to_bam:
    message:"convert sam to bam to extract host genome"
    input:"host_rm_reads/bam/{sample}.x.catfish_genome.bam"
    output:
        unmapped_bam= "host_rm_reads/unmapped_bam/{sample}.x.catfish_genome_unmapped.bam",
        sorted_bam = "host_rm_reads/sorted_filtered_bam/{sample}.x.catfish_genome_unmapped.bam.sorted"
    shell: """
        samtools view -f 12 -F 256 {input} > {output.unmapped_bam}
        samtools sort {output.unmapped_bam} > {output.sorted_bam}
    """
rule bam_to_fastq:
    message: "convert sorted bam file to separate fastq file"
    input: "host_rm_reads/sorted_filtered_bam/{sample}.x.catfish_genome_unmapped.bam.sorted"
    output:
        R1="host_rm_reads/splited_reads/{sample}.host_rm_R1.fq",
        R2="host_rm_reads/splited_reads/{sample}.host_rm_R2.fq"
    shell: """
        samtools fastq {input} -1 {output.R1} -2 {output.R2}
    """
        