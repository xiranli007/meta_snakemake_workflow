SAMPLES,=glob_wildcards("amr_detection/{sample}.amr.alignment.sam")
rule all:
    input:expand("amr_detection/arg_contigs/{sample}.arg_contigs.fa", sample=SAMPLES)
    
rule extract_arg_reads:
    input:"amr_detection/{sample}.amr.alignment.sam"
    output:
          bam="amr_detection/alignment_bam/{sample}.amr.alignment.bam",
          reads="amr_detection/arg_reads/{sample}.amr.alignment.fq"
    shell:"""
        samtools view -bS {input} > {output.bam}
        samtools bam2fq -F 4 {output.bam} > {output.reads}
    """
rule map_arg_reads_x_contigs:
    input:
        arg_reads="amr_detection/arg_reads/{sample}.amr.alignment.fq",
        contigs="assembly_per_sample/contigs/{sample}.contigs.fa"
    output:
        sam="amr_detection/arg_reads.x.contigs/{sample}.arg_reads.x.contigs.sam",
        bam="amr_detection/arg_reads.x.contigs/{sample}.arg_reads.x.contigs.bam",
        sorted_bam="amr_detection/arg_reads.x.contigs/{sample}.arg_reads.x.contigs.bam.sorted"
    shell:"""
        minimap2 -ax sr {input.contigs} {input.arg_reads} -o {output.sam}
        samtools view -b -F 4 {output.sam} > {output.bam}
        samtools sort {output.bam} > {output.sorted_bam}
     """
rule extract_arg_contigs:
    input:
        sorted_bam="amr_detection/arg_reads.x.contigs/{sample}.arg_reads.x.contigs.bam.sorted",
        contigs="assembly_per_sample/contigs/{sample}.contigs.fa"
    output:
        arg_contigs_names="amr_detection/arg_contigs/{sample}.arg.x.contigs.txt",
        arg_contigs="amr_detection/arg_contigs/{sample}.arg_contigs.fa"
    shell:"""
       set +o pipefail
       samtools depth {input.sorted_bam}  | cut -f1 | sort | uniq > {output.arg_contigs_names}
       seqtk subseq {input.contigs} {output.arg_contigs_names} > {output.arg_contigs}
    """

    # Note: snakemake default bash pipefail, to deactivate that us set +o pipefail
    # it should be used in caution (i.e. you should test your pipe to make sure it is correct)
    # otherwise it will mak succeed also a broken command