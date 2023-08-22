(SAMPLES,) = glob_wildcards("host_rm_reads/splited_reads/{sample}.host_rm_R1.fq")
rule all:
     input:"amr_detection/AMR_analytic_matrix.csv"
rule align_to_amr:
     input:
        "/home/xiran007/AMRplusplus/data/amr/megares_database_v3.00.fasta.amb",
        "/home/xiran007/AMRplusplus/data/amr/megares_database_v3.00.fasta.ann",
        "/home/xiran007/AMRplusplus/data/amr/megares_database_v3.00.fasta.bwt",
        "/home/xiran007/AMRplusplus/data/amr/megares_database_v3.00.fasta.pac",
        "/home/xiran007/AMRplusplus/data/amr/megares_database_v3.00.fasta.sa",
        amr = "/home/xiran007/AMRplusplus/data/amr/megares_database_v3.00.fasta",
        R1= "host_rm_reads/splited_reads/{sample}.host_rm_R1.fq",
        R2= "host_rm_reads/splited_reads/{sample}.host_rm_R2.fq"
     output:
        "amr_detection/{sample}.amr.alignment.sam"
     params:
        rg = r"@RG\tID:${sample}\tSM:${sample}"
     threads:24
     shell: """
        bwa mem -t {threads} -R '{params.rg}' {input.amr} {input.R1} {input.R2} > {output}
     """

rule run_resistome:
     input:
        sam = "amr_detection/{sample}.amr.alignment.sam",
        amr = "/home/xiran007/AMRplusplus/data/amr/megares_database_v3.00.fasta",
        annotation = "/home/xiran007/AMRplusplus/data/amr/megares_annotations_v3.00.csv"
     output:
        gene_fp = "amr_detection/{sample}.gene.tsv",
        group_fp = "amr_detection/{sample}.group.tsv",
        mech_fp = "amr_detection/{sample}.mechanism.tsv",
        class_fp = "amr_detection/{sample}.class.tsv",
        type_fp = "amr_detection/{sample}.type.tsv"
     envmodules: "python/3.11"
     shell:
        "/home/xiran007/AMRplusplus/bin/resistome "
        "-ref_fp {input.amr} "
        "-annot_fp {input.annotation} "
        "-sam_fp {input.sam} "
        "-gene_fp {output.gene_fp} "
        "-group_fp {output.group_fp} "
        "-mech_fp {output.mech_fp} "
        "-class_fp {output.class_fp} "
        "-type_fp {output.type_fp} "
        "-t 80 "

rule resistome_results:
     input:
        expand("amr_detection/{sample}.gene.tsv", sample = SAMPLES)
     output:"amr_detection/AMR_analytic_matrix.csv"
     envmodules: "python/3.11"
     shell:
        "/home/xiran007/AMRplusplus/bin/amr_long_to_wide.py -i {input} -o {output}"
