
rule metaplasmidSPAdes:
    input:
        R1="host_rm_reads/splited_reads/{sample}.host_rm_R1.fq",
        R2="host_rm_reads/splited_reads/{sample}.host_rm_R2.fq",
        contig="assembly_per_sample/long_contigs/{sample}-1000.contigs.fa"
    params:
        outdir="spades_plasmid"
    threads: 24
    shell:"""
        python ~/mambaforge/envs/snakemake/bin/metaplasmidspades.py -o {params.outdir} -1 {input.R1} -2 {input.R2} --trusted-contigs {input.contig} --checkpoints all -t {threads}
    """
python ~/mambaforge/envs/snakemake/bin/metaplasmidspades.py \
-o spades_plasmid \
-1 host_rm_reads/splited_reads/LW_DA_1.host_rm_R1.fq \
-2 host_rm_reads/splited_reads/LW_DA_1.host_rm_R2.fq \
--trusted-contigs assembly_per_sample/long_contigs/LW_DA_1-1000.contigs.fa \
--checkpoints all -t 24