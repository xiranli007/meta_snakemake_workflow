(SAMPLES,) = glob_wildcards("assembly_per_sample/long_contigs/{sample}-1000.contigs.fa")
rule all:
    input: expand( "plasflow_prediction/{sample}.plasflow_predictions.tsv", sample=SAMPLES)
rule plasflow:
    input:
        "assembly_per_sample/long_contigs/{sample}-1000.contigs.fa"
    output:
        "plasflow_prediction/{sample}.plasflow_predictions.tsv"
    shell:
        "python /home/xiran007/mambaforge/envs/plasflow/bin/PlasFlow.py --input {input} --output {output}" 


#  use updated plasflow for segmentation error Plasflow v 1.1 -- solved 08/20/2023 -xl