# for direct your self to the folder you would like to generate the .sh script and run the code below on terminal
path="/home/xiran007/davis_water/host_rm_reads/splited_reads"
for R1 in $path/*.host_rm_R1.fq
do
   base="$(basename "$R1" .host_rm_R1.fq)"
   R2=$path/*.host_rm_R2.fq
   cat > run-$base.sh <<EOF
   source ~/.bashrc
   conda activate snakemake
   sourmash sketch dna $R1 $R2 -o "$base".sig.gz --name "$base"
   sourmash gather "$base".sig.gz /group/ctbrowngrp/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k31.zip -o "$base".gather.csv --threshold-bp=1
   sourmash tax metagenome -g "$base".gather.csv -t /group/ctbrowngrp/sourmash-db/gtdb-rs207/gtdb-rs207.taxonomy.csv  -F kreport > "$base".kreport.txt
EOF
echo bash run-$base.sh
done > smash_list.txt