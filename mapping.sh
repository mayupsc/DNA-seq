bw_index_path="/cluster/sharedata/langzhaobo/pscngs/genome/SL3.0_bw_index/S_lycopersicum_chromosomes.3.00.fa"

#input
r1_p=${SRR_dir}/${taxon}_${acc}_R1_paired.fastq.gz
r2_p=${SRR_dir}/${taxon}_${acc}_R2_paired.fastq.gz

sam=${SRR_dir}/${taxon}_${acc}.sam
/cluster/apps/bwa/0.7.17/bwa mem -M -t $threads $bw_index_path $r1_p $r2_p -o $sam

