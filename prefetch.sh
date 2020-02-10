source activate methylKit

while [ ! -f $SRR_dir/$SRR/${SRR}.sra ]
do
prefetch $SRR -O $SRR_dir --ascp-path "/cluster/home/yuma/miniconda2/envs/methylKit/bin/ascp|/cluster/home/yuma/miniconda2/pkgs/aspera-cli-3.9.1-0/etc/asperaweb_id_dsa.putty"
done
sleep 15s

if [ ! -f $SRR_dir/${SRR}_1.fastq.gz || ! -f $SRR_dir/${SRR}_2.fastq.gz ]
then
fastq-dump -I --split-3 -Q 33 --defline-seq "@\$sn" --defline-qual "+" $SRR_dir/$SRR/${SRR}.sra -O $SRR_dir --gzip
fi
