trim_adapter="/cluster/home/yuma/miniconda2/share/trimmomatic-0.38-1/adapters/TruSeq3-PE-2.fa"


# 0. merge fq.gz
gz_1=`ls $SRR_dir|grep SRR.*_1.fastq.gz|sed "s#^#$SRR_dir\/#"`
gz_2=`ls $SRR_dir|grep SRR.*_2.fastq.gz|sed "s#^#$SRR_dir\/#"`

# output
fq_1=$SRR_dir/${taxon}_${acc}_1.fastq.gz
fq_2=$SRR_dir/${taxon}_${acc}_2.fastq.gz

cat $gz_1 > $fq_1
cat $gz_2 > $fq_2

# 1. fastqc_raw

fastqc --quiet -o ${datapath}/qc_raw ${fq_1}
fastqc --quiet -o ${datapath}/qc_raw ${fq_2}

# 2.trim
r1_p=${SRR_dir}/${taxon}_${acc}_R1_paired.fastq.gz
r2_p=${SRR_dir}/${taxon}_${acc}_R2_paired.fastq.gz
r1_s=${SRR_dir}/${taxon}_${acc}_R1_unpaired.fastq.gz
r2_s=${SRR_dir}/${taxon}_${acc}_R2_unpaired.fastq.gz

#shell:
trimmomatic PE -threads 3 -phred33 $fq_1 $fq_2 $r1_p $r1_s $r2_p $r2_s ILLUMINACLIP:${trim_adapter}:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30

# 3. fastqc clean
fastqc --quiet -o ${datapath}/qc_clean ${r1_p}
fastqc --quiet -o ${datapath}/qc_clean ${r2_p}

