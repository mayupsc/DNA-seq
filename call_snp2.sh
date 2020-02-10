#SRR_dir=$1
#taxon=$2
#acc=$3
#datapath=$4

trim_adapter="/cluster/home/yuma/miniconda2/share/trimmomatic-0.38-1/adapters/TruSeq3-PE-2.fa"
bw_index_path="/cluster/sharedata/langzhaobo/pscngs/genome/SL3.0_bw_index/S_lycopersicum_chromosomes.3.00.fa"
gatk_ref="/cluster/sharedata/langzhaobo/pscngs/genome/S_lycopersicum_chromosomes.3.00.fa"


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
# 4. mapping

#output
sam=${SRR_dir}/${taxon}_${acc}.sam
#bwa_map="bwa mem -M  -t 8 -R '@RG\\\tID:${SRR}\\\tLB:${acc}\\\tSM:${taxon}\\\tPL:ILLUMINA' ${bw_index_path} $r1_p $r2_p -o $sam\\n"
/cluster/apps/bwa/0.7.17/bwa mem -M -t 1 ${bw_index_path} $r1_p $r2_p -o $sam
echo -e "bwa mem -M -t 1 ${bw_index_path} $r1_p $r2_p -o $sam"

# 5. reads filter (q20;remove supplementary reads;keep proper paired reads)
bam=${SRR_dir}/${taxon}_${acc}.q1.uniq.sort.bam
grep -v XA:Z $sam|grep -v SA:Z| samtools view -Shub -q 1 -F 2048 -F 256 -f 0x02 | samtools sort - -T ${SRR_dir}/${taxon}_${acc} -o ${bam}

# 6. deduplication

#rmdup_bam=${SRR_dir}/${SRR}.q1.uniq.sort.rmdup.bam
#rmdup="samtools rmdup $bam $rmdup_bam\\n"

# 7. index
#bam_bai=${SRR_dir}/${SRR}.q1.uniq.sort.rmdup.bam.bai
#bam_index="samtools index $rmdup_bam $bam_bai\\n"


# 6. gatk markdup
gatk_mark_bam=${SRR_dir}/${taxon}_${acc}.q1.uniq.sort.markdup.bam
gatk_metrics=${SRR_dir}/${taxon}_${acc}.q1.uniq.sort.markdup.metrics.txt
/cluster/apps/gatk/4.0.3.0-4/gatk MarkDuplicates -I $bam -O $gatk_mark_bam -M $gatk_metrics

# next step is Base quality score recalibration, which needs the species snp database , skip this step 

# 7. add ReadGroupInfo

bam_w_rdinfo=${SRR_dir}/${taxon}_${acc}.q1.uniq.sort.markdup.rdinfo.bam
/cluster/apps/gatk/4.0.3.0-4/gatk AddOrReplaceReadGroups -I $gatk_mark_bam -O $bam_w_rdinfo --RGLB=${taxon} --RGPL=illumina --RGPU=${acc} --RGSM=${taxon}_${acc}
# 8. bam index
bam_idx=$SRR_dir/${taxon}_${acc}.q1.uniq.sort.markdup.rdinfo.bam.bai
samtools index $bam_w_rdinfo

# 9. gatk HaplotypeCaller
vcf=$SRR_dir/${taxon}_${acc}.vcf
/cluster/apps/gatk/4.0.3.0-4/gatk HaplotypeCaller -R $gatk_ref -I $bam_w_rdinfo -O $vcf

# 9-2. gatk HaplotypeCaller to generate gvcf
gvcf=$SRR_dir/${taxon}_${acc}.gvcf
/cluster/apps/gatk/4.0.3.0-4/gatk HaplotypeCaller --emit-ref-confidence GVCF -R $gatk_ref -I $bam_w_rdinfo -O $gvcf

# 10. extract SNPs
snp_vcf=$SRR_dir/${taxon}_${acc}.snps.vcf
gatk SelectVariants -V $vcf --select-type-to-include SNP -O $snp_vcf

# 11. extract Indels
indel_vcf=$SRR_dir/${taxon}_${acc}.indels.vcf
gatk SelectVariants -V $vcf --select-type-to-include INDEL -O $indel_vcf

# 12. filter SNPs
snp_filter_vcf=$SRR_dir/${taxon}_${acc}.snps.filter.vcf
gatk VariantFiltration -R $gatk_ref -V $snp_vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 ||ReadPosRankSum < -8.0" --filter-name "SNP_FILTER" -O $snp_filter_vcf

# 13. filter indels
indel_filter_vcf=$SRR_dir/${taxon}_${acc}.indel.filter.vcf
gatk VariantFiltration -R $gatk_ref -V $indel_vcf --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -20.0" --filter-name "INDEL_FILTER" -O $indel_filter_vcf

# 14. merge snp/indels
#merge_vcfs=$SRR_dir/${SRR}.merge.vcf
