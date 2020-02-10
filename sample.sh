project=dgtang
part=acc_101_200
datapath=/cluster/sharedata/langzhaobo/pscngs/$project/$part
if [ ! -d $datapath ]
then
mkdir -p $datapath
fi 
#mkdir -p $datapath/qc_raw $datapath/qc_clean
SRR_list=/cluster/home/yuma/my_github/DNA-Seq/samples/part_tomato_${part}.txt
last_index=`wc -l ${SRR_list}|awk '{print $1}' -`
prefetch_script=/cluster/home/yuma/my_github/DNA-Seq/script/prefetch.sh
fastqc_trim_script=/cluster/home/yuma/my_github/DNA-Seq/script/fastqc_and_trim.sh
mapping_script=/cluster/home/yuma/my_github/DNA-Seq/script/mapping.sh
gatk_callsnp_script=/cluster/home/yuma/my_github/DNA-Seq/script/gatk_callsnp.sh
fastqc_Rscript=/cluster/home/yuma/my_github/DNA-Seq/script/reads_stat_by_fastqcr.R
## about gatk reference : gatk CreateSequenceDictionary -R /cluster/sharedata/langzhaobo/pscngs/genome/S_lycopersicum_chromosomes.3.00.fa -O /cluster/sharedata/langzhaobo/pscngs/genome/S_lycopersicum_chromosomes.3.00.dict

#for id in `seq 2 $last_index`
#for id in `seq 2 10`
for id in `seq 11 $last_index`

#for sample in SRR7647267
do
taxon=`awk 'BEGIN{FS="\t"}NR=="'$id'"{print $1}' $SRR_list`
acc=`awk 'BEGIN{FS="\t"}NR=="'$id'"{print $2}' $SRR_list`
SRRs=`awk 'BEGIN{FS="\t"}NR=="'$id'"{print $3}' $SRR_list`

#head="#PBS -N ${taxon}_${acc}_${SRR}\\n#PBS -q q2\\n#PBS -l nodes=node33:ppn=1\\n#PBS -j oe\\n"
SRR_dir=$datapath/${taxon}_${acc}
if [ ! -d $SRR_dir ]
then
mkdir -p $SRR_dir
fi

if [ ! -f $datapath/qc_raw ]
then
mkdir -p $datapath/qc_raw
fi

if [ ! -f $datapath/qc_clean ]
then
mkdir -p $datapath/qc_clean
fi

## data download 
SRR_array=($(echo $SRRs | tr "," "\n"))
SRR_jids_arr=()
for SRR in "${SRR_array[@]}"
do
jid=$(qsub $prefetch_script -N $prefetch_${SRR} -q q2 -l nodes=1:ppn=1 -e $SRR_dir/${part}.prefetch.err.log -o $SRR_dir/${part}.prefetch.log -v SRR=$SRR,SRR_dir=$SRR_dir )
SRR_jids_arr+=($jid)
done
#echo ${SRR_jids_arr[@]}

## fastqc and trim
fastqc_and_trim_jid=$(qsub $fastqc_trim_script -N fastqc_${taxon}_${acc} -q q2 -l nodes=1:ppn=1 -e $SRR_dir/${part}.fastqc.trim.err.log -o $SRR_dir/${part}.fastqc.trim.log -W depend=afterok`printf ':%s' ${SRR_jids_arr[@]}` -v SRR_dir=$SRR_dir,taxon=$taxon,acc=$acc,datapath=$datapath)

## mapping
threads=1
mapping_jid=$(qsub $mapping_script  -N mapping_${taxon}_${acc} -q q2 -l nodes=1:ppn=$threads -e $SRR_dir/${part}.mapping.err.log -o $SRR_dir/${part}.mapping.log -W depend=afterok:$fastqc_and_trim_jid  -v SRR_dir=$SRR_dir,taxon=$taxon,acc=$acc,datapath=$datapath,threads=$threads)

## gatk callsnp
gatk_snp_jid=$(qsub $gatk_callsnp_script  -N snp_${taxon}_${acc} -q q2 -l nodes=1:ppn=1 -e $SRR_dir/${part}.snp.err.log -o $SRR_dir/${part}.snp.log -W depend=afterok:$mapping_jid -v SRR_dir=$SRR_dir,taxon=$taxon,acc=$acc,datapath=$datapath)

done
###############################################   Part3 reads count in various steps   #################################
#echo -e "/usr/bin/Rscript $fastqc_Rscript $datapath $part\\n"> $datapath/${part}.fastqc.pbs
#fastqc_stat_jid=$(qsub $datapath/${part}.fastqc.pbs -N ${part}_fastqc -q q2 -l nodes=1:ppn=1 -W depend=afterok:$snp_jid -e $datapath/${part}.fastqc.err.log -o $datapath/${part}.fastqc.log)




