#!/bin/bash
#cmd: bash ERVdelCaller_merge.sh /PATH/TO/YOUR.config
#please use mincpus= n_sample if ≤20 else 20, mem= n_sample GB if ≤20GB else 20GB.


. $1
#---
ref_genome_fasta=$refgenome

if [[ ! -v input_bam_dir ]];then
	input_bam_dir=${EC_out_dir}/h_bam
	bam_suffix="_h"
fi

my_analysis_dir=$EdC_wd

cohortName=$cohort
#---



mkdir -p ${my_analysis_dir}/manta/manta_logs
mkdir -p ${my_analysis_dir}/merge_filter
cd ${my_analysis_dir}/merge_filter

exec 1> ${my_analysis_dir}/merge_filter/merge_filter.log
exec 2>&1

mv -f ${my_analysis_dir}/*no_dup_bed_Manta.log ${my_analysis_dir}/manta/manta_logs


n=0
>vcf_list
>bam_list
for vcf_in in ${my_analysis_dir}/filter_10000/*/*_ERVdelSV_lenflt.vcf.gz
do
if [ -e $vcf_in ];then
	ls ${vcf_in} >> vcf_list
	n=$((n+1))
else
	bgzip ${vcf_in/.gz/} && bcftools index -t ${$vcf_in} && ls ${vcf_in} >> vcf_list
	n=$((n+1))
fi
id0=$(basename ${vcf_in})
id=${id0/_ERVdelSV_lenflt.vcf.gz/}
ls ${input_bam_dir}/${id}${bam_suffix}\.bam >> bam_list
done


if [[ $n -lt 20 ]];then
	nt=$n
else
	nt=20
fi


#Get the list of merged positions
echo -e "\n###############svimmer begins###############\n"
python $svimmer --output ${cohortName}.merge_unsorted.vcf vcf_list chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
bcftools sort -o ${cohortName}.merge.vcf ${cohortName}.merge_unsorted.vcf && rm ${cohortName}.merge_unsorted.vcf
echo -e "\n###############svimmer ends###############\n"


#Get the info of merged positions
grep -v '#' ${cohortName}.merge.vcf|awk '{print $1,$2,$3}'|sed 's/ /\t/g'> ERV.SV.chr.start.SVID
grep -v '#' ${cohortName}.merge.vcf|awk '{print $8}'|awk -F ';' '{print $1}'|sed 's/END=//g'> SV_end
grep -v '#' ${cohortName}.merge.vcf|awk '{print $8}'|awk -F ';' '{print $3}'|sed 's/SVLEN=-//g'> SV_length
paste ERV.SV.chr.start.SVID SV_end SV_length > filter_10000.DEL.chr.start.SVID.end.length.bed
awk '{print $1,$2-1,$4-1,$3,$5}' filter_10000.DEL.chr.start.SVID.end.length.bed |sed 's/ /\t/g'> DEL.chr.start-1.end-1.length.bed
bedtools intersect -a ${ERV_bed_ref} -b DEL.chr.start-1.end-1.length.bed -wa -wb > bedtools.hg38_HERVd_entities_ERV1.DEL.chr.start-1.end-1.length.table
bedtools intersect -a ${ERV_bed_ref} -b DEL.chr.start-1.end-1.length.bed -wa -c|awk '{print $NF}' > bedtools.hg38_HERVd_entities_ERV1.count.table
less bedtools.hg38_HERVd_entities_ERV1.DEL.chr.start-1.end-1.length.table|awk '{print $2,$3,$(NF-3),$(NF-2),$(NF-4),$(NF-1),$4,$NF}'|sed 's/ /\t/g' > reference_interval

python $ERVdelCaller/scripts/calculate_overlap_percentage.py $my_analysis_dir
awk '{print $3-$1}' $PWD/reference_interval_filtered_results.txt > left.distance
awk '{print $4-$2}' $PWD/reference_interval_filtered_results.txt > right.distance
awk '{print $2 - $1}' $PWD/reference_interval_filtered_results.txt > reference.interval.length
paste left.distance right.distance | awk '{abs1 = $1 < 0 ? -$1 : $1; abs2 = $2 < 0 ? -$2 : $2; print abs1 + abs2}' > sum_left_right
paste $PWD/reference_interval_filtered_results.txt left.distance right.distance reference.interval.length > all.bed

awk 'function abs(x) {return x < 0 ? -x : x} {if (abs($(NF-2)) < 1000) print}' all.bed > all.bed1
awk 'function abs(x) {return x < 0 ? -x : x} {if (abs($(NF-1)) < 1000) print}' all.bed1 > all.bed2
awk '{print $5,$3,$(NF-3),$7,$1,$NF,$(NF-2),$(NF-1)}' all.bed2|sed 's/ /\t/g' > ${cohortName}.txt
awk '{print $1,$2+1,$3}' ${cohortName}.txt|sed 's/ /_/g' > SV_ID
paste SV_ID ${cohortName}.txt|awk '{print $2,$1,$3+1,$4,$5,$6+1,$7,$8,$9}'|sed 's/ /\t/g' > ${cohortName}.txt1

awk '{print $1,$3}' ${cohortName}.txt1 |sed 's/ /\t/g' >chr.pos
bcftools view -T chr.pos ${cohortName}.merge.vcf -O v > ${cohortName}_ERVdelSV_noGT.vcf #这个数量对不上是因为同一个SV对应有不同的ERVid 导致这两者数量不一致，vcf里面的数量要较少一点


#Genotype the merged positions using graphtyper
bgzip -c ${cohortName}_ERVdelSV_noGT.vcf > ${cohortName}_ERVdelSV_noGT.vcf.gz
bcftools index -t ${cohortName}_ERVdelSV_noGT.vcf.gz
#for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY;do
#	echo $i >> region_list
#done
awk '{bg=(($3));ed=(($3));if(bg<0){bg=0};print $1":"bg"-"ed}' ${cohortName}.txt1 > region_list

echo -e "\n###############graphtyper begins###############\n"
$graphtyper genotype_sv $ref_genome_fasta ${cohortName}_ERVdelSV_noGT.vcf.gz --sams=bam_list --region_file=region_list --threads=$nt --output=ERVdelSV_results && echo "genotyping done"
ls $PWD/ERVdelSV_results/*/*.vcf.gz > concat_vcf_list
bcftools concat --naive --file-list concat_vcf_list -o ${cohortName}_ERVdelSV_unsorted.vcf
bcftools sort -o ${cohortName}_ERVdelSV.vcf.gz ${cohortName}_ERVdelSV_unsorted.vcf && rm ${cohortName}_ERVdelSV_unsorted.vcf
mv ${cohortName}_ERVdelSV.vcf.gz ${cohortName}_ERVdelSV_raw.vcf.gz
tabix -f -p vcf ${cohortName}_ERVdelSV_raw.vcf.gz
bcftools filter -i 'INFO/SVMODEL="AGGREGATED" && FILTER="PASS"' -o ${cohortName}_ERVdelSV.vcf.gz ${cohortName}_ERVdelSV_raw.vcf.gz
tabix -f -p vcf ${cohortName}_ERVdelSV.vcf.gz
echo -e "\n###############graphtyper ends###############\n"


echo -e "\n###############vcftools get frequency begins###############\n"
vcftools --gzvcf ${cohortName}_ERVdelSV.vcf.gz --freq --out ${cohortName}_ERVdelSV_freq
awk 'NR==1{sub(/:/,"\t",$5);gsub(/[\{\}]/,"",$5);allel2=$5;sub(/ALLELE/,"alt_ALLELE",allel2);sub(/FREQ/,"alt_FREQ",allel2);print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"allel2; exit}' ${cohortName}_ERVdelSV_freq.frq > ${cohortName}_ERVdelSV_freq.txt;\
awk 'NR>1{sub(/:/,"\t",$5);gsub(/[\<\>]/,"",$6);split($6,arr,":");arrlen=length(arr);new6=arr[2];i=3;while (i<arrlen){new6=new6":"arr[i];i++};print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"new6"\t"arr[arrlen]}' ${cohortName}_ERVdelSV_freq.frq >> ${cohortName}_ERVdelSV_freq.txt
echo "SVlen" > ${cohortName}_ERVdelSV_len
egrep -o "(SVSIZE=)[0-9]+" ${cohortName}_ERVdelSV_freq.frq|egrep -o "[0-9]+" >> ${cohortName}_ERVdelSV_len
paste ${cohortName}_ERVdelSV_len ${cohortName}_ERVdelSV_freq.txt|awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$2"_"$3"_"$1}' > ${cohortName}_ERVdelSV_freq_svid.txt
echo -e "\n###############vcftools get frequency ends###############\n"


#Get info of ERVdelSV positions
awk 'NR>1{if(FNR==NR){a[$7]=$6;next};if($2 in a){print $0"\t"a[$2]}}' ${cohortName}_ERVdelSV_freq_svid.txt ${cohortName}.txt1 > ${cohortName}.txt #&& rm ${cohortName}.txt1
echo -e "chromosome\tSV_ID\tSV_start_position\tSV_length\treference_interval_ID\treference_interval_start_position\treference_interval_length\tleft_distance\tright_distance\talt_freq" > ${cohortName}_ERVdelSV_info.txt
cat ${cohortName}.txt >> ${cohortName}_ERVdelSV_info.txt


#modify ID in vcf
mv ${cohortName}_ERVdelSV.vcf.gz ${cohortName}_ERVdelSV_oriid.vcf.gz
mv ${cohortName}_ERVdelSV.vcf.gz.tbi ${cohortName}_ERVdelSV_oriid.vcf.gz.tbi
bcftools view -H ${cohortName}_ERVdelSV_oriid.vcf.gz > ${cohortName}_ERVdelSV_oriid.vcf_records
awk 'NR>1{if(FNR==NR){a[((FNR-1))]=$7;next};$3=a[FNR];gsub(" ","\t", $0);print $0}'\
 ${cohortName}_ERVdelSV_freq_svid.txt ${cohortName}_ERVdelSV_oriid.vcf_records > ${cohortName}_ERVdelSV.vcf_records

bcftools view -h ${cohortName}_ERVdelSV_oriid.vcf.gz > ${cohortName}_ERVdelSV.vcf
cat ${cohortName}_ERVdelSV.vcf_records >> ${cohortName}_ERVdelSV.vcf


#cleaning
#rm -f vcf_list
mkdir -p backup
for i in $(ls|egrep -v "${cohortName}_ERVdelSV_info.txt|${cohortName}_ERVdelSV.vcf|ERVdelSV_results|merge_filter.log|backup");do mv $i backup/;done
cd ${my_analysis_dir}

#info and filter
$python3 ${ERVdelCaller}/scripts/ERVdelCaller_postcall_cohort.py $1
