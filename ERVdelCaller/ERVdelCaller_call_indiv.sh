#!/usr/bin/bash
#per sample.
#cmd: bash ERVdelCaller_coh3.sh /PATH/TO/YOUR.config sampleID
#please use mincpus=8, mem=8GB


sample=$2
bam_suffix=""
. $1
#---
#ERVdelCaller=$ERVdelCaller
#gatk=$gatk
manta_install_path=$manta
#python2=$python2

ref_genome_fasta=$refgenome
if [[ ! -v input_bam_dir ]];then
	manta_input_bam=${EC_out_dir}/h_bam/${sample}_h.bam || echo "Directory of input .bam not found!"
else
	manta_input_bam=${input_bam_dir}/${sample}${bam_suffix}.bam
fi
my_analysis_dir=$EdC_wd
#---


manta_analysis_dir=${my_analysis_dir}/manta/${sample}
manta_log=${manta_analysis_dir}/logs/${sample}.no_dup_bed_Manta.log
mkdir -p ${manta_analysis_dir}
mkdir -p ${manta_analysis_dir}/logs
#****Manta
 #* Configuration Manta
 #* Single Diploid Sample Analysis
${manta_install_path}/bin/configManta.py \
--bam ${manta_input_bam} \
--referenceFasta ${ref_genome_fasta} \
--callRegions ${ERV_bed_search} \
--runDir ${manta_analysis_dir} \
1>${manta_log} 2>&1
if [ $? != 0 ];then exit;fi
/usr/bin/time -v ${python2} ${manta_analysis_dir}/runWorkflow.py -j 4 1>>${manta_log} 2>&1


###filter SVtype=DEL;SVlength<10000bp
mkdir -p ${my_analysis_dir}/filter_10000/${sample}
cd ${my_analysis_dir}/manta/${sample}/results/variants
cp diploidSV.vcf.gz ${my_analysis_dir}/filter_10000/${sample}/${sample}.vcf.gz

cd ${my_analysis_dir}/filter_10000/${sample}
bcftools view -i 'SVTYPE="DEL" && abs(INFO/SVLEN)<10000' ${sample}.vcf.gz -Oz > ${sample}.del.1000.vcf.gz
bcftools index -t ${sample}.del.1000.vcf.gz #/share/apps/gene/bcftools-1.14/bin/
#gunzip ${sample}.del.1000.vcf.gz
#grep -v '#' ${sample}.del.1000.vcf|awk '{print $1,$2,$3}' |sed 's/ /\t/g'> ${sample}.ERV.SV.chr.start.SVID
#grep -v '#' ${sample}.del.1000.vcf|awk '{print $8}'|awk -F ';' '{print $1}'|sed 's/END=//g'> ${sample}.SV_end
#grep -v '#' ${sample}.del.1000.vcf|awk '{print $8}'|awk -F ';' '{print $3}'|sed 's/SVLEN=-//g' > ${sample}.SV_length
bcftools view -H ${sample}.del.1000.vcf.gz|awk '{print $1,$2,$3}' |sed 's/ /\t/g'> ${sample}.ERV.SV.chr.start.SVID
bcftools view -H ${sample}.del.1000.vcf.gz|awk '{print $8}'|awk -F ';' '{print $1}'|sed 's/END=//g'> ${sample}.SV_end
bcftools view -H ${sample}.del.1000.vcf.gz|awk '{print $8}'|awk -F ';' '{print $3}'|sed 's/SVLEN=-//g' > ${sample}.SV_length
paste ${sample}.ERV.SV.chr.start.SVID ${sample}.SV_end ${sample}.SV_length > filter_10000.DEL.chr.start.SVID.end.length.bed
awk '{print $1,$2-1,$4-1,$3,$5}' filter_10000.DEL.chr.start.SVID.end.length.bed|sed 's/ /\t/g'> filter_10000.DEL.chr.start-1.end-1.SVID.length.bed
bedtools intersect -a ${ERV_bed_ref} -b filter_10000.DEL.chr.start-1.end-1.SVID.length.bed -wa -wb > ${sample}.table
bedtools intersect -a ${ERV_bed_ref} -b filter_10000.DEL.chr.start-1.end-1.SVID.length.bed -wa -c|awk '{print $NF}' > ${sample}.count.table
awk '{print $(NF-3)-$2}' ${sample}.table > left.distance
awk '{print $(NF-2)-$3}' ${sample}.table > right.distance
awk '{print $3-$2}' ${sample}.table > reference.interval.length
paste left.distance right.distance | awk '{abs1 = $1 < 0 ? -$1 : $1; abs2 = $2 < 0 ? -$2 : $2; print abs1 + abs2}' > sum_left_right
less ${sample}.table |awk '{print $(NF-4),$(NF-1),$(NF-3),$(NF),$4,$2}'|sed 's/ /\t/g' > ${sample}.intersect.table
paste ${sample}.intersect.table reference.interval.length left.distance right.distance|awk '{print $1,$2,$3+1,$4,$5,$6+1,$7,$8,$9}'|sed 's/ /\t/g' > ${sample}.all.table
echo -e "chromosome\tSV_ID\tSV_start_position\tSV_length\treference_interval_ID\treference_interval_start_position\treference_interval_length\tleft_distance\tright_distance" > ${sample}_ERVdelSV_lenflt_info.txt
cat ${sample}.all.table >> ${sample}_ERVdelSV_lenflt_info.txt
awk '{print $1,$3}' ${sample}.all.table |sed 's/ /\t/g' >chr.pos
bzgip ${sample}.del.1000.vcf
bcftools view -T chr.pos ${sample}.del.1000.vcf.gz -O z > ${sample}_ERVdelSV_lenflt.vcf.gz #这个数量对不上是因为同一个SV对应有不同的ERVid 导致这两者数量不一致，vcf里面的数量要较少一点
bcftools index -t ${sample}_ERVdelSV_lenflt.vcf.gz