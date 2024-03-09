#!/usr/bin/bash
#cmd: bash s1m2.sh /PATH/TO/YOUR.config


configPath=$1
. $configPath


cd $postEC_wd/reads_support_noIns

#get gathered .nonTE and do nonTE GT
cat *.nonTE |awk '$0 !~ /\/\d:[^\d]/{print $0}' > ${cohort}.nonTE_gather && \
perl ${ERVcaller}/Scripts/Distinguish_nonTE_from_missing_genotype.pl -n ${cohort}.nonTE_gather -v ../${cohort}_ERVcaller_preflt_oriPosID.vcf -o ../${cohort}_ERVcaller_preflt_readSupportGT.vcf

cd $postEC_wd/reads_support_noIns_virt

#get gathered .nonTE virt and do nonTE GT
cat *.nonTE |awk '$0 !~ /\/\d:[^\d]/{print $0}' > ${cohort}.nonTE_gather && \
perl ${ERVcaller}/Scripts/Distinguish_nonTE_from_missing_genotype.pl -n ${cohort}.nonTE_gather -v ../${cohort}_ERVcaller_preflt_virtualPos_for_readSupportInfo.vcf -o ../${cohort}_ERVcaller_preflt_virtualPos_readSupportGT.vcf

cd $postEC_wd

#concatenate gathered readSupportGT.vcf files
infile_feature=${cohort}_ERVcaller_preflt
grep ":NA" ${infile_feature}_readSupportGT.vcf|egrep -o "chr[0-9XY]+:[0-9]+:\S+" > ${infile_feature}_readSupportGT_extremGQdel_ID_list.txt
grep ":NA" ${infile_feature}_virtualPos_readSupportGT.vcf|egrep -o "chr[0-9XY]+:[0-9]+:\S+" > ${infile_feature}_virtualPos_readSupportGT_extremGQdel_ID_list.txt

read -a alen1 <<< $(wc -l ${infile_feature}_readSupportGT_extremGQdel_ID_list.txt);len1=${alen1[0]}
if [ $len1 -gt 0 ];then
	mv ${infile_feature}_readSupportGT.vcf ${infile_feature}_readSupportGT_ori_before_extremGQlineDel.vcf && \
	awk '{if($1 ~ /^chr/){gsub("NA","0",$0)};print $0}' ${infile_feature}_readSupportGT_ori_before_extremGQlineDel.vcf | \
	bcftools filter -e ID=@${infile_feature}_readSupportGT_extremGQdel_ID_list.txt -o ${infile_feature}_readSupportGT.vcf
fi
read -a alen2 <<< $(wc -l ${infile_feature}_virtualPos_readSupportGT_extremGQdel_ID_list.txt);len2=${alen2[0]}
if [ $len2 -gt 0 ];then
	mv ${infile_feature}_virtualPos_readSupportGT.vcf ${infile_feature}_virtualPos_readSupportGT_ori_before_extremGQlineDel.vcf && \
	awk '{if($1 ~ /^chr/){gsub("NA","0",$0)};print $0}' ${infile_feature}_virtualPos_readSupportGT_ori_before_extremGQlineDel.vcf | \
	bcftools filter -e ID=@${infile_feature}_virtualPos_readSupportGT_extremGQdel_ID_list.txt -o ${infile_feature}_virtualPos_readSupportGT.vcf
fi

bcftools sort -O z -o ${infile_feature}_readSupportGT.vcf.gz ${infile_feature}_readSupportGT.vcf && \
bcftools index ${infile_feature}_readSupportGT.vcf.gz && \
rm ${infile_feature}_readSupportGT.vcf && \
bcftools sort -O z -o ${infile_feature}_virtualPos_readSupportGT.vcf.gz ${infile_feature}_virtualPos_readSupportGT.vcf && \
rm ${infile_feature}_virtualPos_readSupportGT.vcf && \
bcftools index ${infile_feature}_virtualPos_readSupportGT.vcf.gz && \
bcftools concat -a -O u -o ${infile_feature}_withVirtual_readSupportGT.vcf ${infile_feature}_readSupportGT.vcf.gz ${infile_feature}_virtualPos_readSupportGT.vcf.gz && \
rm ${infile_feature}*_readSupportGT_extremGQdel_ID_list.txt
unset infile_feature

