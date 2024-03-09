#!/usr/bin/bash
#cmd: bash sm02.sh PATH/TO/YOUR.s.config


configPath=$1
. $configPath


cd $postEC_wd

#prepare .bcf list file
IDprefix=""
IDsufffix=_adpt
ls $EC_out_dir/${IDprefix}*${IDsufffix}_sorted\.bcf > ${cohort}_ERVcaller_bcf_list.txt

#merge
bcftools merge --threads -0 --info-rules INFOR:join -l ${cohort}_ERVcaller_bcf_list.txt -O v -o ${cohort}_ERVcaller_merge.vcf

#get positions and alleles info
#	_ERVcaller_merge_siteinfo.txt
pivcf=${cohort}_ERVcaller_merge.vcf
awk 'BEGIN{print "chr\tpos\tdirection\tref_allele\talt_allele\tmax_aln_len\taln_bg\taln_ed\tmax_status\tmax_left_status\tmax_right_status\tTE_acc"} {if ($0 !~ /^#/) {split($5,arr1,":");split(arr1[2],arr2,">");split($8,arr3,";");split(arr3[6],arr4,",");split(arr4[6],arr5,"_");\
sub("INFOR=","",arr4[1]);alnbg=alnbg2=arr4[2];alned=alned2=arr4[3];maxAlnLen=arr4[4];direction=arr4[5];maxstatus=arr5[1];lmaxstatus=arr5[2];rmaxstatus=arr5[3];\
if(length(arr4)>6){for(i=7;i<=length(arr4);i++){if((i%6)==4){if(arr4[i]!="NULL" && (maxAlnLen=="NULL" || arr4[i]>maxAlnLen)) {maxAlnLen=arr4[i];alnbg=arr4[i-2];alned=arr4[i-1]} \
else if(arr4[i]=="NULL" && arr4[i-2]!="NULL" && (alnbg2=="NULL" || arr4[i-2]<alnbg2)){alnbg2=arr4[i-2]} \
else if(arr4[i]=="NULL" && arr4[i-1]!="NULL" && (alned2=="NULL" || arr4[i-1]>alned2)){alned2=arr4[i-1]}};\
if((i%6)==5 && arr4[i]!=arr4[i-6]){direction="NULL"};if((i%6)==0) {if(arr5[1]>maxstatus){maxstatus=arr5[1]};\
if(arr5[2]>lmaxstatus){lmaxstatus=arr5[2]};\
if(arr5[3]>rmaxstatus){rmaxstatus=arr5[3]};\
}}};\
if(maxAlnLen=="NULL"){alnbg=alnbg2;alned=alned2};\
print $1"\t"$2"\t"direction"\t"$4"\t"arr2[1]"\t"maxAlnLen"\t"alnbg"\t"alned"\t"maxstatus"\t"lmaxstatus"\t"rmaxstatus"\t"arr4[1]}}' $pivcf > ${pivcf/.vcf/_siteinfo.txt}

#get _addID.vcf and freq
bash ${postERVcaller}/s/s1.sh $configPath

mkdir -p reads_support_noIns
mkdir -p reads_support_noIns_virt
