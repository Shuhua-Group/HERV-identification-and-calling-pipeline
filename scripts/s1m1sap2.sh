#!/usr/bin/bash
#cmd: bash s1m1sap2.sh /PATH/TO/YOUR.s.config


configPath=$1
. $configPath


cd $postEC_wd

#get preflt_oriPosID.vcf
#---prepare for reads-supported_nonInsertion procedure bg
#files needed in advance: reference genome, dictionary of the reference genome, the 2 files of header lines and variant by spliting the original xxx_addID.vcf.
vcftools --snps ${cohort}_ERVcaller_preflt_posID_list.txt --vcf ${cohort}_ERVcaller_addID.vcf --recode --recode-INFO-all --out ${cohort}_ERVcaller_preflt_oriPosID
mv ${cohort}_ERVcaller_preflt_oriPosID.recode.vcf ${cohort}_ERVcaller_preflt_oriPosID.vcf
#---prepare for reads-supported_nonInsertion procedure ed

#make vcf of the virt positions
infile_feature=${cohort}_ERVcaller_preflt
sample_num=$(grep "#CHROM" ${infile_feature}_oriPosID.vcf|awk '{split($0,arr,"\t"); for (elements in arr) sn++; sn-=9; print sn}')
awk -v sample_num=$sample_num \
'NR>1{arr[1]=300;arr[2]=600;arr[3]=900;arr[4]=1200;arr[5]=1500; \
for (i=1; i<=5; i++) {up[i]=$7-arr[i];down[i]=$8+arr[i]; GTs=""; for (j=1;j<=sample_num;j++){GTs=GTs"\t0/0:.:.:.:."}; sub("\t","",GTs); print $4"\t"up[i]"\t"$3"\t.\tINS\t.\t.\t.\tGT:GQ:GL:DPN:DPI\t"GTs; print $4"\t"down[i]"\t"$3"\t.\tINS\t.\t.\t.\tGT:GQ:GL:DPN:DPI\t"GTs}}' \
${infile_feature}_clustermerged_cluster_info.txt > ${infile_feature}_virtualPos_for_readSupportInfo.vcf_records
check_fsize=$(ls -l ${infile_feature}_virtualPos_for_readSupportInfo.vcf_records | egrep -o "\s[0-9]{3,}\s")
if [ check_fsize = "" ];then
	echo "Errors occurred and exited.";exit
fi
bcftools view -h ${infile_feature}_oriPosID.vcf > ${infile_feature}_virtualPos_for_readSupportInfo.vcf && \
cat ${infile_feature}_virtualPos_for_readSupportInfo.vcf_records >> ${infile_feature}_virtualPos_for_readSupportInfo.vcf && \
rm ${infile_feature}_virtualPos_for_readSupportInfo.vcf_records && unset infile_feature

