#!/usr/bin/bash
#cmd: bash s2m2sap4.sh PATH/TO/YOUR.s.config
#echo "The script is under revision. Please wait for the release in near future.";exit


configPath=$1
. $configPath


cd $postEC_wd

awk 'NR>1{sub(":","\t",$2);print $2}' ${cohort}_siteAnal_panERV_finflt_clustermerged_cluster_info.txt > ${cohort}_siteAnal_panERV_finflt_poslist.txt
bcftools view -o ${cohort}_ERVcaller_fltSetO1_represenID_clustermergedReGT.bcf ${cohort}_ERVcaller_fltSetO1_represenID_clustermergedReGT.vcf
bcftools index ${cohort}_ERVcaller_fltSetO1_represenID_clustermergedReGT.bcf
bcftools view -R ${cohort}_siteAnal_panERV_finflt_poslist.txt -o ${cohort}_siteAnal_panERV_finflt.vcf ${cohort}_ERVcaller_fltSetO1_represenID_clustermergedReGT.bcf

rm ${cohort}_ERVcaller_fltSetO1_represenID_clustermergedReGT.bcf*

