#!/usr/bin/bash
#cmd: bash s1p2.sh /PATH/TO/YOUR.config sampleID


configPath=$1
. $configPath
i=$2
if [[ $i = "" ]]; then echo "Errors in arguments."; exit;fi


cd $postEC_wd

infile_feature=${cohort}_ERVcaller_preflt
in_dir=$postEC_wd
$python3 ${postERVcaller}/s/represenPOS_GTreplace.py \
$i \
${in_dir}/${infile_feature}_withVirtual_readSupportGT.vcf \
${in_dir}/${infile_feature}_pos_cluster_info.txt \
${in_dir}/${infile_feature}_clustermerged_cluster_info.txt


