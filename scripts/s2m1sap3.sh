#!/usr/bin/bash
#cmd: bash s2m1sap3.sh /PATH/TO/YOUR.s.config
#echo "The script is under revision. Please wait for the release in near future.";exit


configPath=$1
. $configPath


cd $postEC_wd

if [ ! -f ${cohort}_ERVcaller_fltSetO1_step1_pos_cluster_info.txt ]; then
	ln -s ${cohort}_ERVcaller_preflt_pos_cluster_info.txt ${cohort}_ERVcaller_fltSetO1_step1_pos_cluster_info.txt
fi
if [ ! -f ${cohort}_ERVcaller_fltSetO1_step1_clustermerged_cluster_info.txt ]; then
	ln -s ${cohort}_ERVcaller_preflt_clustermerged_cluster_info.txt ${cohort}_ERVcaller_fltSetO1_step1_clustermerged_cluster_info.txt
fi
bash ${postERVcaller}/s/s3.sh $configPath


#END
