#!/usr/bin/bash
#cmd: bash step07_cohort.sh /PATH/TO/YOUR.config
#echo "The script is under revision. Please wait for the release in near future.";exit


if [[ $1 = "" ]];then
	echo "Usage: bash step07_cohort.sh /PATH/TO/YOUR.config";exit
fi


configPath=$1
check_config=$(awk -v configPath=$configPath '{sub(/#[^\n]+/,"",$0);$1=$1;print $0}' $configPath | egrep "\s")
if [ "$check_config" != "" ];then
	echo "Failed. Please correct the content in your config file.";exit
fi
. $configPath

failed=0
confb=$(basename $configPath);confd=${configPath/"/${confb}"/};cp $configPath ${confd}/tpc.py && \
bash ${postERVcaller}/scripts/s1m3s.sh $configPath && \
$python3 ${postERVcaller}/p/ppart3.py $confd && \
bash ${postERVcaller}/scripts/s2m1sap3.sh $configPath && \
#bedtools intersect -loj -a ${cohort}_siteAnal_panERV_good1_cluster_range.bed -b ${postERVcaller}/p/HERVdEnt.bed > ${cohort}_siteAnal_panERV_good1_inERV.bed
$python3 ${postERVcaller}/p/ppart4.py $confd && \
bash ${postERVcaller}/scripts/s2m2sap4.sh $configPath || failed=1
rm ${confd}/tpc.py
unset confd
if [ ! $failed -eq 0 ]; then
	echo "Step07 failed."
else
	#cleaning
	mkdir -p backup
	mv ${cohort}_*vcf* backup/
	mv ${cohort}_*pkl backup/
	mv ${cohort}_*list* backup/
	mv ${cohort}_*merge* backup/
	mv ${cohort}_ERVcaller_*flt* backup/
	mv ${cohort}_siteAnal_* backup/
	for i in ${cohort}_siteAnal_panERV_finflt_clustermerged_cluster_info.txt ${cohort}_siteAnal_panERV_finflt.vcf ${cohort}_siteAnal_panERV_finflt_pos_geno.pkl;do
		mv -f backup/$i ./
	done
fi
