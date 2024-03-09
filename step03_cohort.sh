#!/usr/bin/bash
#cmd: bash step03.sh /PATH/TO/YOUR.config
#echo "The script is under revision. Please wait for the release in near future.";exit


if [ $1 = "" ];then
	echo "Usage: bash step03.sh /PATH/TO/YOUR.config";exit
fi


configPath=$1
check_config=$(awk -v configPath=$configPath '{sub(/#[^\n]+/,"",$0);$1=$1;print $0}' $configPath | egrep "\s")
if [ "$check_config" != "" ];then
	echo "Failed. Please correct the content in your config file.";exit
fi
. $configPath

failed=0
confb=$(basename $configPath);confd=${configPath/"/${confb}"/};cp $configPath ${confd}/tpc.py && \
bash ${postERVcaller}/scripts/s0m2bp1.sh $configPath && \
$python3 ${postERVcaller}/p/ppart1.py $confd && \
$python3 ${postERVcaller}/p/ppart2.py $confd && \
bash ${postERVcaller}/scripts/s1m1sap2.sh $configPath || failed=1
rm ${confd}/tpc.py
unset confd
rm -f ${postEC_wd}/${cohort}_ERVcaller_post_part1_indivmerge_site_info.pkl
if [ ! $failed -eq 0 ]; then
	echo "Step03 failed."
fi
