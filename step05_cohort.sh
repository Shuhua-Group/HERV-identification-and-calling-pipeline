#!/usr/bin/bash
#cmd: bash step05_cohort.sh /PATH/TO/YOUR.config
#echo "The script is under revision. Please wait for the release in near future.";exit


if [ $1 = "" ];then
	echo "Usage: bash step05_cohort.sh /PATH/TO/YOUR.config";exit
fi


configPath=$1
check_config=$(awk -v configPath=$configPath '{sub(/#[^\n]+/,"",$0);$1=$1;print $0}' $configPath | egrep "\s")
if [ "$check_config" != "" ];then
	echo "Failed. Please correct the content in your config file.";exit
fi
. $configPath


bash ${postERVcaller}/scripts/s1m2.sh $configPath || echo "Step05 failed."
