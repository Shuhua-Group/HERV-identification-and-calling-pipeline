#!/usr/bin/bash
#cmd: bash step02_indiv.sh /PATH/TO/YOUR.config sampleID
#echo "The script is under revision. Please wait for the release in near future.";exit


if [ $2 = "" ];then
	echo "Usage: bash step02_indiv.sh /PATH/TO/YOUR.config sampleID";exit
fi


configPath=$1
check_config=$(awk -v configPath=$configPath '{sub(/#[^\n]+/,"",$0);$1=$1;print $0}' $configPath | egrep "\s")
if [ "$check_config" != "" ];then
	echo "Please correct the content in your config file.";exit
fi
. $configPath


bash ${postERVcaller}/scripts/s0p1.sh $configPath $2 || echo "Step02 failed."
