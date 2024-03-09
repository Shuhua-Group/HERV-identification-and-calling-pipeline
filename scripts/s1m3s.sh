#!/usr/bin/bash
#cmd: bash s1m3s.sh /PATH/TO/YOUR.s.config
#echo "The script is under revision. Please wait for the release in near future.";exit


configPath=$1
. $configPath


cd $postEC_wd

bash ${postERVcaller}/s/s2.sh $configPath

