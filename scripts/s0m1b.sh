#!/usr/bin/bash
#cmd: bash sm01.sh /PATH/TO/YOUR.config


configPath=$1
. $configPath

to_extit=0
#prepare files and folders for next steps
cd $EC_out_dir
mkdir -p size
for i in $EC_out_dir/job_io/*\.out; do
	line=$(grep -a "# sample_ID:" $i)
	if [[ $line != "" ]];then
		id=$(echo $line|egrep -o "sample_ID: \S+");id=${id/sample_ID: /};echo $line > $EC_out_dir/size/${id}_h.bamInsertSizeMetrics
	else
		echo "Content of $i is incomplete."; to_extit=1
	fi
done
if [ "$to_extit" -eq 1 ];then
	echo "Errors occur.";exit
fi
mkdir -p adpt_vcf
#mkdir -p rlen
