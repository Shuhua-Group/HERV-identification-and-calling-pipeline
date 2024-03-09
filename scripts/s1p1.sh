#!/usr/bin/bash
#cmd: bash s1p1.sh /PATH/TO/YOUR.config sampleID


configPath=$1
. $configPath
i=$2
if [[ $i = "" ]]; then echo "Errors in arguments."; exit;fi


inBamDir=$EC_out_dir/h_bam
inBamMetricsDir=$EC_out_dir/size
if [ -f $EC_out_dir/rlen/${i}.rl ]; then readLen=$(cat $EC_out_dir/rlen/${i}.rl); fi
read -ra insertMetrics <<< $(awk 'FNR==1{sub("\\.[0-9]+","",$6); sub("\\.[0-9]+","",$8);print $6" "(($8 + 1))}' ${inBamMetricsDir}/${i}_h.bamInsertSizeMetrics)
insertLenAvg=${insertMetrics[0]}
insertLenStd=${insertMetrics[1]}

#reads-supported non-insertion typing for preflt set
cd $postEC_wd/reads_support_noIns
perl ${ERVcaller}/Scripts/Calculate_reads_among_nonTE_locations.pl -i ../${cohort}_ERVcaller_preflt_oriPosID.vcf -S $i -o $i.nonTE -b ${inBamDir}/${i}_h.bam -s paired-end -l $insertLenAvg -L $insertLenStd -r $readLen -t 16

#reads-supported non-insertion typing for virt set
cd $postEC_wd/reads_support_noIns_virt
perl ${ERVcaller}/Scripts/Calculate_reads_among_nonTE_locations.pl -i ../${cohort}_ERVcaller_preflt_virtualPos_for_readSupportInfo.vcf -S $i -o $i.nonTE -b ${inBamDir}/${i}_h.bam -s paired-end -l $insertLenAvg -L $insertLenStd -r $readLen -t 16
