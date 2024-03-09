import pandas as pd
import numpy as np
import sys,re
confd=sys.argv[1]
sys.path.append(confd)
from tpc import *


data_dir=postEC_wd


df_position_preflt_cluster_info=pd.read_pickle(data_dir+'/'+cohort+'_ERVcaller_post_part2_preflt_clustermerged_cluster_info.pkl')


#make df_position_fltSetO1_cluster_info step1
df_siteMissGT=pd.read_csv(data_dir+'/'+cohort+\
                          '_ERVcaller_fltSetO1_step1_represenID_clustermergedReGT_missing.lmiss',
                         sep='\t')#.sort_values('CHR').reset_index(drop=True)
df_hwe=pd.read_csv(data_dir+'/'+cohort+\
                          '_ERVcaller_fltSetO1_step1_represenID_clustermergedReGT_hwe.hwe',
                         sep='\t')
df_freq=pd.read_csv(data_dir+'/'+cohort+\
                          '_ERVcaller_fltSetO1_step1_represenID_clustermergedReGT_frq.tsv',
                         sep='\t')
df_site_oriOrder=pd.read_csv(
    data_dir+'/'+cohort+'_ERVcaller_fltSetO1_step1_represenID_clustermergedReGT_IDoriOrder.txt',sep='\t',header=None)
    
#make df_position_fltSetO1_cluster_info step2
df_siteMissGT['pos_id']=df_site_oriOrder[0]
df_hwe['pos_id']=df_site_oriOrder[0]
df_freq['pos_id']=df_site_oriOrder[0]


df_position_fltSetO1_step1_cluster_info=df_position_preflt_cluster_info.merge(df_siteMissGT[['pos_id','F_MISS']],on='pos_id',how='inner')
df_position_fltSetO1_step1_cluster_info=df_position_fltSetO1_step1_cluster_info.merge(df_hwe[['pos_id','P_HET_EXCESS']],on='pos_id',how='inner')
df_position_fltSetO1_step1_cluster_info=df_position_fltSetO1_step1_cluster_info.merge(df_freq[['pos_id','alt_FREQ']],on='pos_id',how='inner')


#filter positions. X and Y no filtering, leave them for future.
df_position_fltSetO1_cluster_info=\
df_position_fltSetO1_step1_cluster_info.loc[
    (df_position_fltSetO1_step1_cluster_info['alt_FREQ']>0) & 
    (((df_position_fltSetO1_step1_cluster_info['F_MISS']<=0.5) &
     (df_position_fltSetO1_step1_cluster_info['P_HET_EXCESS']>0.05)) | 
    (df_position_fltSetO1_step1_cluster_info['chr'].isin(['chrX','chrY'])))].reset_index(drop=True)

#make a list of represenID_toDelete from prefilt set to fltSetO1
ser_represenID_toDelete=df_position_preflt_cluster_info['pos_id'].loc[~df_position_preflt_cluster_info['pos_id'].isin(df_position_fltSetO1_cluster_info['pos_id'])]
len(ser_represenID_toDelete)

#for sh script fltSetO1 fliltering vcf and more
ser_represenID_toDelete.to_csv(
    data_dir+'/'+cohort+"_ERVcaller_fltSetO1_represenID_toDelete_list.txt",header=None,index=None,sep='\t')

#export a .bed and full info txt of clusters of filtered Set O1 final
df_position_fltSetO1_cluster_info.to_pickle(
    data_dir+'/'+cohort+"_ERVcaller_post_part2_fltSetO1_clustermerged_cluster_info.pkl")

df_position_fltSetO1_cluster_info[['chr','cluster_bed_bg','cluster_bed_ed']].to_csv(
    data_dir+'/'+cohort+"_ERVcaller_fltSetO1_clustermerged_cluster_pos.bed",header=None,index=None,sep='\t')

#df_position_fltSetO1_cluster_info.drop(columns=['cluster']).to_csv(
#    data_dir+'/'+cohort+"_ERVcaller_fltSetO1_clustermerged_cluster_info.txt",header=1,index=1,sep='\t')

df_position_fltSetO1_cluster_info[['pos_id','alt_FREQ']].to_csv(
    data_dir+'/'+cohort+"_ERVcaller_fltSetO1_freq.txt",header=1,index=None,sep='\t')
