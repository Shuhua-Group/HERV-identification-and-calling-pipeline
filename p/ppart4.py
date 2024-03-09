import pandas as pd
import numpy as np
import os,sys,re
confd=sys.argv[1]
sys.path.append(confd)
from tpc import *


data_dir=postEC_wd


df_fltSetO1_poscluster=pd.read_pickle(data_dir+'/'+cohort+'_ERVcaller_post_part2_fltSetO1_clustermerged_cluster_info.pkl')
df_classification=pd.read_csv(postERVcaller+'/'+'p/element.info',sep='\t')


df_siteAnal_ERV=df_fltSetO1_poscluster.loc[
    (df_fltSetO1_poscluster['both_ends_evident']==True),:]
df_siteAnal_ERV_good1=df_siteAnal_ERV.loc[
    (df_siteAnal_ERV['cluster_len']<200)|
    ((df_siteAnal_ERV['cluster_len']>=200)&
     (df_siteAnal_ERV['cluster_len']<400)&
     (df_siteAnal_ERV['status_max']>=4))]

df_siteAnal_partialERV=df_fltSetO1_poscluster.loc[
    (df_fltSetO1_poscluster['ERVfalseSig1_oneEnd_only']==True),:]
df_siteAnal_partialERV_good1=df_siteAnal_partialERV.loc[
    (df_fltSetO1_poscluster['cluster_len']<200),:]

df_siteAnal_panERV_good1=df_fltSetO1_poscluster.loc[
    (df_fltSetO1_poscluster['pos_id'].isin(df_siteAnal_ERV_good1['pos_id']))|
    (df_fltSetO1_poscluster['pos_id'].isin(df_siteAnal_partialERV_good1['pos_id'])),:].copy()

df_siteAnal_panERV_good1['panERV_type']=\
df_siteAnal_panERV_good1['ERVfalseSig1_oneEnd_only'].apply(lambda x: 'partial-ERV' if x else 'ERV')


df_siteAnal_panERV_good1[['chr','cluster_bed_bg','cluster_bed_ed','element_maj','pos_id','element_set']].to_csv(
    data_dir+'/'+cohort+
    '_siteAnal_panERV_good1_cluster_range.bed',index=None,header=None,sep='\t')
os.system("bedtools intersect -loj -a "+cohort+"_siteAnal_panERV_good1_cluster_range.bed -b "+postERVcaller+"/p/HERVdEnt.bed > "+cohort+"_siteAnal_panERV_good1_inERV.bed")

df_inERV=pd.read_csv(
    data_dir+'/'+cohort+'_siteAnal_panERV_good1_inERV.bed',sep='\t',header=None).\
iloc[:,:10]
df_inERV.columns=['chr','cluster_bed_bg','cluster_bed_ed','element_maj','pos_id','element_set',
                  'in_chr','in_bed_bg','in_bed_ed',
                  'in_entity']

def get_in_ERV_len(rowser_in):
    if rowser_in['in_bed_bg']==-1:
        return None
    else:
        avgPos=(rowser_in['cluster_bed_bg']+rowser_in['cluster_bed_ed'])/2
        avgBackgroundPos=(rowser_in['in_bed_bg']+rowser_in['in_bed_ed'])/2
        cluster_len=rowser_in['cluster_bed_ed']-rowser_in['cluster_bed_bg']
        background_len=rowser_in['in_bed_ed']-rowser_in['in_bed_bg']
        bias=avgPos-avgBackgroundPos
        bias_ratio=abs(bias/((cluster_len+background_len)/2))
        return bias_ratio
        
df_inERV['inERV_pos_bias_ratio']=\
df_inERV[['cluster_bed_bg','cluster_bed_ed','in_bed_bg','in_bed_ed']].apply(lambda x: get_in_ERV_len(x), axis=1)

#df_inERV_keyinfo=\
#df_inERV[['pos_id','in_chr','in_entity','in_bed_bg','in_bed_ed','inERV_pos_bias_ratio']]\
#.reset_index()\
#.sort_values(by='inERV_pos_bias_ratio')\
#.groupby('pos_id').\
#agg({'index':'first','in_chr':'count','in_entity':list,
#     'in_bed_bg':list,'in_bed_ed':list,
#    'inERV_pos_bias_ratio':list})\
#.sort_values(by='index').reset_index().set_index('index')\
#.rename(columns={'in_chr':'in_entity_count'})

df_inERV_used=\
df_inERV[['pos_id','in_chr','in_entity','in_bed_bg','in_bed_ed','inERV_pos_bias_ratio']]\
.reset_index()\
.sort_values(by='inERV_pos_bias_ratio')\
.groupby('pos_id').\
agg({'index':'first','in_chr':'count','in_entity':'first',
     'in_bed_bg':'first','in_bed_ed':'first',
    'inERV_pos_bias_ratio':'first'})\
.sort_values(by='index').reset_index().set_index('index')\
.rename(columns={'in_chr':'in_entity_count'})


df_finflt=df_siteAnal_panERV_good1.copy()
df_finflt=df_finflt.merge(df_inERV_used[['pos_id','in_entity']],
                 on='pos_id',how='left')
df_finflt['in_entity']=df_finflt['in_entity'].apply(lambda x: x!='.')
df_finflt=df_finflt.rename(columns={'in_entity':'in_ERV'})

def get_panERV_type2(col_ser_in):
    if col_ser_in['panERV_type'] =='ERV':
        if col_ser_in['in_ERV']==False:
            type2='free ERV'
        else:
            type2='nested ERV'
    else:
        if col_ser_in['in_ERV']==False:
            type2='free partial-ERV'
        else:
            type2='nested partial-ERV'
    return type2
        
df_finflt['panERV_type2']=df_finflt[['panERV_type','in_ERV']].apply(lambda x: get_panERV_type2(x),axis=1)


df_finflt.drop(columns=['status_dedup', 'status_max', 'in_ERV', 'panERV_type', 'ERVfalseSig1_oneEnd_only', 'ERVfalseSig2_otherRpt_largely', 'ERVfalseSig3_simpleRpt_largely', 'left_status_max', 'right_status_max','notOtherRepeat_score_max','notSimpRpt_score_max','alt_FREQ_max','alt_FREQ_sum'], inplace=True)

#Export
df_finflt.to_csv(
    data_dir+'/'+cohort+'_siteAnal_panERV_finflt_clustermerged_cluster_info.txt',
    header=1,index=1,sep='\t')
print('Got',len(df_finflt),'loci after final filtering.',sep=' ')



#Export
#df_siteAnal_panERV_good1.to_csv(
#    data_dir_o+popCodeName+'_siteAnal_panERV_good1_clustermerged_cluster_info.txt',
#    header=1,index=1,sep='\t')
#print('Got',len(df_siteAnal_panERV_good1),'loci after final filtering.',sep=' ')