import pandas as pd
import numpy as np
import sys,re
confd=sys.argv[1]
sys.path.append(confd)
from tpc import *


data_dir=postEC_wd
df_position=pd.read_pickle(data_dir+'/'+cohort+'_ERVcaller_post_part1_indivmerge_site_info.pkl')

#Make clustering ground information for all positions.
ser_last_chr_num=pd.Series(np.array([0]+df_position['chr_num'].iloc[:-1].tolist()))
ser_last_pos=pd.Series(np.array([0]+df_position['pos'].iloc[:-1].tolist()))

df_posDistance=pd.concat([(df_position['pos']-ser_last_pos),(df_position['chr_num']-ser_last_chr_num==0)],axis=1)\
                .apply(lambda x: x[0] if x[1]>0 else -1, axis=1)\
.reset_index().rename(columns={'index':'order',0:'distance_to_last'})

df_posDistance['D1000_start']=((df_position['pos']-ser_last_pos>1000) + (df_position['chr_num']-ser_last_chr_num>0))
df_posDistance['D1000cluster']=df_posDistance.apply(lambda x: df_posDistance.iloc[:x['order']+1]['D1000_start'].sum(),axis=1)

df_posDistance['D100_start']=((df_position['pos']-ser_last_pos>100) + (df_position['chr_num']-ser_last_chr_num>0))
df_posDistance['D100cluster']=df_posDistance.apply(lambda x: df_posDistance.iloc[:x['order']+1]['D100_start'].sum(),axis=1)

df_posDistance['D50_start']=((df_position['pos']-ser_last_pos>50) + (df_position['chr_num']-ser_last_chr_num>0))
df_posDistance['D50cluster']=df_posDistance.apply(lambda x: df_posDistance.iloc[:x['order']+1]['D50_start'].sum(),axis=1)

df_posDistance['D20_start']=((df_position['pos']-ser_last_pos>20) + (df_position['chr_num']-ser_last_chr_num>0))
df_posDistance['D20cluster']=df_posDistance.apply(lambda x: df_posDistance.iloc[:x['order']+1]['D20_start'].sum(),axis=1)

ser_last_312class=pd.Series(['non']+df_position['312_class'].iloc[:-1].tolist())
df_posDistance['312class_start']=((df_position['312_class']!=ser_last_312class) + (df_position['chr_num']-ser_last_chr_num>0))
df_posDistance['312class_cluster']=df_posDistance.apply(lambda x: df_posDistance.iloc[:x['order']+1]['312class_start'].sum(),axis=1)

ser_last_TE_acc=pd.Series(['non']+df_position['TE_acc'].iloc[:-1].tolist())
df_posDistance['TE_acc_start']=((df_position['TE_acc']!=ser_last_TE_acc) + (df_position['chr_num']-ser_last_chr_num>0))
df_posDistance['TE_acc_cluster']=df_posDistance.apply(lambda x: df_posDistance.iloc[:x['order']+1]['TE_acc_start'].sum(),axis=1)

ser_last_TE_acc=pd.Series(['non']+df_position['direction'].iloc[:-1].tolist())
df_posDistance['direction_start']=((df_position['direction']!=ser_last_TE_acc) + (df_position['chr_num']-ser_last_chr_num>0))
df_posDistance['direction_cluster']=df_posDistance.apply(lambda x: df_posDistance.iloc[:x['order']+1]['direction_start'].sum(),axis=1)

df_pos_posDistance=df_position[['pos_id','chr_num','pos','max_status']].join(df_posDistance) #,'alt_FREQ','P_HET_EXCESS'
df_pos_posDistance.drop(columns=['order'],inplace=True)
df_pos_posDistance.index.rename('order',inplace=True)
df_pos_posDistance.reset_index(inplace=True)

#add cluster info. time consuming (>3 hours).
df_pos_posDistance['cluster_bg']=\
df_pos_posDistance.apply(lambda x: (x['D1000_start']) | 
(x['direction_start']) |
((x['TE_acc_start']) & (x['D100_start'])) |
(x['312class_start']),axis=1)
df_pos_posDistance['cluster']=df_pos_posDistance.apply(lambda x: df_pos_posDistance.iloc[:x['order']+1]['cluster_bg'].sum(),axis=1)


#export
df_pos_posDistance.to_pickle(data_dir+'/'+cohort+'_pos_clustering_ground_info.pkl')


#add clustering bg-ed info
df_cluster0=pd.DataFrame(df_pos_posDistance['cluster'].\
                         groupby(df_pos_posDistance['cluster']).first())
df_cluster0['cluster_bg']=df_pos_posDistance[['cluster','pos']].groupby('cluster').min()
df_cluster0['cluster_ed']=df_pos_posDistance[['cluster','pos']].groupby('cluster').max()
df_cluster0['cluster_region_len']=df_cluster0['cluster_ed']-df_cluster0['cluster_bg']+1
df_cluster0=df_cluster0.reset_index(drop=True)

if 'cluster' not in df_position.columns:
    df_position=df_position.join(df_pos_posDistance['cluster'])
    df_position=df_position.merge(df_cluster0,on='cluster',how='left')


df_position_view=df_position[['cluster','chrpos','pos_id','chr','chr_num',
                              'pos','cluster_bg','cluster_ed',
                              'element','organism','clades','alt_FREQ',
                              'max_status','max_left_status','max_right_status',
                              'max_aln_len', 'notOtherRepeat_score', 'notSimpRpt_score',
                              '312_class','ERV_part']].copy()
df_position_view['clades']=df_position_view['clades'].fillna('noInfo')

df_position_view['cluster_element']=df_position_view[['cluster','element']].apply(lambda x:str(x[0])+'_'+x[1],axis=1)

df_temp=\
df_position_view['cluster_element'].groupby(df_position_view['cluster_element'])\
.count().rename('cluster_element_count').reset_index()

df_position_view=df_position_view.merge(df_temp,on='cluster_element',how='left')


#do the merge
df_position_view['element_set']=df_position_view['element']
df_position_view['element_maj']=df_position_view['element']
df_position_view_cluster=df_position_view.sort_values(
    by=['chr_num','cluster','max_status','cluster_element_count','alt_FREQ','max_aln_len','pos'],
    ascending=[True,True,False,False,False,False,True]).groupby(
    'cluster').agg({
    'chrpos':'first','pos_id':'first','chr':'first','chr_num':'first','pos':'first',
    'cluster_bg':'first','cluster_ed':'first',
    'element_set':lambda x: '|'.join(sorted(list(set(x.to_list())))),
    'element_maj':'first',
    'organism':lambda x: '|'.join(sorted(list(set(x.to_list())))),
    'clades':lambda x: '|'.join(sorted(list(set(x.to_list())))),
    'max_status':lambda x: x.to_list(),
    'max_left_status':lambda x: x.to_list(),
    'max_right_status':lambda x: x.to_list(),
    'alt_FREQ':lambda x: x.to_list(),
    'notOtherRepeat_score':lambda x: x.to_list(),
    'notSimpRpt_score':lambda x: x.to_list(),
    'max_aln_len':lambda x: x.to_list(),
    '312_class':lambda x: '|'.join(sorted(list(set(x.to_list())))),
    'ERV_part':lambda x: '|'.join(sorted(list(set(x.to_list()))))
    }).sort_values(by=['chr_num','cluster'])

#make a table for export
df_position_cluster_info=df_position_view_cluster.copy()
df_position_cluster_info['alt_FREQ_max']=df_position_cluster_info['alt_FREQ'].apply(lambda x: max(x))
df_position_cluster_info['alt_FREQ_sum']=df_position_cluster_info['alt_FREQ'].apply(lambda x: sum(x))
df_position_cluster_info['status_dedup']=df_position_cluster_info['max_status'].apply(
    lambda x: str(list(set(x))).strip('[]').replace(' ',''))
df_position_cluster_info['status_max']=df_position_cluster_info['max_status'].apply(lambda x: max(x))
df_position_cluster_info['left_status_max']=df_position_cluster_info['max_left_status'].apply(lambda x: max(x))
df_position_cluster_info['right_status_max']=df_position_cluster_info['max_right_status'].apply(lambda x: max(x))
df_position_cluster_info['notOtherRepeat_score_max']=df_position_cluster_info['notOtherRepeat_score'].apply(lambda x: max(x))
df_position_cluster_info['notSimpRpt_score_max']=df_position_cluster_info['notSimpRpt_score'].apply(lambda x: max(x))
df_position_cluster_info['aln_len_max']=df_position_cluster_info['max_aln_len'].apply(lambda x: max(x))

df_position_cluster_info=df_position_cluster_info\
.drop(columns=['alt_FREQ','max_status','max_left_status', 'max_right_status',
               'notOtherRepeat_score','notSimpRpt_score'])


#add more info columns
df_position_cluster_info['both_ends_evident']=\
df_position_cluster_info[['left_status_max', 'right_status_max']].apply(lambda x: x[0]>0 and x[1]>0,axis=1)
df_position_cluster_info['ERVfalseSig1_oneEnd_only']=\
df_position_cluster_info[['left_status_max', 'right_status_max']].apply(lambda x: int(x[0]>0) + int(x[1]>0) == 1 ,axis=1)
df_position_cluster_info['ERVfalseSig2_otherRpt_largely']=\
df_position_cluster_info['notOtherRepeat_score_max'].apply(lambda x: x<20)
df_position_cluster_info['ERVfalseSig3_simpleRpt_largely']=\
df_position_cluster_info[['notSimpRpt_score_max']].apply(lambda x: x<20)


#prefiltering
df_position_preflt_cluster_info=df_position_cluster_info.loc[
    (df_position_cluster_info['ERVfalseSig3_simpleRpt_largely']==False) &\
    (df_position_cluster_info['status_max']>=1)
].copy()

#prepare: make position .bed and full info txt of cluster prefilt Set.
df_position_preflt_cluster_info['cluster_len']=df_position_preflt_cluster_info['cluster_ed']-df_position_preflt_cluster_info['cluster_bg']
df_position_preflt_cluster_info['cluster_bed_bg']=df_position_preflt_cluster_info['cluster_bg']-1
df_position_preflt_cluster_info['cluster_bed_ed']=df_position_preflt_cluster_info['cluster_ed']-1


df_position_preflt=df_position.loc[df_position['cluster'].isin(df_position_preflt_cluster_info.index)]
df_position_preflt=df_position_preflt.reset_index(drop=True)


#export prefiltered positions (not merged)
df_position_preflt.to_pickle(data_dir+'/'+cohort+'_ERVcaller_post_part2_preflt_pos_cluster_info.pkl')

df_position_preflt[['chr','pos','pos_id','cluster']].to_csv(
    cohort+"_ERVcaller_preflt_pos_cluster_info.txt",header=1,index=1,sep='\t')

#export prefiltered position merged clusters info (merged)
#print(len(df_position_preflt_cluster_info))
df_position_preflt_cluster_info.to_pickle(data_dir+'/'+cohort+'_ERVcaller_post_part2_preflt_clustermerged_cluster_info.pkl')

df_position_preflt_cluster_info.to_csv(data_dir+'/'+cohort+'_ERVcaller_preflt_clustermerged_cluster_info.txt',
                                          header=1,index=1,sep='\t')

#export prefilt set positions.
#chr' - 'pos not always unique. pos_id are unique.
df_position_preflt[['pos_id']].to_csv(
    data_dir+'/'+cohort+"_ERVcaller_preflt_posID_list.txt",
    header=0,index=None)
print('Got',len(df_position_preflt),'loci after prefiltering.',sep=' ')