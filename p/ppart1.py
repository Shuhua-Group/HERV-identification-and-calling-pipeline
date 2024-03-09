import pandas as pd
import numpy as np
import sys,re
confd=sys.argv[1]
sys.path.append(confd)
from tpc import *


data_dir=postEC_wd

#Make Dfam element family info table
df_family=pd.read_csv(postERVcaller+'/p/element.info',sep='\t',header=0)
df_family_internalRepeat=pd.read_csv(postERVcaller+'/p/masked.info',sep='\t',header=0)

df_family_internalRepeat['acc']=df_family_internalRepeat['query_acc']
df_family_internalRepeat['matchRange_on_element']=\
df_family_internalRepeat[['query_begin','query_end']].apply(lambda x: [x[0],x[1]], axis=1)
df_family_internalRepeat['matchLen_on_element']=df_family_internalRepeat['query_end']-df_family_internalRepeat['query_begin']

df_family_internalRepeat_grouped=\
df_family_internalRepeat.groupby('query_acc').aggregate\
({'acc':'first', 'matching':lambda x: list(x), 'matching_family':lambda x: list(x), 
  'matchRange_on_element':lambda x: list(x), 'matchLen_on_element':lambda x: list(x)})

df_family=df_family.iloc[:,[0,1,2,3,4,5,6,8]].merge(df_family_internalRepeat_grouped,how='outer',on='acc')


# Make df_freqANDhwe
df_freqANDhwe=pd.read_csv(data_dir+'/'+cohort+'_ERVcaller_addID_freq.txt',sep='\t',header=0)
df_freqANDhwe['pos_id']=df_freqANDhwe.loc[:,['CHROM','POS','alt_ALLELE']].apply(
    lambda x: str(x[0])+':'+str(x[1])+':'+str(x[2]),axis=1).astype('string')


df_freqANDhwe['pos_id']=df_freqANDhwe['pos_id'].apply(lambda x: x.replace('GI|','gi|'))
df_freqANDhwe['alt_ALLELE']=df_freqANDhwe['alt_ALLELE'].apply(lambda x: x.replace('GI|','gi|'))
df_freqANDhwe['pos_id'].apply(lambda x: x.replace('EMB|','emb|'))
df_freqANDhwe['pos_id']=df_freqANDhwe['pos_id'].apply(lambda x: x.replace('EMB|','emb|'))


#to quantify the possinility that the predicted insertion of ERVs is not a repeat sequence within the ERV element (low complex repeat or LINEs, etc.)
def quan_not_other_repeats(left,right,masked_ranges,masked_families,masked_matches):
    #---
    sig_min=10
    #---
    bg_ed_both_not_nan=(left>=0 or left<0) and (right>=0 or right<0)
    bg_ed_one_nan=(left>=0 or left<0) + (right>=0 or right<0) == 1 
    
    rg_cv_bg=rg_cv_ed=0 #1-based; virtual, can exceed actual range index.
    notOtherRepeat_score=1
    notSimpRpt_score=1
    ends_residue=0
    
    rep_unit_l=rep_unit_r=None
    same=False
    
    if type(masked_ranges)!=list:
        go_on=False
        notOtherRepeat_score=sig_min*2
        notSimpRpt_score=sig_min*2

    elif not bg_ed_one_nan and not bg_ed_both_not_nan:
        go_on=False
    else:
        go_on=True
    
    if go_on:
        if left >= 0 or left<0:
            rg_cv_bg=0
            left_residue=max(masked_ranges[0][0]-left,0)
            for e,rg in enumerate(masked_ranges):
                if left <= rg[1]:
                    rg_cv_bg=e
                    left_residue=max(rg[0]-left,0)
                    break
        else:
            rg_cv_bg=-1
            left_residue=0
        
        if right >= 0 or right<0:
            rg_cv_ed=len(masked_ranges)-1
            right_residue=max(right-masked_ranges[-1][1],0)
            for e,rg in enumerate(masked_ranges):
                if right < rg[0] and right < rg[1]:
                    if e>0:
                        rg_cv_ed=e-1
                        right_residue=max(right-masked_ranges[e-1][1],0)
                    else:
                        rg_cv_ed=0
                        right_residue=right-rg[0]+1
                    break
            if rg_cv_bg==-1:
                rg_cv_bg=rg_cv_ed
        else:
            if rg_cv_bg==-1:
                rg_cv_bg=rg_cv_ed=None
            else:
                rg_cv_ed=rg_cv_bg
            right_residue=0
        
        non_masked_len=0
        for e,rg in enumerate(masked_ranges[rg_cv_bg:rg_cv_ed+1]):
            if e>0:
                non_masked_len+=rg[0]-last_rg[1]
            last_rg=rg
        
        ends_residue=left_residue+right_residue
        non_masked_len=non_masked_len+ends_residue
        notOtherRepeat_score=non_masked_len
        
        notSimpRpt_score=ends_residue
        if rg_cv_bg==rg_cv_ed:
            if masked_families[rg_cv_bg] in('Simple_repeat'):#,'Low_complexity'
                notSimpRpt_score=ends_residue
        elif masked_families[rg_cv_bg]== masked_families[rg_cv_ed]=='Simple_repeat':
            if masked_matches[rg_cv_bg]==masked_matches[rg_cv_ed]:
                same=True
            else:
                rep_unit_l=masked_matches[rg_cv_bg].replace('(','').replace(')n','')
                rep_unit_r=masked_matches[rg_cv_ed].replace('(','').replace(')n','')
                if rep_unit_l > rep_unit_r:
                    rep_unit_longer=rep_unit_l;rep_unit_shorter=rep_unit_r
                else:
                    rep_unit_longer=rep_unit_r;rep_unit_shorter=rep_unit_l
                same=False
                remd=len(rep_unit_longer)%len(rep_unit_shorter)
                if remd == 0:
                    fold=len(rep_unit_longer)//len(rep_unit_shorter)
                    same = rep_unit_longer==rep_unit_shorter*fold
                    if same==False:
                        for i in [i+1 for i in range(len(rep_unit_shorter))]:
                            rep_unit_longer_i_shift=rep_unit_longer[i:]+rep_unit_longer[0:i]
                            same = rep_unit_longer_i_shift==rep_unit_shorter*fold
                            if same:
                                break
                    
            if same:
                left_inner_bonus_space=0
                right_inner_bonus_space=0
                left_inner_bonus=0
                right_inner_bonus=0
                i=0
                while rg_cv_bg+i<len(masked_matches)-1:
                    i+=1
                    if rg_cv_bg+i == rg_cv_ed or masked_matches[rg_cv_bg+i] in ('Simple_repeat','Low_complexity'):
                        left_inner_bonus_space=masked_ranges[rg_cv_bg+i][0]-masked_ranges[rg_cv_bg][1]
                        break
                i=0
                while rg_cv_ed-i>0:
                    i+=1
                    if rg_cv_ed-i == rg_cv_bg or masked_matches[rg_cv_ed-i] in ('Simple_repeat','Low_complexity'):
                        right_inner_bonus_space=masked_ranges[rg_cv_ed][0]-masked_ranges[rg_cv_ed-i][1]
                        break
                            
                if left_inner_bonus_space >= sig_min:
                    xl=masked_ranges[rg_cv_bg][1]-left
                    if xl < sig_min:
                        left_inner_bonus=min(left_inner_bonus_space,(90-xl))
                if right_inner_bonus_space >= sig_min:
                    xr=right-masked_ranges[rg_cv_ed][0]
                    if xr < sig_min:
                        right_inner_bonus=min(right_inner_bonus_space,(90-xr))

                possible_inner_residue_bonus =left_inner_bonus+right_inner_bonus
                notSimpRpt_score=ends_residue+possible_inner_residue_bonus
            
    return notOtherRepeat_score,notSimpRpt_score


#Make insertion position table
df_position=pd.read_csv(data_dir+'/'+cohort+'_ERVcaller_merge_siteinfo.txt',sep='\t',header=0)
df_position['pos_id']=df_position.loc[:,['chr','pos','alt_allele']].apply(
    lambda x: str(x[0])+':'+str(x[1])+':'+str(x[2]),axis=1).astype('string')
df_position['chrpos']=df_position[['chr','pos']].apply(lambda x: ':'.join([str(x[0]),str(x[1])]),axis=1)

print('Sites number of the original input file:',df_position.shape[0])


#plink outfiles chr name for X, Y, MT are 23, 24, 26. Borrow this for sorting and other potential future usage.
ser_temp=df_position['chr'].str.strip("chr")
ser_temp=ser_temp.where(ser_temp!='X', 23)
ser_temp=ser_temp.where(ser_temp!='Y', 24)
ser_temp=ser_temp.where(ser_temp!='MT', 26)
df_position['chr_num']=ser_temp.astype(int)
del ser_temp

#sort
df_position=df_position.sort_values(by=['chr_num','chr','pos'], 
                              ascending=[True, True, True],
                              ignore_index=True)

#merge tables. for 'aln_bg', 'aln_ed', 'max_aln_len', aln means on element.
df_position=df_position[['pos_id','chrpos','chr','chr_num','pos','direction','aln_bg', 'aln_ed', 
                         'max_aln_len','max_status', 'max_left_status','max_right_status','TE_acc']]
df_position=df_position.merge(df_freqANDhwe[['pos_id','alt_FREQ']],on='pos_id',how='left')#,'P_HET_EXCESS'

df_position=df_position.merge(df_family,left_on='TE_acc',right_on='acc',how='left')

ser_quan_not_other_repeat=df_position[['aln_bg','aln_ed','matchRange_on_element','matching_family','matching']].apply\
(lambda x: quan_not_other_repeats(x[0],x[1],x[2],x[3],x[4]), axis=1)
df_position['notOtherRepeat_score']=ser_quan_not_other_repeat.apply(lambda x: x[0])
df_position['notSimpRpt_score']=ser_quan_not_other_repeat.apply(lambda x: x[1])

df_position['element']=df_position.apply(lambda x: x['TE_acc'] if type(x['element'])!=str else x['element'],axis=1)
df_position['312_class'].fillna('XRV',inplace=True)
df_position['ERV_part'].fillna('XRV',inplace=True)
df_position['organism'].fillna('XRVhost',inplace=True)

print('Sites number before keeping only ERV1-3 classes and XRV:',df_position.shape[0])
#print('Before keeping there are these 312_classes:', df_position['312_class'].unique().tolist())

#filter inappropriate entries of positions
lis_RV_classes=['ERV3','MaLR','ERV1','ERV2','XRV']
df_position=df_position.loc[df_position['312_class'].isin(lis_RV_classes),:]
df_position=df_position.loc[~df_position['direction'].isna(),:]

df_position=df_position.reset_index(drop=True)
print('Sites number after keeping only ERV1-3 classes and XRV:',df_position.shape[0])

#assign data type
for col in ['pos','aln_bg','aln_ed','max_aln_len','max_status']:#,'bedtools_cluster_bg','bedtools_cluster_ed'
    try:
        df_position[col]=df_position[col].astype(int)
    except Exception as e:
        pass#print('Warning of not excuted:', e)

df_position.to_pickle(data_dir+'/'+cohort+'_ERVcaller_post_part1_indivmerge_site_info.pkl')

