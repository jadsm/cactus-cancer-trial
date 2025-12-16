# other_utils
# Author: Juan Delgado, PhD 
# Date: March 2024
# Documentation: This script contains supporting functions

import pandas as pd
from lifelines import KaplanMeierFitter
import numpy as np
from datetime import date,timedelta
from scipy.interpolate import interp1d
import statsmodels.api as sm



def count_week_days(d1,d2):
        daygenerator = (d1 + timedelta(x + 1) for x in range((d2 - d1).days)) # generate all days from d1 to d2
        days = sum(1 for day in daygenerator if day.weekday() < 5)
        return days

def get_first_descend(df,varname='pct_VAF_from_baseline',index = False):
        d,idx_out = [],[]
        for li,l in df.groupby('PID'):
                try:
                        iidx = np.where(l.loc[:,varname].diff()>=0)[0]
                        red = l.loc[np.array(l.index)[:iidx[-1]],:].reset_index(drop=True)
                        idx_out.append(np.array(l.index)[:iidx[-1]])
                except:
                        red = l.copy()
                d.append(red)
        if not index:
                return  pd.concat(d,axis=0,ignore_index = True)
        else:
                return np.concatenate(idx_out)
        
def get_apperances(df,sort_by = ['PID','days_from_randomization'],groupvar = 'PID',eventvar='event_outcome'):
        aa = df.sort_values(sort_by)
        bb = aa.copy()
        L = list(aa.groupby(groupvar))
        for li,l in L:
                previousri = l.index[0]
                for ri,row in l.iterrows():
                        if l.index[0] == ri:
                                continue
                        if row[eventvar] == l.loc[previousri,eventvar]:
                                bb.drop(index=ri,inplace=True)
                        previousri = ri
                
                # bb = aa.drop_duplicates(subset=['PID','event_outcome']).groupby(by=['PID','event_outcome']).cumcount()+1
                cc = bb.groupby(by=[groupvar,eventvar]).cumcount()+1
        aa = pd.concat([aa,cc],axis=1).rename(columns={0:'count'}).ffill()
        return aa


def kaplanmeier_fit(df, duration_col='duration', event_col='dead'):
        kmf = KaplanMeierFitter()
        kmf.fit(durations=df[duration_col], event_observed=df[event_col])

        # kmf.plot_survival_function()  # this is the function we're unpacking

        df_plot = kmf.survival_function_.copy(deep=True)
        df_plot['lower_bound'] = np.round(kmf.confidence_interval_['KM_estimate_lower_0.95'].values,2)
        df_plot['upper_bound'] = np.round(kmf.confidence_interval_['KM_estimate_upper_0.95'].values,2)
        df_plot.reset_index(inplace=True)
        # df_plot.head()
        return df_plot


def postproc_outcome(df):

        df['event_name'] = df['event_name'].str.capitalize()
        df = df.merge(df.query('event_name=="Drug"').groupby('PID')['event_date'].min().reset_index().rename(columns={'event_date':'Treatment_start_date'}),
                      on='PID')
        df['days_from_randomization'] = (pd.to_datetime(df['event_date'],utc=True)-pd.to_datetime(df['Randomisation_date'],utc=True)).dt.days
        df['days_from_treatment'] = (pd.to_datetime(df['event_date'],utc=True)-pd.to_datetime(df['Treatment_start_date'],utc=True)).dt.days
        # df['event_date'] = pd.to_datetime(df['event_date']).dt.strftime('%Y-%m-%d')

        mappings = {'Recist':'scan','Drug':'drug', 'Participant died':'Dead','Participant decision':'Withdrawal',
                    'Decrease in performance status':'Trial end','As per sa4 (fup period reduced and eot definition updated)':'Trial end','As per sa4 (fup period reduced and eot definition) updated':'Trial end'}
        df['event_type'] = df['event_name'].map(mappings).fillna('Adverse event').astype(str)
        df['event_outcome'] = df['event_outcome'].astype(str)
        df['death_reason'] = df['death_reason'].map({'Related to DP':'Melanoma'}).fillna('Other')
        idx = df['event_type'] == 'Dead'
        df.loc[idx,'event_type'] += ":"+df.loc[idx,'death_reason']

        idx = df.loc[:,'event_type'] == 'scan'
        df.loc[idx,'event_outcome'] = df.loc[idx,'event_outcome'].replace('nan','NE')

        # mapping2 = {"null":"Unknown or stable disease","CR":"Complete response","PR":"Partial response","SD":"Unknown or stable disease",
        # "PD":"Progressive disease","NE":"Unknown or stable disease"}# not evaluable

        # df.loc[idx,'event_outcome'] = df.loc[idx,'event_outcome'].map(mapping2).fillna("Unknown or stable disease")

        # simplify drugs
        drug_mapping = {'Dab':'Targeted', 'Tram':'Targeted', 'Niv':'IO', 'Ipi':'IO'}
        idx = df['event_type'] == 'drug'
        df.loc[idx,'event_outcome'] = df.loc[idx,'event_outcome'].map(drug_mapping).astype(str)
        return df

def get_treatment_classes(df):
        from constants import treatment_schedule_classes,first_line_treatment_classes
        # reverse the dictionary
        treatment_schedule_classes = {v: k for k in treatment_schedule_classes for v in treatment_schedule_classes[k]}
        df['treatment_classification']  = df['PID'].map(treatment_schedule_classes)
        # df['treatment_classification_num']  = df['treatment_classification'].map({'Targeted_to_immuno':"2",
        #                                                                         'Immuno_to_targeted':"3",
        #                                                                         'Targeted_only':"4",
        #                                                                         'Immuno_only':"5",
        #                                                                         'double_switch':"1"})
        first_line_treatment_classes = {v: k for k in first_line_treatment_classes for v in first_line_treatment_classes[k]}
        df['first_line_treatment']  = df['PID'].map(first_line_treatment_classes)


        df['treatment_classification_num']  = df['treatment_classification'].map({'T->Io':"2",
                                                                                'Io->T':"5",
                                                                                'T only':"1",
                                                                                'I only':"4",
                                                                                'T->Io->T':"3"})
        

def get_unique_drugs(df,df_drug,dates = ['Sample_Date','dispensed_date'],drugvar = 'drug'):
        df[drugvar] = ''

        for pid in df['PID'].unique():
                drug_now = df_drug.query('PID == @pid',inplace=False)
                df_now = df.query('PID == @pid')
                for ri,row in df_now.iterrows():
                        min_idx = (row[dates[0]]-drug_now[dates[1]]).abs().argmin()
                        df.loc[ri,drugvar] = drug_now.loc[drug_now.index[min_idx],drugvar]
        df[drugvar] = df[drugvar].apply(lambda x:"+".join(np.sort(np.unique(x.split('+')))))

def calculate_cutoff(df,timevar = 'days_from_baseline',variable='pct_VAF_from_baseline',threshold=20,othervars=['PID','Percentage_Mutant_VAF']):
        days = pd.DataFrame([],columns=['PID','days_cutoff'])
        for li,l in df.groupby('PID'):
                
                try:
                # if True:
                        # print(l.loc[:,['PID','Percentage_Mutant_VAF','pct_VAF_from_baseline','days_from_baseline','days_cutoff']])
                        # iidx = np.where(l.loc[:,'pct_VAF_from_baseline'].diff()>=0)[0]
                        iidx = np.where(l.loc[:,variable]<=threshold)[0][0]
                        # print(l.loc[:,'pct_VAF_from_baseline'].diff())
                        red = l.loc[np.array(l.index)[:iidx+1],[variable,timevar]+othervars].reset_index(drop=True)
                        # print(red)
                        # print("cutoff",np.interp(20, red[  'pct_VAF_from_baseline'], df.loc[red.index,'pct_VAF_from_baseline']))
                        d = red.loc[iidx-1,timevar]+(red.loc[iidx,timevar]-red.loc[iidx-1,timevar])*(red.loc[iidx-1,variable] - threshold)/(red.loc[iidx-1,variable] - red.loc[iidx,variable])
                        # d = red.loc[iidx-1,variable]+(-red.loc[iidx,timevar]+red.loc[iidx-1,timevar])*(red.loc[iidx-1,variable]-threshold)/(red.loc[iidx-1,variable] - red.loc[iidx,variable])
                        days_now = pd.DataFrame(np.array([li,d]).reshape(1,2),columns=['PID','days_cutoff'])
                        days = pd.concat([days,days_now],axis=0,ignore_index=True)
                        print('PID',li,d//1,'days')
                        # print(li,days)
                except:
                        print(li,'failed')
                        # days.append(np.nan)
                        days = pd.concat([days,pd.DataFrame(np.array([li,np.nan]).reshape(1,2),
                                          columns=['PID','days_cutoff'])],axis=0,ignore_index=True)
        
        return days


def get_baseline_or_screening(df):
    D = []
    g = list(df.groupby("PID"))
    for gi,gg in g:
        idx = (gg['Schedule'] == 'Baseline')
        if idx.any():
            D.append(gg.loc[idx,['PID','sample date','BRAF VAF ctDNA result']])
        else:
            idx = (gg['Schedule'] == 'Screening')
            D.append(gg.loc[idx,['PID','sample date','BRAF VAF ctDNA result']])
        # print(gi,np.where(idx))
    return pd.concat(D,axis=0,ignore_index=True).rename(columns={'sample date':'min_date','BRAF VAF ctDNA result':'Baseline_VAF'})

def add_dummy_points(df,idx,varcol='event_outcome',varcolidx = 6):
        L = list(df.loc[idx,:].groupby('PID'))
        for li,l in L:
                # if li != 6:
                #         continue
                l = l.sort_values('days_from_treatment').reset_index(drop=False)
                # df.loc[l.loc[0,'index'],'event_outcome'] = 'Baseline'
                # l.loc[0,'event_outcome'] = 'Baseline'
        
                aux = pd.DataFrame([],columns=l.columns)
                for i in range(1,l.shape[0]):
                        
                        if l.loc[i-1,varcol] != l.loc[i,varcol]: # if there is a change add a new point
                                new_point = l.loc[i-1,:].copy()
                                new_point['days_from_treatment'] = l.loc[i,'days_from_treatment']
                                new_point['dummy'] = True
                                aux = pd.concat([aux,pd.DataFrame(new_point).T],axis=0,ignore_index=True)
                                
                        # if dead or withdrawal last recist
                        if l.loc[i,'event_type'] in ('Withdrawal','Dead:Melanoma','Dead:Other','Trial end'):
                                df.loc[l.loc[i,'index'],varcol] = l.loc[i-1,varcol]
                                l.loc[i,varcol] = l.loc[i-1,varcol]
                                aux.iloc[-1,varcolidx] = l.loc[i-1,varcol]
                                # continue
                df = pd.concat([df,aux],axis=0,ignore_index=True)
        return df


def calculate_simple_visits(df,bins = [0,20,35,50,70,90,110,130,190,250,350]):
    df['simple_visits'] = pd.cut(df['days_from_baseline'], bins=bins, labels=[str('0'*(2-len(str(i))))+str(i) for i in range(1,len(bins))]).astype(str)
    df.loc[df['days_from_baseline']<=bins[0],'simple_visits'] = '0'
    df.loc[df['days_from_baseline']>=bins[-1],'simple_visits'] = '0'
    return df

def enrich_sld(df_sld,dfaux,df_scan,bins):
    df_sld = df_sld.merge(dfaux,on='PID')
    df_sld['days_from_baseline'] = (pd.to_datetime(df_sld['date_combined'],utc=True) - pd.to_datetime(df_sld['min_date'],utc=True)).dt.days
    df_sld = calculate_simple_visits(df_sld,bins=bins)
    return df_sld

def enrich_group_sld(df_sld,dfaux,df_scan,bins):
    df_sld = enrich_sld(df_sld,dfaux,df_scan,bins)
    grouvars = ['Treatment_arm','simple_visits']
    return df_sld.groupby(grouvars)['mm_SLD'].describe().reset_index().round(1)

def enrich_scan(df_scan,dfaux,bins):
    df_scan = df_scan.merge(dfaux,on='PID')
    df_scan['days_from_baseline'] = (pd.to_datetime(df_scan['scan_date'],utc=True) - pd.to_datetime(df_scan['min_date'],utc=True)).dt.days
    df_scan = calculate_simple_visits(df_scan,bins=bins)
    aux = df_scan.groupby(['simple_visits','Overall_response'])['PID'].nunique().reset_index()
    aux2 = aux.groupby('simple_visits')['PID'].sum().reset_index()
    aux = aux.merge(aux2,on='simple_visits')
    aux['patient_pct'] = (aux['PID_x']/aux['PID_y']*100).round(0)
    return df_scan,aux



def find_nn(df_out):
    print(df_out['PID'].unique())
    f = interp1d(df_out.dropna(subset=['LDH']).loc[:,'days_from_baseline'],
            df_out.dropna(subset=['LDH']).loc[:,'LDH'],
            kind = 'nearest',fill_value="extrapolate")

    f2 = interp1d(df_out.dropna(subset=['mm_SLD']).loc[:,'days_from_baseline'],
            df_out.dropna(subset=['mm_SLD']).loc[:,'mm_SLD'],
            kind = 'nearest',fill_value="extrapolate")
    newdf = df_out.dropna(subset=['Percentage_Mutant_VAF']).reset_index(drop=True)
    newdf['LDH_interp'] = f(newdf.loc[:,'days_from_baseline'])
    newdf['mm_SLD_interp'] = f2(newdf.loc[:,'days_from_baseline'])
    return newdf


def shift_first_scan_to_treatment_start(df):
        idx = df['event_type'] == 'scan'
        min_treatment_date = df.loc[idx,:].groupby('PID')['days_from_treatment'].min().reset_index().rename(columns={'days_from_treatment':'first_day_treatment'})
        min_treatment_date = df.loc[idx,:].merge(min_treatment_date,on='PID')
        iidx = min_treatment_date['first_day_treatment'] == min_treatment_date['days_from_treatment']
        min_treatment_date.loc[iidx,'days_from_treatment'] = 0.0
        min_treatment_date.drop(columns=['first_day_treatment'],inplace=True)
        df = pd.concat([df.loc[~idx,:],min_treatment_date],axis=0,ignore_index=True)
        return df

def transform_adverse_events(df_ae_full,df_drugs):
    
# remove dates with switched onset and resolved dates
    
    idx = np.logical_or(df_ae_full.loc[:,'Date_resolved'].isna(), df_ae_full.loc[:,'Date_resolved'] >= df_ae_full.loc[:,'Onset_date'])
    df_ae_full = df_ae_full.loc[idx,:].reset_index(drop=True)
    
    # decode
    df_ae_full.sort_values(by=['Onset_date'],inplace=True)
    df_ae_full['Event_grade'] = df_ae_full['Event_grade'].fillna(1.)
    df_ae_full['AE_name'] = df_ae_full['Event_name'].fillna(df_ae_full['AE_Code'])
    df_ae_full['Dead'] = df_ae_full['Outcome'] == 'Death'
    
    dfaux = df_drugs.groupby('PID')['event_date'].min().reset_index().rename(columns={'event_date':'treatment_start_date'})
    df_ae_full = dfaux.merge(df_ae_full,on='PID')
    df_ae_full['Onset_days'] = (pd.to_datetime(df_ae_full['Onset_date'],utc=True) - pd.to_datetime(df_ae_full['treatment_start_date'],utc=True)).dt.days
    df_ae_full['Resolved_days'] = (pd.to_datetime(df_ae_full['Date_resolved'],utc=True) - pd.to_datetime(df_ae_full['treatment_start_date'],utc=True)).dt.days

    # Onset days need to be positive
    df_ae_full = df_ae_full.query('Onset_days >= 0').reset_index(drop=True)

    # linearise but wthout the crosstab
    # d = pd.crosstab(df_ae_full['Onset_date'],columns=df_ae_full['Event_grade']).cumsum()
    # d = d.stack().reset_index().rename(columns={0:'count','Onset_date':'date'}).sort_values(by=['Event_grade','date'])
    d = df_ae_full.loc[:,['Onset_days','Event_grade','PID']].rename(columns={'Onset_days':'days'})
    d['count'] = 1
    d['operation'] = 'add'

    # d2 = pd.crosstab(df_ae_full['Date_resolved'],columns=df_ae_full['Event_grade']).cumsum()
    # d2 = d2.stack().reset_index().rename(columns={0:'count','Date_resolved':'date'}).sort_values(by=['Event_grade','date'])
    d2 = df_ae_full.loc[:,['Resolved_days','Event_grade','PID']].rename(columns={'Resolved_days':'days'})
    d2['count'] = -1
    d2['operation'] = 'remove'

    dd = pd.concat([d,d2],axis=0,ignore_index=True).dropna(subset='days')
    dd['weeks'] = dd['days'].astype(float)//7
    
    dd = dd.groupby(['Event_grade','days','operation','PID'])['count'].sum().reset_index()

    dd = dd.sort_values(by=['Event_grade','days']).reset_index(drop=True)#.query('days >=0')
    dd['Event_grade'] = dd['Event_grade'].astype(int)
    
#     dd = dd.sort_values(by=['days']).reset_index(drop=True)
    
#     dd['cumsum'] = dd['count'].cumsum()
#     dd.groupby(['Event_grade','days'])['PID'].nunique().cumsum().reset_index()

    d = []
    for gi,g in list(dd.groupby('Event_grade')):
        patient = []
        g['cumsum'] = g['count'].cumsum()
        g['cumpid'] = 0
        for ri,row in g.iterrows():
               if row['operation'] == 'add' and not row['PID'] in patient: # only if it doesnt exist count it - if not, pass
                        patient.append(row['PID'])
               elif row['operation'] == 'remove' and row['PID'] in patient: # if exists remove it
                        patient.remove(row['PID'])
                        #     g.loc[ri,'cumpid'] = - 1 + g.loc[ri-1,'cumpid']                
                      
               g.loc[ri,'cumpid'] = len(patient)
                                              
        d.append(g.loc[:,['Event_grade', 'days', 'cumsum','cumpid']])
    d = pd.concat(d,axis=0,ignore_index = True)
    d = d.query('days >=0').reset_index(drop=True)
    d = d.query('cumsum > 0').reset_index(drop=True)
    d['days'] = d['days'].astype(str)

    df_drugs['drug_simple'] = df_drugs['event'].map({'Dab':'Targeted','Tram':'Targeted','Niv':'IO','Ipi':'IO'})
    df_drugs = dfaux.merge(df_drugs,on='PID')
    df_drugs['days_from_treatment_start'] = (pd.to_datetime(df_drugs['event_date'],utc=True) - pd.to_datetime(df_drugs['treatment_start_date'],utc=True)).dt.days
    df_drugs = df_drugs.loc[:,['drug_simple','days_from_treatment_start','PID']].drop_duplicates().reset_index(drop=True)

    bins = range(0,750,30)
    df_drugs['simple_days'] = pd.cut(df_drugs['days_from_treatment_start'], 
                                bins=bins, 
                                labels=range(30,750,30)).astype(float).fillna(750)
    df_drugs.loc[df_drugs['days_from_treatment_start']<=0,'simple_days'] = 0
    df_drugs['simple_days'] = df_drugs['simple_days'].astype(int)

    df_drugs = df_drugs.loc[:,['drug_simple','simple_days','PID']].drop_duplicates().reset_index(drop=True)
    # prolong treatment if has not changed until last date available for each patient
    dx = pd.pivot_table(df_drugs, values='drug_simple', index='PID', columns='simple_days',aggfunc='max')
    
    ci = 0
    for ri,row in dx.iterrows():
        print(ri,row)
        idx = np.where(row.notna())[0][-1]
        dx.iloc[ci,:idx] = dx.iloc[ci,:idx].fillna(method='ffill')
        dx.iloc[ci,idx+1:] = 'End'
        ci +=1
    dx = dx.reset_index() 
    dx = dx.melt(id_vars = ['PID'],value_vars = dx.loc[:,0:].columns,value_name = 'drug')
    return dx,d

def calculate_CC_summary_stats(D2):    
    D2['BRAF_SLD_CC'] = sm.tsa.stattools.ccf(D2.loc[:,'Percentage_Mutant_VAF'].values, D2.loc[:,'mm_SLD_interp'].values, adjusted=False)
    D2['BRAF_LDH_CC'] = sm.tsa.stattools.ccf(D2.loc[:,'Percentage_Mutant_VAF'].values, D2.loc[:,'LDH_interp'].values, adjusted=False)
    G = D2.groupby('PID')
    aa = G['BRAF_SLD_CC'].describe().reset_index(drop=False)
    bb = G['BRAF_LDH_CC'].describe().reset_index(drop=False)
    aa['var'] = 'ctDNA-SLD'
    bb['var'] = 'ctDNA-LDH'
    # aa.rename(columns={k:k+"_SLD" for k in aa.keys()[1:]},inplace=True)
    # bb.rename(columns={k:k+"_LDH" for k in bb.keys()[1:]},inplace=True)
    # df_cc = aa.merge(bb,on='PID')
    df_cc = pd.concat([aa,bb],axis=0,ignore_index = True)

    df_cc = df_cc.merge(D2.loc[:,['PID','Treatment_arm']].drop_duplicates(),on='PID')
    # df_cc.loc[:,['PID','mean_SLD','mean_LDH']]
    return df_cc
