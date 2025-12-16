# adherence explanation
# Author: Juan Delgado, PhD 
# Date: March 2024
# Documentation: This script generates a swimmer plot for systems forecasting.

import pandas as pd
import pandas_gbq as pdg

# post calulation of study_drugs compliance
study_drugs = pdg.read_gbq('select * from cactus.study_drugs', project_id='systemsforecasting')

# remove rows with missing dispensed date
study_drugs.dropna(subset=['dispensed_date'],inplace=True)
study_drugs = study_drugs.query('dispensed != "0"').reset_index(drop=True)

# remove rows with missing counted date
study_drugs['discrepancy'] = (study_drugs['discrepancy'].replace("nan",0)).astype(float)

# define units per day
units_per_day_dict = {'Dab':2,'Tram':1,'Niv':1,'Ipi':1}

# compute adherence
G = []
for gi,g in list(study_drugs.groupby(['PID', 'drug'])):
    if g.shape[0]>1 and not g['dispensed_date'].isna().all():
        g = g.sort_values(by='dispensed_date').reset_index(drop=True)
        g['days_between_visits'] = pd.to_datetime(g['dispensed_date']).diff().dt.days        
        g['tablets_per_day'] = g['drug'].map(units_per_day_dict)
        g['adherence'] = (g['discrepancy']/(g['days_between_visits'] * g['tablets_per_day'])).clip(0,1)
        print(g)
    G.append(g)
study_drugs = pd.concat(G,axis=0,ignore_index=True)
# this removes nans
study_drugs.query('counted_date==counted_date')

# save to bigquery
pdg.to_gbq(study_drugs, 'cactus.adherence', project_id='systemsforecasting', if_exists='fail')