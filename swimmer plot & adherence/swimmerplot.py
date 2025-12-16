# swimmer plot
# Author: Juan Delgado, PhD 
# Date: March 2024
# Documentation: This script generates a swimmer plot for systems forecasting.

import altair as alt
import pandas as pd
from other_utils import *

def plot_outcomes(df, title, legend=True, colors_outcome=['green','orange','orange','red','red','gray','gray']):
    """
    Generate a swimmer plot for outcomes.

    Parameters:
    - df (pandas.DataFrame): The input DataFrame containing the data.
    - title (str): The title of the plot.
    - legend (bool): Whether to include a legend in the plot. Default is True.
    - colors_outcome (list): The list of colors for the outcome categories. Default is ['green','orange','orange','red','red','gray','gray'].

    Returns:
    - chart (altair.Chart): The swimmer plot chart.
    """
        
    base = alt.Chart(df).transform_calculate(
        order='datum.treatment_classification_num'
    ).encode(
        alt.X(f"days_from_treatment:Q").title("days from treatment start"),
        alt.Y("PID:N", sort=alt.SortField('order')).axis(offset=0, ticks=False, minExtent=0, domain=False).title('patient ID')
    )
        
    palette = ["#138086",'#DC8665'] if df.query('event_type == "drug"')['event_outcome'].nunique() == 2 else ['#DC8665']
    colourset = alt.Color('event_outcome').scale(range=palette).legend(title='Treatment', orient='bottom') if legend else alt.Color('event_outcome').scale(range=palette).legend(None)
        
    line0 = base.mark_line().encode(
        detail="PID",
        color=colourset,
        opacity=alt.value(.7),
        yOffset=alt.value(15),
        size=alt.value(7)
    ).transform_filter(alt.FieldEqualPredicate(field='event_type', equal='drug'))
        
    colourset = alt.Color('event_outcome').scale(range=colors_outcome).legend(title='Outcome', orient='bottom') if legend else alt.Color('event_outcome').scale(range=colors_outcome).legend(None)
    line = base.mark_line().encode(
        detail="PID",
        size=alt.value(17),
        color=colourset,
        opacity=alt.value(.6),
    ).transform_filter(alt.FieldEqualPredicate(field='event_type', equal='scan'))

    point = base.mark_point(filled=False).encode(
        color=alt.value("black"),
        size=alt.value(100),
        opacity=alt.value(.8),
        shape=alt.Shape('event_type').legend(title='Events', orient='bottom'),
        tooltip=['Treatment_arm', 'PID', 'event_name', 'days_from_treatment','days_from_randomization', 'event_grade',
        'event_outcome', 'event_type']
    ).transform_filter(alt.FieldOneOfPredicate(field='event_type', oneOf=['Dead:Melanoma','Dead:Other','Withdrawal','Trial end'])).transform_filter(alt.FieldEqualPredicate(field='dummy', equal=False))
        
    ticks = base.mark_point(filled=True).encode(
        color=alt.value("black"),
        size=alt.value(25),
        shape=alt.value("circle"),  # stroke
        opacity=alt.value(.8),
        tooltip=['Treatment_arm', 'PID', 'event_name', 'days_from_treatment','days_from_randomization',
        'event_outcome', 'event_type']
    ).transform_filter(alt.FieldEqualPredicate(field='event_type', equal='scan'))

    reference_line = alt.Chart(pd.DataFrame({'x':[0]})).mark_rule(color='black', size=1).encode(x='x')

    chart = (line + point + ticks + reference_line + line0).properties(
        width=650,
        height=850,
        title=title
    ).resolve_scale(color='independent')
    
    return chart

### main code
# load data
df = pd.read_csv('outcome.csv')

# perform some data cleaning
df['event_date'] = df['event_date'].apply(lambda x: str(x).split(" ")[0])
df = df.drop_duplicates().reset_index(drop=True)
df = postproc_outcome(df)

# add dummy points for events
df['dummy'] = False
df = df.merge(df.groupby('PID')['days_from_treatment'].apply(lambda x: x.max() - x.min()).reset_index(drop=False).rename(columns={'days_from_treatment':'total_duration'}),
            on='PID', how='left')

aux = get_apperances(df.query('event_type == "scan"'), sort_by=['PID','days_from_treatment'], groupvar='PID', eventvar='event_outcome')
df = pd.concat([df, aux.loc[:, 'count']], axis=1)
df['event_outcome'] = (df['event_outcome'] + " " + df['count'].astype(str).replace('nan','')).str.replace('.0','')

idx = df.loc[:, 'event_type'].apply(lambda x: x in ('scan','Withdrawal','Dead:Melanoma','Dead:Other','Trial end')) 
df = add_dummy_points(df, idx, varcol='event_outcome', varcolidx=6)

idx = df.loc[:, 'event_type'].apply(lambda x: x in ('drug','Withdrawal','Dead:Melanoma','Dead:Other','Trial end')) 
df = add_dummy_points(df, idx, varcol='event_outcome', varcolidx=6)

# get treatment classes
get_treatment_classes(df)

# remove trailing _
df['event_outcome'] = df['event_outcome'].str.rstrip('_')

# remove columns
df.drop(columns=['event_date','Randomisation_date','Randomisation_date','level_0','index'], inplace=True)

# make first scan coincide with start of treatment
df = shift_first_scan_to_treatment_start(df)

# remove duplicates
idx = df['event_type'] == 'drug'
aux = df.loc[idx, :].drop_duplicates(subset=['Treatment_arm','PID','event_outcome','days_from_treatment']).dropna(subset=['days_from_treatment'])
df = pd.concat([df.loc[~idx, :], aux], axis=0).reset_index()

# generate swimmer plot
chart = alt.hconcat(plot_outcomes(df.query('Treatment_arm == "ARM A"'), "ARM A", colors_outcome=['green','gray','red','red','orange','orange','lightblue','lightblue']),
                    plot_outcomes(df.query('Treatment_arm == "ARM B"'), "ARM B", legend=False, colors_outcome=['gray','red','orange','orange','lightblue'])).configure_axis(labelFontSize=18,
                                                                                                titleFontSize=18
        ).configure_title(fontSize=20, 
        anchor="middle").configure_legend(titleColor='black', 
        titleFontSize=16, labelFontSize=16).interactive()
chart.save("timeline_swimmer.html")
