-- recist response
-- event
-- withdrawal
-- this is the data behind the swimmer plot

create or replace view `systemsforecasting.cactus.outcome`
as
(select T.*,c.Randomisation_date from
(select * from
 (SELECT
    Treatment_arm,
    PID,
    Main_Reason event_name,
    Death_reason death_reason, 
    cast(Withdrawal_date as string) event_date,
    "end" event_grade,
    "end" event_outcome
  FROM
    `systemsforecasting.cactus.patient`
  LEFT JOIN
    `systemsforecasting.cactus.withdrawal`
  USING
    (PID)
  WHERE
    Treatment_arm IS NOT NULL ) 
    union all
(SELECT Treatment_arm,
    PID,
    Event_name event_name,
    "" death_reason,
    cast(Onset_date as string) event_date,
    Event_grade event_grade,
    Outcome event_outcome
  FROM
    `systemsforecasting.cactus.events`
    JOIN 
    `systemsforecasting.cactus.patient`
    USING 
        (PID)
  WHERE
    SAE_Note_No IS NOT NULL
    AND Serious = 'Yes' 
    and Event_name is not null) 
union all
(select Treatment_arm,
PID,
"recist" event_name,
"" death_reason,
cast(ifnull(ifnull(CT_scan_date,MRI_scan_date),PET_CT_scan) as string) event_date,
Null event_grade,
Overall_response event_outcome
from `systemsforecasting.cactus.recist`
JOIN 
    `systemsforecasting.cactus.patient`
    USING 
        (PID)
 )
 union all
 (select Treatment_arm,
PID,
"drug" event_name,
"" death_reason,
cast(dispensed_date as string) event_date,
Null event_grade,
drug event_outcome
from `systemsforecasting.cactus.study_drugs`
JOIN 
    `systemsforecasting.cactus.patient`
    USING 
        (PID)
 )
 order by PID
 ) T
join `systemsforecasting.cactus.consent` c
using(PID)
)
