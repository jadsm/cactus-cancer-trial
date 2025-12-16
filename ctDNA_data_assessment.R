# load libraries
library(tidyverse)  # read and elaborate data
library(haven)      # read .dta files
library(lubridate)  # operations on dates
library(ggplot2)
library(survminer)
library(survival)
library(dplyr)
# This script will go through collection of baseline characteristics 
# and do the plots relating to the primary end-point around
# ctDNA and the KM plots as per the secondary analysis

# Demographics -----------------------------------------------------------

# we will go through each table and select what we need
d1<-read.csv("../Data/CSV_export/CAcTUS.Table 1.csv")
colnames(d1)<-d1[1,]
d1<-d1[-1,]
colnames(d1)
# this data-set tells us the eligibility of patients to the trial
# from this file one can see who got onto the trail and reasons why others
# were not recruited - this is the file that would be sued to create the consort
# our interest is in only those that entered the study
d1$`failed Inc number`[d1$`Eligible for study`=="00|No"]
d1$`Failed Ex number`[d1$`Eligible for study`=="00|No"]
# patients were excluded mainly due to failure on inclusion criteria and reasons
# were criteria 8 and 9 - from protocol:
# 8.	Screening ctDNA (as defined by the mutant BRAF VAF in plasma) ≥1.5% 
# apparently it started off with 5% and then changed to 1.5%
# 9.	Adequate organ function (see table 2.)
  

# interestingly all had ctDNA measurements
d1$`BRAF VAF ctDNA result`
# 888.88 are NA
d1$`BRAF VAF ctDNA result`[d1$`BRAF VAF ctDNA result`=="888.88"]<-NA
d1$`If no reason ctDNA BRAF VAF less 5`[d1$`Eligible for study`=="00|No"]
# so quite a few were rejected due to VAF levels being too low so <5%
d1$`If no reason ctDNA BRAF VAF less 5`[is.na(d1$`BRAF VAF ctDNA result`)]
# so teh two who did not have ctDNA the reason for not including was not ctDNA
d1$`failed Inc number`[is.na(d1$`BRAF VAF ctDNA result`)]

d1$`BRAF VAF ctDNA result`<-as.numeric(d1$`BRAF VAF ctDNA result`)
d1[,c("PID","BRAF VAF ctDNA result")]
ggplot(d1,aes(`Eligible for study`,`BRAF VAF ctDNA result`,
              col=`Eligible for study`))+geom_jitter()+
  theme_bw(base_size=14)+scale_y_log10()+ylab("BRAF VAF ctDNA [%]")+
  geom_hline(yintercept = 1.5,col=1,lty=2)+
  geom_hline(yintercept = 5,col=1,lty=3)

# Plotting Cumulative Density
ggplot(d1, aes(x = `BRAF VAF ctDNA result`)) +
  stat_ecdf(geom = "step", pad = FALSE) + # pad=FALSE ensures the plot starts at the first data point
  theme_bw(base_size=14) +
  labs(x = "BRAF VAF ctDNA [%]",
       y = "Cumulative Fraction")+
  scale_x_log10(breaks=c(0.1,0.2,0.5,1.5,5,10,100),lim=c(0.05,100))+
  scale_y_continuous(breaks=c(seq(0,1,by=0.1)))



d1$`failed Inc number`[d1$`Eligible for study`=="00|No" & 
                         d1$`ctDNA blood sample taken`=="01|Yes"]
# 5 had organ function issues
d1$`BRAF VAF ctDNA result`[d1$`Eligible for study`=="00|No" &
                            c(is.na(d1$`failed Inc number`))]
# the 9 NAs all had ctDNA<1.5
d1$`BRAF VAF ctDNA result`[d1$`Eligible for study`=="00|No" &
                             c(is.na(d1$`failed Inc number`)==F & 
                                 c(d1$`failed Inc number`=="8"|
                                     d1$`failed Inc number`=="INCLUSION #8"|
                                     d1$`failed Inc number`=="NO.8"))]
# those numbered as ctDNA failures 2 did have ctDNA >1.5
d1$`BRAF VAF ctDNA result`[d1$`Eligible for study`=="00|No" &
                             c(is.na(d1$`failed Inc number`)==F & 
                                 c(d1$`failed Inc number`=="9"|
                                     d1$`failed Inc number`=="INCLUSION 9"|
                                     d1$`failed Inc number`=="9 CARDIAC LVEF"|
                                     d1$`failed Inc number`=="INCLUSION CRITERIA 9"))]
d1$`BRAF VAF ctDNA result`[d1$`Eligible for study`=="00|No" &
                             c(is.na(d1$`failed Inc number`))]

# so 4 out of the 6 did meet ctDNA criteria but had organ issues
# so we go from 46 to 21 in the following way
# 2 had organ issues and so never got ctDNA
# 6 more had organ issues but did get ctDNA
# 8 were classified as failing ctDNA of those 2 did have >1.5 based on refiend criteria
# we will still note them as ctDNA fails
# 9 have no reason for exclusion but did have very low ctDNA<1
# so thats where the 25 come from

# we shall now keep the 21
d1<-d1[d1$`Eligible for study`=="01|Yes",]
# columns we shall keep
colnames(d1)
d1<-d1[,c("PID",
          "Dab and Tram immediately",
          "Rapidly progressing disease",
          "Date of randomisation",
          "Treatment arm",
          "Gender",
          "Date of birth",
          "BRAF VAF ctDNA result",
          "Point mutation",
          "T Classification entry",
          "N Classification entry",
          "M Classification entry",
          "Stage",
          "Unresectable",
          "No of disease sites",
          "Liver",
          "Lung",
          "Subcutaneous",
          "Brain",
          "Nodal",
          "Abdominal",
          "Bone",
          "Other")]


# we will keep these in case they start asking for adjusted analyses
# we know where the labs screening data is and can go and get those as required
# next step is to collect progression and death times - juan has shown us where 
# these are so will go get them...


# KM plots ----------------------------------------------------------------


d19<-read.csv("../Data/CSV_export/CAcTUS.Table 19.csv")
# has death dates and right-censoring for OS
d6<-read.csv("../Data/CSV_export/CAcTUS.Table 6_RECIST.csv")
# this has progression times so a 1st and 2nd one
# we are meant to compare arms for 1st PFS and then 2nd PFS

# lets do OS first its easier
colnames(d19)<-d19[1,]
d19<-d19[-1,]
colnames(d19)

# its the 1st 21 patients we want data for
d19<-d19[c(1:21),]
# 4 column are most useful
d19<-d19[,c(1:3,7)]

# merge with d1
d1<-merge(d1,d19,by=c("PID"),all.x=T,all.y=F)

# calculate OS times and outcomes
d1$`Date of randomisation`<-as.Date(d1$`Date of randomisation`,"%d/%m/%Y")
d1$`Date of withdrawal`<-as.Date(d1$`Date of withdrawal`,"%d/%m/%Y")

d1$OS<-as.numeric(d1$`Date of withdrawal`-d1$`Date of randomisation`)
d1$DEATH<-0
d1$DEATH[d1$`Main Reason`=="03|Participant died"]<-1


# set options for plotting survival curves
font.lab     <- 14
font.title   <- 16
font.tick    <- 12
conf.int     <- T
risk.table   <- T
ncensor.plot <- T
palette      <- c('orange','#558DAE')
ggtheme      <- theme_bw()
risk.table.height   <- 0.2
ncensor.plot.height <- 0.15

d1$ARM<-d1$`Treatment arm`
p.value <- surv_pvalue(fit = survfit(Surv(OS, DEATH) ~ ARM, data=d1), data = d1)
# plot OS first
ggsurv<-ggsurvplot(
  survfit(Surv(OS, DEATH) ~ ARM, data=d1),
  conf.int = conf.int,
  pval = round(p.value$pval,3),
  risk.table = risk.table,
  ncensor.plot = ncensor.plot,
  palette = palette,
  risk.table.col = "strata",
  legend.labs =
    c("ARM A","ARM B"),
  ggtheme = ggtheme,
  title = "Overall Survival",
  ylab = "Survival Fraction",
  xlab = "Time from Randomisation (Days)",
  font.title = c(font.title),
  font.x = c(font.lab),
  font.y = c(font.lab),
  break.time.by = 100,
  font.tickslab = c(font.tick),
  risk.table.height = risk.table.height,
  ncensor.plot.height = ncensor.plot.height
)
ggsurv$plot <- ggsurv$plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 70), "points"))
ggsurv$table <- ggsurv$table + 
  theme(
    plot.margin = unit(c(10, 10, 10, 5), "points"))
ggsurv$ncensor.plot <- ggsurv$ncensor.plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 95), "points"))
ggsurv

# we will look at cancer specific death and this is done 
# by right censoring the COVID death
d1$CDEATH<-d1$DEATH
d1$CDEATH[d1$`Death reason`=="77|Other specify"]<-0
p.value <- surv_pvalue(fit = survfit(Surv(OS, CDEATH) ~ ARM, data=d1), data = d1)
ggsurv<-ggsurvplot(
  survfit(Surv(OS, CDEATH) ~ ARM, data=d1),
  conf.int = conf.int,
  pval = round(p.value$pval,3),
  risk.table = risk.table,
  ncensor.plot = ncensor.plot,
  palette = palette,
  risk.table.col = "strata",
  legend.labs =
    c("ARM A","ARM B"),
  ggtheme = ggtheme,
  title = "Cancer Specific Survival",
  ylab = "Survival Fraction",
  xlab = "Time from Randomisation (Days)",
  font.title = c(font.title),
  font.x = c(font.lab),
  font.y = c(font.lab),
  break.time.by = 100,
  font.tickslab = c(font.tick),
  risk.table.height = risk.table.height,
  ncensor.plot.height = ncensor.plot.height
)



ggsurv$plot <- ggsurv$plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 70), "points"))
ggsurv$table <- ggsurv$table + 
  theme(
    plot.margin = unit(c(10, 10, 10, 5), "points"))
ggsurv$ncensor.plot <- ggsurv$ncensor.plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 95), "points"))
ggsurv
# next we do PFS
# we have two PFS options to consider the first PFS and the second
# lets look at the 1st PFS initially - so this is the 1st time a patient 
# experiences a PD
colnames(d6)<-d6[1,]
d6<-d6[-1,]
colnames(d6)

# we have 3 scan methods - lets quickly look to see what method was used for each
unique(d6$`TL1 Scan method`)
unique(d6$`TL2 Scan method`)
unique(d6$`TL3 Scan method`)
unique(d6$`TL4 Scan method`)
unique(d6$`TL5 Scan method`)
# so MRI date is no use we shall look at CT and PET CT
d6$`CT scan date`
# 8888-88-88 is missing value
d6$`CT scan date`[d6$`CT scan date`=="8888-88-88"]<-NA
d6$`PET CT scan`[d6$`PET CT scan`=="8888-88-88"]<-NA
# lets check the dates of CT scan v PET-CT
d6[,c("CT scan date","PET CT scan")]
# looking at this we can see row 39 is odd as we have a CT scan date and PET-CT
# date that are almost 3 months apart
# if we look at them...
d6[38,]
# although they had PET CT the data used was from CT so we shall stick with CT
# date

d6$Date<-d6$`CT scan date`
# when these are missing we use PET CT
d6$Date[is.na(d6$Date)]<-d6$`PET CT scan`[is.na(d6$Date)]
# all i want is...
d6<-d6[,c("PID","Date","Overall response")]
d6$Date<-as.Date(d6$Date,"%d/%m/%Y")

# collect 1st PD or right-censored
result <- d6 %>%
  arrange(PID, Date) %>%
  group_by(PID) %>%
  summarise(
    First_PD_or_Last_Date = if(any(`Overall response` == "04|PD")) {
      min(Date[`Overall response` == "04|PD"])
    } else {
      last(Date)
    },
    PD_Occurred = any(`Overall response` == "04|PD") 
  )

# now add to the main data-set
d1<-merge(d1,result,by=c("PID"),all.x=T,all.y=F)
colnames(d1)[32]<-"PFS"
d1$PFS<-as.numeric(d1$PFS)
d1$PDAY<-as.numeric(d1$First_PD_or_Last_Date-d1$`Date of randomisation`)
# for those that did not have PD we takes the OS times#
# note death is a progression event so we add death
d1$PDAY[d1$PFS==0]
d1$OS[d1$PFS==0]
d1$DEATH[d1$PFS==0]
# two died before radiologicaly progressing so they join PFS
d1$PDAY[d1$PFS==0]<-d1$OS[d1$PFS==0]
d1$PFS[d1$PFS==0]<-d1$DEATH[d1$PFS==0]
# right we now have 1st progression
p.value <- surv_pvalue(fit = survfit(Surv(PDAY, PFS) ~ ARM, data=d1), data = d1)

ggsurv<-ggsurvplot(
  survfit(Surv(PDAY, PFS) ~ ARM, data=d1),
  conf.int = conf.int,
  pval = round(p.value$pval,3),
  risk.table = risk.table,
  ncensor.plot = ncensor.plot,
  palette = palette,
  risk.table.col = "strata",
  legend.labs =
    c("ARM A","ARM B"),
  ggtheme = ggtheme,
  title = "Progression-Free Survival 1",
  ylab = "Fraction Progression Free",
  xlab = "Time from Randomisation (Days)",
  font.title = c(font.title),
  font.x = c(font.lab),
  font.y = c(font.lab),
  break.time.by = 100,
  font.tickslab = c(font.tick),
  risk.table.height = risk.table.height,
  ncensor.plot.height = ncensor.plot.height
)
ggsurv$plot <- ggsurv$plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 70), "points"))
ggsurv$table <- ggsurv$table + 
  theme(
    plot.margin = unit(c(10, 10, 10, 5), "points"))
ggsurv$ncensor.plot <- ggsurv$ncensor.plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 95), "points"))
ggsurv

# they did ask for PFS rate at 1 year too
km<-survfit(Surv(PDAY, PFS) ~ ARM, data=d1)
summary(km,times = c(365))

survfit(formula = Surv(PDAY, PFS) ~ ARM, data = d1)

# ARM=01|ARM A 
# time       n.risk      n.event     survival      std.err lower 95% CI upper 95% CI 
# 365.000        3.000        4.000        0.600        0.155        0.362        0.995 
# 
# ARM=02|ARM B 
# time       n.risk      n.event     survival      std.err lower 95% CI upper 95% CI 
# 365.0000       2.0000       7.0000       0.2694       0.1511       0.0897       0.8087 

# we also have PFS2 defined as:
# iii.	Second Progression Free Survival (PFS) defined as the time interval from 
# randomisation until confirmed second disease progression according to RECIST v1.1 
# or death, whichever occurs first.



result <- d6 %>%
  arrange(PID, Date) %>%
  group_by(PID) %>%
  mutate(
    PD_Flag = `Overall response` == "04|PD", # Create a logical flag for PD events
    PD_Cumulative = cumsum(PD_Flag) # Cumulative sum of PD events
  ) %>%
  filter(PD_Flag & PD_Cumulative == 2) %>% # Filter for the second PD event
  summarise(
    PD_Occurred = any(`Overall response` == "04|PD"),
    Second_PD_Date = first(Date) # Get the date of the second PD event if it exists
  )

colnames(result)[2]<-"PFS2"
result$PFS2<-as.numeric(result$PFS2)

d1<-merge(d1,result, by=c("PID"),all.x=T,all.y=F)
d1$PFS2[is.na(d1$PFS2)]<-0
d1$Second_PD_Date[d1$PFS2==0]<-d1$`Date of withdrawal`[d1$PFS2==0]
d1$PDAY2<-as.numeric(d1$Second_PD_Date-d1$`Date of randomisation`)
d1$PFS2[d1$PFS2==0]<-d1$DEATH[d1$PFS2==0]

p.value <- surv_pvalue(fit = survfit(Surv(PDAY2, PFS2) ~ ARM, data=d1), data = d1)
km<-survfit(Surv(PDAY2, PFS2) ~ ARM, data=d1)
summary(km,times = c(365))

ggsurv<-ggsurvplot(
  survfit(Surv(PDAY2, PFS2) ~ ARM, data=d1),
  conf.int = conf.int,
  pval = round(p.value$pval,3),
  risk.table = risk.table,
  ncensor.plot = ncensor.plot,
  palette = palette,
  risk.table.col = "strata",
  legend.labs =
    c("ARM A","ARM B"),
  ggtheme = ggtheme,
  title = "Progression-Free Survival 2",
  ylab = "Fraction Progression Free",
  xlab = "Time from Randomisation (Days)",
  font.title = c(font.title),
  font.x = c(font.lab),
  font.y = c(font.lab),
  break.time.by = 100,
  font.tickslab = c(font.tick),
  risk.table.height = risk.table.height,
  ncensor.plot.height = ncensor.plot.height
)
ggsurv$plot <- ggsurv$plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 70), "points"))
ggsurv$table <- ggsurv$table + 
  theme(
    plot.margin = unit(c(10, 10, 10, 5), "points"))
ggsurv$ncensor.plot <- ggsurv$ncensor.plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 95), "points"))
ggsurv

summary(coxph(Surv(PDAY,PFS)~(`BRAF VAF ctDNA result`),data=d1))
summary(coxph(Surv(PDAY2,PFS2)~(`BRAF VAF ctDNA result`),data=d1))
summary(coxph(Surv(OS,DEATH)~(`BRAF VAF ctDNA result`),data=d1))

d1$ctDNA[d1$`BRAF VAF ctDNA result`<=10]<-"<=10"
d1$ctDNA[d1$`BRAF VAF ctDNA result`>10]<-">10"
d1$ctDNA[d1$`BRAF VAF ctDNA result`>20]<-">20"

p.value <- surv_pvalue(fit = survfit(Surv(OS, DEATH) ~ ctDNA, data=d1), data = d1)
ggsurv<-ggsurvplot(
  survfit(Surv(OS, DEATH) ~ ctDNA, data=d1),
  conf.int = conf.int,
  pval = round(p.value$pval,3),
  risk.table = risk.table,
  ncensor.plot = ncensor.plot,
  palette = palette,
  risk.table.col = "strata",
  legend.labs =
    c("ctDNA<=10","ctDNA>10"),
  ggtheme = ggtheme,
  title = "Ovarall Survival",
  ylab = "Survival Fraction",
  xlab = "Time from Randomisation (Days)",
  font.title = c(font.title),
  font.x = c(font.lab),
  font.y = c(font.lab),
  break.time.by = 100,
  font.tickslab = c(font.tick),
  risk.table.height = risk.table.height,
  ncensor.plot.height = ncensor.plot.height
)
print(ggsurv)

d1$ARM

p.value <- surv_pvalue(fit = survfit(Surv(OS, DEATH) ~ ctDNA, data=d1), data = d1[d1$ARM=="01|ARM A",])
ggsurv<-ggsurvplot(
  survfit(Surv(OS, DEATH) ~ ctDNA, data=d1[d1$ARM=="01|ARM A",]),
  conf.int = conf.int,
  pval = round(p.value$pval,3),
  risk.table = risk.table,
  ncensor.plot = ncensor.plot,
  palette = palette,
  risk.table.col = "strata",
  legend.labs =
    c("ctDNA<=10","ctDNA>10"),
  ggtheme = ggtheme,
  title = "Ovarall Survival: ARM A",
  ylab = "Survival Fraction",
  xlab = "Time from Randomisation (Days)",
  font.title = c(font.title),
  font.x = c(font.lab),
  font.y = c(font.lab),
  break.time.by = 100,
  font.tickslab = c(font.tick),
  risk.table.height = risk.table.height,
  ncensor.plot.height = ncensor.plot.height
)
print(ggsurv)


p.value <- surv_pvalue(fit = survfit(Surv(OS, DEATH) ~ ctDNA, data=d1), data = d1[d1$ARM=="02|ARM B",])
ggsurv<-ggsurvplot(
  survfit(Surv(OS, DEATH) ~ ctDNA, data=d1[d1$ARM=="02|ARM B",]),
  conf.int = conf.int,
  pval = round(p.value$pval,3),
  risk.table = risk.table,
  ncensor.plot = ncensor.plot,
  palette = palette,
  risk.table.col = "strata",
  legend.labs =
    c("ctDNA<=10","ctDNA>10"),
  ggtheme = ggtheme,
  title = "Ovarall Survival: ARM B",
  ylab = "Survival Fraction",
  xlab = "Time from Randomisation (Days)",
  font.title = c(font.title),
  font.x = c(font.lab),
  font.y = c(font.lab),
  break.time.by = 100,
  font.tickslab = c(font.tick),
  risk.table.height = risk.table.height,
  ncensor.plot.height = ncensor.plot.height
)
print(ggsurv)

# Toxicity #####################################################################

d13<-read.csv("../Data/CSV_export/CAcTUS.Table 13.csv")
colnames(d13)<-d13[1,]
d13<-d13[-1,]

# we first one to remove dates of toxicities before treatment start
d13$`Event at study entry`
d13<-d13[d13$`Event at study entry`=="00|No",]
# we want to collect the worst grade seen for each patient and the date
d13[,c("PID","AE No","Max CTCAE grade")]

unique(d13$`Max CTCAE grade`)
d13<-d13[d13$`Max CTCAE grade`!="88|Not applicable",]
# Splitting the Max CTCAE grade into two columns if needed (optional)
d13 <- d13 %>%
  mutate(Max_CTCAE_grade_1 = as.numeric(sub("^(\\d+)\\|\\d+$", "\\1", `Max CTCAE grade`)),
         Max_CTCAE_grade_2 = as.numeric(sub("^\\d+\\|(\\d+)$", "\\1", `Max CTCAE grade`)))

# Find the highest grade per patient
highest_grade_per_patient <- d13 %>%
  group_by(PID) %>%
  summarise(Highest_Grade = max(Max_CTCAE_grade_1, na.rm = TRUE))

# merge with d1
colnames(d1)
tox.df<-merge(d1[,c("PID","ARM","BRAF VAF ctDNA result")],highest_grade_per_patient,
              by=c("PID"))
summary(factor(tox.df$Highest_Grade[tox.df$ARM=="01|ARM A"]))
summary(factor(tox.df$Highest_Grade[tox.df$ARM=="02|ARM B"]))

tox.df$Highest_Grade<-factor(tox.df$Highest_Grade)
ggplot(tox.df,aes(Highest_Grade,`BRAF VAF ctDNA result`,
                  col=Highest_Grade))+geom_jitter(height=0,width=0.3,
                                                  size=4,alpha=0.3)+
  theme_bw(base_size=14)+xlab("Max. CTCAE Grade")+
  ylab("Mutant VAF [%]")

tox.df$ARM[tox.df$ARM=="02|ARM B"]<-"ARM B"
tox.df$ARM[tox.df$ARM=="01|ARM A"]<-"ARM A"
ggplot(tox.df,aes(Highest_Grade,`BRAF VAF ctDNA result`,
                  col=Highest_Grade))+geom_jitter(height=0,width=0.3,
                                                  size=4,alpha=0.3)+
  theme_bw(base_size=14)+xlab("Max. CTCAE Grade")+
  ylab("Mutant VAF [%]")+facet_wrap(~ARM)

# lets look at cumulative incidence time from randomisation

# Splitting the Max CTCAE grade into two columns if needed (optional)
d13$`Start date`<-as.Date()
d13 <- d13 %>%
  mutate(Max_CTCAE_grade_num = as.numeric(sub("^(\\d+)\\|\\d+$", "\\1", `Max CTCAE grade`)),
         EventDate = as.Date(EventDate, format = "%Y-%m-%d")) # Adjust the format as per your actual data

# Finding the first occurrence of the worst grade for each patient
worst_grade_first_occurrence <- d13 %>%
  arrange(PID, desc(Max_CTCAE_grade_num), EventDate) %>%
  group_by(PID) %>%
  summarise(First_Worst_Grade = first(Max_CTCAE_grade_num),
            Date_of_First_Worst_Grade = first(EventDate),
            .groups = 'drop')