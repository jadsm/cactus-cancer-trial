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
d1
d1$`If no reason ctDNA BRAF VAF less 5`[d1$`Eligible for study`=="00|No"]
# so quite a few were rejected due to VAF levels being too low so <5%
d1$`If no reason ctDNA BRAF VAF less 5`[is.na(d1$`BRAF VAF ctDNA result`)]
# so teh two who did not have ctDNA the reason for not including was not ctDNA
d1$`failed Inc number`[is.na(d1$`BRAF VAF ctDNA result`)]

d1$`BRAF VAF ctDNA result`<-as.numeric(d1$`BRAF VAF ctDNA result`)
d1[,c("PID","BRAF VAF ctDNA result")]

d1$`Date result received by Investigator`[d1$`Date result received by Investigator`=="8888-88-88"]<-NA
d1$`Date result received by Investigator`<-as.Date(d1$`Date result received by Investigator`, "%d/%m/%Y")
d1 <- d1 %>% mutate(`Eligible for study` = case_when(
  `Eligible for study` == "00|No" ~ "No",
  `Eligible for study` == "01|Yes" ~ "Yes"))

pdf("Fig1A.pdf", width = 8, height = 6.6)
ggplot(d1,aes(`Date result received by Investigator`,`BRAF VAF ctDNA result`,
              col=`Eligible for study`))+geom_point(size=3)+
  theme_bw(base_size=14)+scale_y_log10()+ylab("BRAF VAF ctDNA [%]")+
  geom_hline(yintercept = 1.5,col=1,lty=2)+
  geom_hline(yintercept = 5,col=1,lty=3)+
  geom_vline(xintercept = as.Date("2020-07-02"), linetype = "dashed", color = "#008080", linewidth = 1) +
  scale_color_manual(values = c("Yes" = "blue", "No" = "red"))
dev.off()

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
d1<-d1[d1$`Eligible for study`=="Yes",]
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
          "Other",
          "Other if yes specify",
          "Screening ID")]

summary(d1$`BRAF VAF ctDNA result`[d1$`Treatment arm`=="01|ARM A"])
summary(d1$`BRAF VAF ctDNA result`[d1$`Treatment arm`=="02|ARM B"])

d1[d1$`Treatment arm`=="01|ARM A",c("Liver",
                                    "Lung",
                                    "Subcutaneous",
                                    "Brain",
                                    "Nodal",
                                    "Abdominal",
                                    "Bone",
                                    "Other",
                                    "Other if yes specify")]

summary(c(4,6,5,2,5,2,3,3,5,2))


d1[d1$`Treatment arm`=="02|ARM B",c("Liver",
                                    "Lung",
                                    "Subcutaneous",
                                    "Brain",
                                    "Nodal",
                                    "Abdominal",
                                    "Bone",
                                    "Other",
                                    "Other if yes specify")]

summary(c(5,9,4,3,3,2,3,3,5,2,3))
# we will keep these in case they start asking for adjusted analyses
d1$`BRAF VAF ctDNA result`[d1$`Treatment arm`=="02|ARM B"]
d1$`BRAF VAF ctDNA result`[d1$`Treatment arm`=="01|ARM A"]



# KM plots ----------------------------------------------------------------

# we've been told by email that 3 patients had been defined as clinically
# progressing at a certain dates:
# 04-009 10/11/21
# 04-014 21/2/22 
# 01-004 23/12/19
# we will alter these manually when we go through the data



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

# Arm A #FFA741
# Arm B #3684BC

# #087383 Targeted
# #F37B6A IO


# CR #34C59A
# PR #FFA741
# SD #81B4DA
# NE #7F93A6
# set options for plotting survival curves
font.lab     <- 14
font.title   <- 16
font.tick    <- 12
conf.int     <- T
risk.table   <- T
ncensor.plot <- T
palette      <- c('#FFA741','#3684BC')
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
    plot.margin = unit(c(10, 10, 10, 64), "points"))
ggsurv$ncensor.plot <- ggsurv$ncensor.plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 85), "points"))

survfit(Surv(OS, DEATH) ~ ARM, data=d1)


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
    plot.margin = unit(c(10, 10, 10, 64), "points"))
ggsurv$ncensor.plot <- ggsurv$ncensor.plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 85), "points"))
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

# it hear we introduce the progression dates
# 04-009 10/11/21
# 04-014 21/2/22 
# 01-004 23/12/19
d6<-rbind(d6,c(PID="9",Date="10/11/2021",`Overall response`="04|PD"))
d6<-rbind(d6,c(PID="14",Date="21/02/2022",`Overall response`="04|PD"))
d6<-rbind(d6,c(PID="4",Date="23/12/2019",`Overall response`="04|PD"))


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
    plot.margin = unit(c(10, 10, 10, 64), "points"))
ggsurv$ncensor.plot <- ggsurv$ncensor.plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 65), "points"))
ggsurv

# they did ask for PFS rate at 1 year too
km<-survfit(Surv(PDAY, PFS) ~ ARM, data=d1)
summary(km,times = c(365))

survfit(formula = Surv(PDAY, PFS) ~ ARM, data = d1)

# ARM=01|ARM A 
# time       n.risk      n.event     survival      std.err lower 95% CI upper 95% CI 
# 365.000        2.000        4.000        0.560        0.171        0.308        1.000 
# 
# ARM=02|ARM B 
# time       n.risk      n.event     survival      std.err lower 95% CI upper 95% CI 
# 365.000        3.000        4.000        0.606        0.154        0.368        0.998 

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
# ARM=01|ARM A 
# time       n.risk      n.event     survival      std.err lower 95% CI upper 95% CI 
# 365.000        3.000        4.000        0.600        0.155        0.362        0.995 
# 
# ARM=02|ARM B 
# time       n.risk      n.event     survival      std.err lower 95% CI upper 95% CI 
# 365.000        3.000        5.000        0.420        0.176        0.184        0.956
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
    plot.margin = unit(c(10, 10, 10, 64), "points"))
ggsurv$ncensor.plot <- ggsurv$ncensor.plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 65), "points"))
ggsurv

# summary(coxph(Surv(PDAY,PFS)~(`BRAF VAF ctDNA result`),data=d1))
# summary(coxph(Surv(PDAY2,PFS2)~(`BRAF VAF ctDNA result`),data=d1))
# summary(coxph(Surv(OS,DEATH)~(`BRAF VAF ctDNA result`),data=d1))
# 
# d1$ctDNA[d1$`BRAF VAF ctDNA result`<=10]<-"<=10"
# d1$ctDNA[d1$`BRAF VAF ctDNA result`>10]<-">10"
# d1$ctDNA[d1$`BRAF VAF ctDNA result`>20]<-">20"
# 
# p.value <- surv_pvalue(fit = survfit(Surv(OS, DEATH) ~ ctDNA, data=d1), data = d1)
# ggsurv<-ggsurvplot(
#   survfit(Surv(OS, DEATH) ~ ctDNA, data=d1),
#   conf.int = conf.int,
#   pval = round(p.value$pval,3),
#   risk.table = risk.table,
#   ncensor.plot = ncensor.plot,
#   palette = palette,
#   risk.table.col = "strata",
#   legend.labs =
#     c("ctDNA<=10","ctDNA>10"),
#   ggtheme = ggtheme,
#   title = "Ovarall Survival",
#   ylab = "Survival Fraction",
#   xlab = "Time from Randomisation (Days)",
#   font.title = c(font.title),
#   font.x = c(font.lab),
#   font.y = c(font.lab),
#   break.time.by = 100,
#   font.tickslab = c(font.tick),
#   risk.table.height = risk.table.height,
#   ncensor.plot.height = ncensor.plot.height
# )
# print(ggsurv)
# 
# d1$ARM
# 
# p.value <- surv_pvalue(fit = survfit(Surv(OS, DEATH) ~ ctDNA, data=d1), data = d1[d1$ARM=="01|ARM A",])
# ggsurv<-ggsurvplot(
#   survfit(Surv(OS, DEATH) ~ ctDNA, data=d1[d1$ARM=="01|ARM A",]),
#   conf.int = conf.int,
#   pval = round(p.value$pval,3),
#   risk.table = risk.table,
#   ncensor.plot = ncensor.plot,
#   palette = palette,
#   risk.table.col = "strata",
#   legend.labs =
#     c("ctDNA<=10","ctDNA>10"),
#   ggtheme = ggtheme,
#   title = "Ovarall Survival: ARM A",
#   ylab = "Survival Fraction",
#   xlab = "Time from Randomisation (Days)",
#   font.title = c(font.title),
#   font.x = c(font.lab),
#   font.y = c(font.lab),
#   break.time.by = 100,
#   font.tickslab = c(font.tick),
#   risk.table.height = risk.table.height,
#   ncensor.plot.height = ncensor.plot.height
# )
# print(ggsurv)
# 
# 
# p.value <- surv_pvalue(fit = survfit(Surv(OS, DEATH) ~ ctDNA, data=d1), data = d1[d1$ARM=="02|ARM B",])
# ggsurv<-ggsurvplot(
#   survfit(Surv(OS, DEATH) ~ ctDNA, data=d1[d1$ARM=="02|ARM B",]),
#   conf.int = conf.int,
#   pval = round(p.value$pval,3),
#   risk.table = risk.table,
#   ncensor.plot = ncensor.plot,
#   palette = palette,
#   risk.table.col = "strata",
#   legend.labs =
#     c("ctDNA<=10","ctDNA>10"),
#   ggtheme = ggtheme,
#   title = "Ovarall Survival: ARM B",
#   ylab = "Survival Fraction",
#   xlab = "Time from Randomisation (Days)",
#   font.title = c(font.title),
#   font.x = c(font.lab),
#   font.y = c(font.lab),
#   break.time.by = 100,
#   font.tickslab = c(font.tick),
#   risk.table.height = risk.table.height,
#   ncensor.plot.height = ncensor.plot.height
# )
# print(ggsurv)


# waterfall plots of imaging ----------------------------------------------
# Arm A v #FFA741
# Arm B #3684BC

# #087383 Targeted
# #F37B6A IO


# CR #34C59A
# PR #FFA741
# SD #81B4DA
# NE #7F93A6
# we will 
d6<-read.csv("../Data/CSV_export/CAcTUS.Table 6_RECIST.csv")

colnames(d6)<-d6[1,]
d6<-d6[-1,]
colnames(d6)
d6$`CT scan date`[d6$`CT scan date`=="8888-88-88"]<-NA
d6$`PET CT scan`[d6$`PET CT scan`=="8888-88-88"]<-NA

d6$Date<-d6$`CT scan date`
# when these are missing we use PET CT
d6$Date[is.na(d6$Date)]<-d6$`PET CT scan`[is.na(d6$Date)]

# we need to set some valeus to NA  basically the 88 code
d6$`TL1 mm`[d6$`TL1 mm`=="888"]<-NA
d6$`TL1 mm`[d6$`TL1 mm`=="991"]<-NA
d6$`TL2 mm`[d6$`TL2 mm`=="888"]<-NA
d6$`TL2 mm`[d6$`TL2 mm`=="991"]<-NA
d6$`TL3 mm`[d6$`TL3 mm`=="888"]<-NA
d6$`TL3 mm`[d6$`TL3 mm`=="991"]<-NA
d6$`TL4 mm`[d6$`TL4 mm`=="888"]<-NA
d6$`TL4 mm`[d6$`TL4 mm`=="991"]<-NA
d6$`TL5 mm`[d6$`TL5 mm`=="888"]<-NA
d6$`TL5 mm`[d6$`TL5 mm`=="991"]<-NA

colnames(d6)
d6<-d6[,c("PID","Schedule","TL1 mm","TL2 mm","TL3 mm","TL4 mm","TL5 mm","Overall response","Treatment","Date")]%>%
  filter(`Overall response`!="05|NE")
# in case you want to exclude the 'Not evaluable label'
d6<-merge(d6,d1[,c("PID","Treatment arm","Dab and Tram immediately","First_PD_or_Last_Date","Second_PD_Date")],by=c("PID"),all.x=T,all.y=T)
d6$SUM<-NA
d6$`TL1 mm`<-as.numeric(d6$`TL1 mm`)
d6$`TL2 mm`<-as.numeric(d6$`TL2 mm`)
d6$`TL3 mm`<-as.numeric(d6$`TL3 mm`)
d6$`TL4 mm`<-as.numeric(d6$`TL4 mm`)
d6$`TL5 mm`<-as.numeric(d6$`TL5 mm`)
d6$SUM<-rowSums(d6[,c("TL1 mm","TL2 mm","TL3 mm","TL4 mm","TL5 mm")],na.rm=T)


dummy<-d6[d6$Schedule=="00|Baseline",c("PID","SUM","Treatment arm")]

colnames(dummy)[2]<-"BSL"
summary(dummy$BSL[dummy$`Treatment arm`=="01|ARM A"])
summary(dummy$BSL[dummy$`Treatment arm`=="02|ARM B"])

d6<-merge(d6,dummy[,c("PID","BSL")],by=c("PID"),all.x=T,all.y=T)
d6$CHG<-d6$SUM-d6$BSL
d6$PCHG<-100*d6$CHG/d6$BSL
unique(d6$Treatment)
d6$Treatment[d6$Treatment=="88|Not applicable"]<-"00|Baseline"
sum(d6$Treatment=="03|Mono")


# -	Arm A vs B % change first line treatment (for arm B this would be considered 
# the targeted run in + the immune therapy) to 1st progression

d6$Date<-as.Date(d6$Date,"%d/%m/%Y")
d6$First_PD_or_Last_Date<-as.Date(d6$First_PD_or_Last_Date,"%Y-%m-%d")
d6$Second_PD_Date<-as.Date(d6$Second_PD_Date,"%Y-%m-%d")

d6$BPCH1<-NA
id<-unique(d6$PID)
for (ii in 1:length(id)){
  dummy<-d6$PCHG[d6$PID==id[ii] & d6$Schedule!="00|Baseline" & d6$Date<=d6$First_PD_or_Last_Date]
  d6$BPCH1[d6$PID==id[ii]]<-min(dummy)
  # if (sum(dummy$`Treatment arm`=="01|ARM A" & c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="03|Mono"))>0){
  #   d6$BPCH[d6$PID==id[ii]]<-min(dummy$PCH[c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="03|Mono")])
  # }
  # if(sum(dummy$`Treatment arm`=="02|ARM B" & c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="02|Immune"))>0){
  #   d6$BPCH[d6$PID==id[ii]]<-min(dummy$PCH[c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="02|Immune")])
  # }
}

df1<-unique(d6[,c("PID","BPCH1","Treatment arm")])
b1<-df1[df1$`Treatment arm`=="01|ARM A",]
b1<-arrange(b1,BPCH1)
b2<-df1[df1$`Treatment arm`=="02|ARM B",]
b2<-arrange(b2,BPCH1)

# now get the colours by recist best - # red - PD, Orange - SD, blue - PR, green - CR
mapping <- c(
  "CR"="#34C59A","PR"="#FFA741","SD"="#81B4DA","PD"="red","NE"='gray'
)

df_separated <- d6 %>%
  # filter(`Overall response`!="5|NE") %>%
  separate(`Overall response`, into = c("response_code", "response"), sep = "\\|") %>%
  filter(BPCH1 == PCHG) %>%
  # arrange(PID, response_code) %>% # Sort: id ascending, timestamp descending
  distinct(PID, .keep_all = TRUE) %>%
  select(PID, "response","Treatment arm")
colours1<-df_separated[df_separated$`Treatment arm`=="01|ARM A",]
colours2<-df_separated[df_separated$`Treatment arm`=="02|ARM B",]

# subFig2.2. S
par(mar=c(5,5,5,5))
sorted_PIDs <- b1$PID[order(b1$BPCH1, decreasing = TRUE)]
colours1<-colours1 %>%
  mutate(PID = factor(PID, levels = sorted_PIDs))%>%
  arrange(PID) %>%
  mutate(PID = as.character(PID))
colours1<- mapping[colours1$response]
pdf("Subfig22S.pdf", width = 8, height = 6)
bar_positions <- barplot(sort(b1$BPCH1, decreasing = TRUE), 
                         col = colours1, #"#FFA741"
                         main = "ARM A", 
                         cex.axis = 1.5, 
                         ylab = "Best % Change SLD", 
                         cex.lab = 1.5,ylim=c(-113,100))  # Increases the y-axis label size
text(x = bar_positions, 
     y = sort(b1$BPCH1, decreasing = TRUE), 
     labels = sorted_PIDs, 
     pos = 1,   # Position labels above the bars
     cex = 1.2, # Adjust the size of the labels
     col = "black")  # Color of the labels
dev.off()

# subFig2.2. T
sorted_PIDs <- b2$PID[order(b2$BPCH1[is.finite(b2$BPCH1)], decreasing = TRUE)]
colours2<-colours2 %>%
  mutate(PID = factor(PID, levels = sorted_PIDs))%>%
  arrange(PID) %>%
  mutate(PID = as.character(PID))
colours2<- mapping[colours2$response]
pdf("Subfig22T.pdf", width = 8, height = 6)
bar_positions <- barplot(sort(b2$BPCH1[is.finite(b2$BPCH1)], decreasing = TRUE), 
                         col = colours2,#"#3684BC" 
                         main = "ARM B", 
                         cex.axis = 1.5, 
                         ylab = "Best % Change SLD", 
                         cex.lab = 1.5,ylim=c(-113,100))  # Increases the y-axis label size
text(x = bar_positions, 
     y = sort(b2$BPCH1[is.finite(b2$BPCH1)], decreasing = TRUE), 
     labels = sorted_PIDs, 
     pos = 1,   # Position labels above the bars
     cex = 1.2, # Adjust the size of the labels
     col = "black")  # Color of the labels
dev.off()

# export the lengend only
pdf('legend2_2.pdf', width = 3, height = 2)
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
plot.new() # Creates a new, empty plot. Simpler than plot(1, type="n",...) for this purpose.
legend("center", # Position the legend (e.g., "center", "topleft", "right")
       legend = names(mapping), # Labels for legend items (CR, PR, SD, PD)
       fill = mapping,         # Colors for the legend keys
       bty = "n",                      # No box around the legend
       cex = 1.2,                      # Character expansion factor (text size)
       title = "Response Status"       # Optional title for the legend
)

# 5. Close the graphics device to save the file
dev.off()


# 43.396226  -8.474576 -17.543860 -26.562500 -41.176471 -44.444444 -47.058824 -60.759494 -63.235294

# -	Arm A vs B % change second line treatment
d6$BPCH2<-NA
id<-unique(d6$PID)
for (ii in 1:length(id)){
  dummy<-d6$PCHG[d6$PID==id[ii] & d6$Schedule!="00|Baseline" & d6$Date<=d6$Second_PD_Date & d6$Date>d6$First_PD_or_Last_Date]
  d6$BPCH2[d6$PID==id[ii]]<-min(dummy)
  # if (sum(dummy$`Treatment arm`=="01|ARM A" & c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="03|Mono"))>0){
  #   d6$BPCH[d6$PID==id[ii]]<-min(dummy$PCH[c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="03|Mono")])
  # }
  # if(sum(dummy$`Treatment arm`=="02|ARM B" & c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="02|Immune"))>0){
  #   d6$BPCH[d6$PID==id[ii]]<-min(dummy$PCH[c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="02|Immune")])
  # }
}

df1<-unique(d6[,c("PID","BPCH1","BPCH2","Treatment arm")])
b1<-df1[df1$`Treatment arm`=="01|ARM A",]
b1<-arrange(b1,BPCH2)
b2<-df1[df1$`Treatment arm`=="02|ARM B",]
b2<-arrange(b2,BPCH2)
par(mar=c(5,5,5,5))
bar_positions <- barplot(sort(b1$BPCH2[is.finite(b1$BPCH2)], decreasing = TRUE), 
                         col = "#FFA741", 
                         main = "ARM A", 
                         cex.axis = 1.5, 
                         ylab = "Best % Change SLD", 
                         cex.lab = 1.5,ylim=c(-113,100))  # Increases the y-axis label size
sorted_PIDs <- b1$PID[order(b1$BPCH2[is.finite(b1$BPCH2)], decreasing = TRUE)]
text(x = bar_positions, 
     y = sort(b1$BPCH2[is.finite(b1$BPCH2)], decreasing = TRUE), 
     labels = sorted_PIDs, 
     pos = 1,   # Position labels above the bars
     cex = 1.2, # Adjust the size of the labels
     col = "black")  # Color of the labels


pdf('Fig2K.pdf', width = 8, height = 6)
bar_positions <- barplot(sort(b2$BPCH2[is.finite(b2$BPCH2)], decreasing = TRUE), 
                         col = "#3684BC", 
                         main = "", 
                         cex.axis = 1.5, 
                         ylab = "Best % Change SLD", 
                         cex.lab = 1.5,ylim=c(-113,0))  # Increases the y-axis label size
sorted_PIDs <- b2$PID[order(b2$BPCH2[is.finite(b2$BPCH2)], decreasing = TRUE)]
text(x = bar_positions, 
     y = sort(b2$BPCH2[is.finite(b2$BPCH2)], decreasing = TRUE), 
     labels = sorted_PIDs, 
     pos = 1,   # Position labels above the bars
     cex = 1.2, # Adjust the size of the labels
     col = "black")  # Color of the labels
dev.off()

pdf('Fig2L.pdf', width = 8, height = 6)
bar_positions <- barplot(c(-53,-99,-99,-99,-99), 
                         col = "#008080", 
                         main = "", 
                         cex.axis = 1.5, 
                         ylab = "Best % Change Mutant VAF", 
                         cex.lab = 1.5,ylim=c(-113,0))  # Increases the y-axis label size
text(x = bar_positions, 
     y = c(-53,-99,-99,-99,-99),
     labels = sorted_PIDs, 
     pos = 1,   # Position labels above the bars
     cex = 1.2, # Adjust the size of the labels
     col = "black")  # Color of the labels
dev.off()

# now we manually do VAf best % change - read off graphs
# Supp Fig 2.2U
pdf("Subfig22Ualt.pdf", width = 8, height = 6)
bar_positions <- barplot(c(-96,-99,-99,-99,-99,-99,-99,-99,-99,-99), 
                         col = colours1, 
                         main = "ARM A", 
                         cex.axis = 1.5, 
                         ylab = "Best % Change VAF", 
                         cex.lab = 1.5,ylim=c(-113,100))  # Increases the y-axis label size
sorted_PIDs <- c(18,3,9,8,1,17,16,15,6,10)
text(x = bar_positions, 
     y = c(-96,-99,-99,-99,-99,-99,-99,-99,-99,-99), 
     labels = sorted_PIDs, 
     pos = 1,   # Position labels above the bars
     cex = 1.2, # Adjust the size of the labels
     col = "black")  # Color of the labels
dev.off()

# Supp Fig 2.2V
pdf("Subfig22V.pdf", width = 8, height = 6)
bar_positions <- barplot(c(-99,-85,-99,-99,-99,-99,-99,-99,-99,-99), 
                         col = colours2, 
                         main = "ARM B", 
                         cex.axis = 1.5, 
                         ylab = "Best % Change VAF", 
                         cex.lab = 1.5,ylim=c(-113,100))  # Increases the y-axis label size
sorted_PIDs <- c(11,13,4,21,5,19,20,7,12,14)
text(x = bar_positions, 
     y = c(-99,-85,-99,-99,-99,-99,-99,-99,-99,-99), 
     labels = sorted_PIDs, 
     pos = 1,   # Position labels above the bars
     cex = 1.2, # Adjust the size of the labels
     col = "black")  # Color of the labels
dev.off()

# Supp Fig 2.2W - read from graph
pdf("Subfig22W.pdf", width = 8, height = 6)
data22W<-c(-99,-100*(10-2.4)/10,-99,-99,0,-99,-99,-99,-98.3,-99)
bar_positions <- barplot(data22W, 
                         col = colours2, 
                         main = "ARM B", 
                          cex.axis = 1.5, 
                         ylab = "% Change VAF from baseline to the switch", 
                         cex.lab = 1.5,ylim=c(-113,100))  # Increases the y-axis label size
sorted_PIDs <- c(11,13,4,21,5,19,20,7,12,14)
text(x = bar_positions, 
     y = data22W, 
     labels = sorted_PIDs, 
     pos = 1,   # Position labels above the bars
     cex = 1.2, # Adjust the size of the labels
     col = "black")  # Color of the labels
dev.off()

# Supp Fig 2.2X - read from graph
pdf("Subfig22X.pdf", width = 8, height = 6)
data22X<-c(-100*(24.5-7.5)/24.5,-100*(10-2.6)/10,-99,-100*(10-15)/10,0,-99,-99,-99,-99,-99)
bar_positions <- barplot(data22X, 
                         col = colours2, 
                         main = "ARM B", 
                         cex.axis = 1.5, 
                         ylab = "% Change VAF to end of trial after CPI", 
                         cex.lab = 1.5,ylim=c(-113,100))  # Increases the y-axis label size
sorted_PIDs <- c(11,13,4,21,5,19,20,7,12,14)
text(x = bar_positions, 
     y = data22X, 
     labels = sorted_PIDs, 
     pos = 1,   # Position labels above the bars
     cex = 1.2, # Adjust the size of the labels
     col = "black")  # Color of the labels
dev.off()



# -	Arm B the % change from baseline scan to the scan which was performed at switch to immune therapy
# need to get data on 1st period of D+T treatment
d6<-arrange(d6,`Treatment arm`,PID, Date)
d6$BPCH3<-NA
id<-unique(d6$PID[d6$`Treatment arm`=="02|ARM B"])
for (ii in 1:length(id)){
  dummy<-d6$PCHG[d6$PID==id[ii] & d6$Schedule!="00|Baseline" & d6$Treatment=="01|D + T 1"]
  d6$BPCH3[d6$PID==id[ii]]<-min(dummy)
  # if (sum(dummy$`Treatment arm`=="01|ARM A" & c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="03|Mono"))>0){
  #   d6$BPCH[d6$PID==id[ii]]<-min(dummy$PCH[c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="03|Mono")])
  # }
  # if(sum(dummy$`Treatment arm`=="02|ARM B" & c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="02|Immune"))>0){
  #   d6$BPCH[d6$PID==id[ii]]<-min(dummy$PCH[c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="02|Immune")])
  # }
}
df1<-unique(d6[,c("PID","BPCH3","Treatment arm")])
b2<-df1[df1$`Treatment arm`=="02|ARM B",]
b2<-arrange(b2,BPCH3)
barplot(b2$BPCH3[is.finite(b2$BPCH3)],ylim=c(-100,25),
        col="#008000",main="ARM B",cex.axis=1.5)

pdf("fig2I.pdf", width = 8, height = 6)
bar_positions <- barplot(sort(b2$BPCH3[is.finite(b2$BPCH3)], decreasing = TRUE), 
                         col = "#3684BC", 
                         main = "", 
                         cex.axis = 1.5, 
                         ylab = "Best % Change SLD during TT run-in", 
                         cex.lab = 1.5,ylim=c(-80,60),yaxt = "n")  # Increases the y-axis label size
axis(2, at = seq(-80, 60, by = 20), cex.axis = 1.5)
sorted_PIDs <- b2$PID[order(b2$BPCH3[is.finite(b2$BPCH3)], decreasing = TRUE)]
text(x = bar_positions, 
     y = sort(b2$BPCH3[is.finite(b2$BPCH3)], decreasing = TRUE), 
     labels = sorted_PIDs, 
     pos = 1,   # Position labels above the bars
     cex = 1.2, # Adjust the size of the labels
     col = "black")  # Color of the labels
dev.off()
# patient 2 did not have post imaging and patient 11 switched after 2 weeks

# -	Arm B the best percentage change from the scan at switch up to either study end 
# or progression (this is to give an idea of the response to the IO component of the Arm B strategy)

d6$BPCH4<-NA
id<-unique(d6$PID[d6$`Treatment arm`=="02|ARM B"])
for (ii in 1:length(id)){
  dummy<-d6$PCHG[d6$PID==id[ii] & d6$Schedule!="00|Baseline" & d6$Treatment=="02|Immune"]
  d6$BPCH4[d6$PID==id[ii]]<-min(dummy)
  # if (sum(dummy$`Treatment arm`=="01|ARM A" & c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="03|Mono"))>0){
  #   d6$BPCH[d6$PID==id[ii]]<-min(dummy$PCH[c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="03|Mono")])
  # }
  # if(sum(dummy$`Treatment arm`=="02|ARM B" & c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="02|Immune"))>0){
  #   d6$BPCH[d6$PID==id[ii]]<-min(dummy$PCH[c(dummy$Treatment=="01|D + T 1"|dummy$Treatment=="02|Immune")])
  # }
}
df1<-unique(d6[,c("PID","BPCH4","Treatment arm")])
b2<-df1[df1$`Treatment arm`=="02|ARM B",]
b2<-arrange(b2,BPCH4)

pdf("fig2J.pdf", width = 8, height = 6)
bar_positions <- barplot(sort(b2$BPCH4[is.finite(b2$BPCH4)], decreasing = TRUE), 
                         col = "#3684BC", 
                         main = "", 
                         cex.axis = 1.5, 
                         ylab = "Best % Change SLD during CPI", 
                         cex.lab = 1.5,# Increases the y-axis label size
                         ylim=c(-80,60),
                         yaxt = "n")  
axis(2, at = seq(-80, 60, by = 20), cex.axis = 1.5)
sorted_PIDs <- b2$PID[order(b2$BPCH4[is.finite(b2$BPCH4)], decreasing = TRUE)]
text(x = bar_positions, 
     y = sort(b2$BPCH4[is.finite(b2$BPCH4)], decreasing = TRUE), 
     labels = sorted_PIDs, 
     pos = 1,   # Position labels above the bars
     cex = 1.2, # Adjust the size of the labels
     col = "black")  # Color of the labels
dev.off()

df1<-unique(d6[,c("PID","BPCH3","BPCH4","Treatment arm")])
plot(df1$BPCH3,df1$BPCH4,xlab="Targeted",ylab="IO")
abline(a=0,b=1)

# Toxicity #####################################################################

d13<-read.csv("../Data/CSV_export/CAcTUS.Table 13.csv")
colnames(d13)<-d13[1,]
d13<-d13[-1,]

# we do now have groupings of the toxicities
d13$`AE Code`
# we will ignore all 7777s
d13<-d13[d13$`AE Code`!="7777",]

# what we want are those that are related to drug...
unique(d13$`Dab study therapy`)
unique(d13$`Tram study therapy`)
unique(d13$`Nivo study therapy`)
unique(d13$`Iplili study therapy`)
# we don't want 00, 02, or 03

# Define the values to filter out
values_to_remove <- c("88|Not applicable", "00|No", "02|Not Related", "03|Unlikely")

# Filter the dataset
d13 <- d13[!(d13$`Dab study therapy` %in% values_to_remove &
                        d13$`Tram study therapy` %in% values_to_remove &
                        d13$`Related to study pro` %in% values_to_remove &
                        d13$`Nivo study therapy` %in% values_to_remove &
                        d13$`Iplili study therapy` %in% values_to_remove), ]

unique(d13$`AE Code`)
d13$`AE Code`<-as.numeric(d13$`AE Code`)
# we want to collect the worst grade seen for each patient and the date

unique(d13$`Max CTCAE grade`)
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

# we want a summary of toxiicty by each arm

d13<-merge(d13,unique(d1[,c("PID","ARM")]),by=c("PID"),all.x=T,all.y=F)

unique(d13$`AE Code`)
d13$`Max CTCAE grade`
colnames(d13)
d13[,c("PID","AE Code","ARM","Max CTCAE grade","Dab study therapy",
       "Tram study therapy","Nivo study therapy","Iplili study therapy")]

summary <- d13 %>%
  group_by(PID, `AE Code`, ARM) %>%
  summarize(`Max CTCAE grade` = first(`Max CTCAE grade`)) %>%
  ungroup() %>%
  count(`AE Code`, ARM, `Max CTCAE grade`) %>%
  arrange(`AE Code`, ARM, `Max CTCAE grade`)


# Spread the data for better readability
summary_spread <- summary %>%
  pivot_wider(names_from = ARM, values_from = n, values_fill = 0)
summary_spread$`Max CTCAE grade`[summary_spread$`Max CTCAE grade`=="01|01"]<-1
summary_spread$`Max CTCAE grade`[summary_spread$`Max CTCAE grade`=="02|02"]<-2
summary_spread$`Max CTCAE grade`[summary_spread$`Max CTCAE grade`=="03|03"]<-3
summary_spread$`Max CTCAE grade`[summary_spread$`Max CTCAE grade`=="04|04"]<-4

summary_spread<-arrange(summary_spread,`AE Code`,`Max CTCAE grade`)

write.csv(summary_spread,"toxicity.csv",row.names = F)

values_to_keep <- c("04|Possible", "01|Yes", "05|Probable", "06|Causal Relationship")
# Filter the dataset
dummy <- d13[(d13$`Dab study therapy` %in% values_to_keep &
               d13$`Tram study therapy` %in% values_to_keep), ]


summary <- dummy %>%
  group_by(PID, `AE Code`, ARM) %>%
  summarize(`Max CTCAE grade` = first(`Max CTCAE grade`)) %>%
  ungroup() %>%
  count(`AE Code`, ARM, `Max CTCAE grade`) %>%
  arrange(`AE Code`, ARM, `Max CTCAE grade`)


# Spread the data for better readability
summary_spread <- summary %>%
  pivot_wider(names_from = ARM, values_from = n, values_fill = 0)
summary_spread$`Max CTCAE grade`[summary_spread$`Max CTCAE grade`=="01|01"]<-1
summary_spread$`Max CTCAE grade`[summary_spread$`Max CTCAE grade`=="02|02"]<-2
summary_spread$`Max CTCAE grade`[summary_spread$`Max CTCAE grade`=="03|03"]<-3
summary_spread$`Max CTCAE grade`[summary_spread$`Max CTCAE grade`=="04|04"]<-4

summary_spread<-arrange(summary_spread,`AE Code`,`Max CTCAE grade`)

write.csv(summary_spread,"toxicity_targeted.csv",row.names = F)




# Filter the dataset
dummy <- d13[(d13$`Nivo study therapy` %in% values_to_keep &
                d13$`Iplili study therapy` %in% values_to_keep), ]


summary <- dummy %>%
  group_by(PID, `AE Code`, ARM) %>%
  summarize(`Max CTCAE grade` = first(`Max CTCAE grade`)) %>%
  ungroup() %>%
  count(`AE Code`, ARM, `Max CTCAE grade`) %>%
  arrange(`AE Code`, ARM, `Max CTCAE grade`)


# Spread the data for better readability
summary_spread <- summary %>%
  pivot_wider(names_from = ARM, values_from = n, values_fill = 0)
summary_spread$`Max CTCAE grade`[summary_spread$`Max CTCAE grade`=="01|01"]<-1
summary_spread$`Max CTCAE grade`[summary_spread$`Max CTCAE grade`=="02|02"]<-2
summary_spread$`Max CTCAE grade`[summary_spread$`Max CTCAE grade`=="03|03"]<-3
summary_spread$`Max CTCAE grade`[summary_spread$`Max CTCAE grade`=="04|04"]<-4

summary_spread<-arrange(summary_spread,`AE Code`,`Max CTCAE grade`)

write.csv(summary_spread,"toxicity_IO.csv",row.names = F)
