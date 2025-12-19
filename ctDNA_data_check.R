# Here we try and find all ctDNA data and compare with...
require(tidyverse)
require(readxl)



# load and clean ctDNA data -----------------------------------------------


orig.df<-read_xlsx("../Data/CAcTUS VAF final report ddPCR Data for sponsor_R1.xlsx",sheet=1)
# note the file contains replicates but these are only for screening samples
# after discussions with Becki it was decided we will use the excel sheet for the bulk of the
# analyses


orig.df$Patient<-as.numeric(str_remove(orig.df$Patient,"-"))

# the IDs are a mix of screening and patient ID, the link between the two is taken from the
# ctDNA report

links<-read_xlsx("../Data/ID_links.xlsx")
links$`Patient ID`<-as.numeric(str_remove(links$`Patient ID`,"-"))
links$`Screening ID`<-as.numeric(str_remove(links$`Screening ID`,"-"))

orig.df$Screening<-"No"
colnames(orig.df)[1]<-"ID"

orig.df$Screening[orig.df$ID %in% links$`Screening ID`]<-"Yes"


# lets now assign red and green labelling

red.df<-read_xlsx("../Data/CAcTUS VAF final report ddPCR Data for sponsor_R1.xlsx",sheet=2)
red.df<-red.df[is.na(red.df$`Patient ID`)==F,]
colnames(red.df)
red.df<-red.df[,c("Patient ID","Visit","Working days from Receipt to report")]
red.df$`Patient ID`<-as.numeric(str_remove(red.df$`Patient ID`,"-"))
colnames(red.df)[1]<-"ID"
orig.df<-merge(orig.df,red.df,by=c("ID","Visit"),all.x=T,all.y=T)
orig.df$Red<-"No"
orig.df$Red[is.na(orig.df$`Working days from Receipt to report`)==F]<-"Yes"

green.df<-read_xlsx("../Data/CAcTUS VAF final report ddPCR Data for sponsor_R1.xlsx",sheet=3)
green.df<-green.df[is.na(green.df$`Patient ID`)==F,]
green.df$`Patient ID`<-as.numeric(str_remove(green.df$`Patient ID`,"-"))
green.df$Green<-"Yes"
green.df<-green.df[,c("Patient ID","Date of Visit","Green")]
colnames(green.df)[1]<-"ID"
colnames(green.df)[2]<-"Sample Date"
orig.df<-merge(orig.df,green.df,by=c("ID","Sample Date"),all.x=T,all.y=T)

orig.df$Green[is.na(orig.df$Green)]<-"No"

# one screening ID was labelled baseline but we shall rename as screening 04-605
orig.df$Visit[orig.df$ID==4605]<-"Screening"
sum(orig.df$Visit=="Screening")
sum(orig.df$`Percentage Mutant VAF`[orig.df$Visit=="Screening"]<orig.df$`%VAF LLOQ`[orig.df$Visit=="Screening"])
sum(orig.df$`Percentage Mutant VAF`[orig.df$Visit=="Screening"]<orig.df$`%VAF LLOD`[orig.df$Visit=="Screening"])


length(unique(orig.df$ID[orig.df$`Percentage Mutant VAF`<orig.df$`%VAF LLOQ` & orig.df$Visit=="Screening"]))
length(unique(orig.df$ID[orig.df$`Percentage Mutant VAF`<orig.df$`%VAF LLOD` & orig.df$Visit=="Screening"]))

# collect all screening visit and describe concordance
dummy<-orig.df[orig.df$Visit=="Screening",c("ID","Percentage Mutant VAF")]
dummy

library(epiR)
library(dplyr)

dummy_wide <- dummy %>%
  group_by(ID) %>%
  summarise(Sample1 = first(`Percentage Mutant VAF`), Sample2 = last(`Percentage Mutant VAF`))
ccc_result <- epi.ccc(dummy_wide$Sample1, dummy_wide$Sample2)

print(ccc_result)
cor.test(dummy_wide$Sample1,dummy_wide$Sample2)

# Now we have all the red and green samples assigned and the screening values
# assigend properly, we will now assign the PID so that we can link up to the
# clinical data

# If we load up Table 1 we will find a link between screening ID and CTU ID 
# screening values:
d1<-read.csv("../Data/CSV_export/CAcTUS.Table 1.csv")
colnames(d1)<-d1[1,]
d1<-d1[-1,]
d1<-d1[,c("PID","Screening ID")]
links<-merge(links,d1,by=c("Screening ID"),all.x=T,all.y=T)

write.csv(links, "data/linkingIDs.csv",row.names = F)
colnames(links)[1]<-"ID"

orig.df<-merge(orig.df,links[,c(1,3)],by=c("ID"),all.x=T,all.y=T)
colnames(links)[1]<-"DID"
colnames(links)[2]<-"ID"
orig.df<-merge(orig.df,links[,c(2,3)],by=c("ID"),all.x=T,all.y=T)

orig.df$PID<-orig.df$PID.x
orig.df$PID[is.na(orig.df$PID)]<-orig.df$PID.y[is.na(orig.df$PID)]
orig.df<-orig.df[is.na(orig.df$`Sample Date`)==F,]

colnames(orig.df)
orig.df<-orig.df[,c(1:15,18)]

write.csv(orig.df,"data/CAcTUS.Table 31.csv",row.names=F)

# how many are screening samples
orig.df[orig.df$Visit=="Screening",c("PID")]
unique(orig.df[orig.df$Visit=="Screening",c("PID")])
sum(orig.df$LLOQ=="<LLOQ")
sum(orig.df$LLOD=="<LLOD")

sum(orig.df$LLOQ[orig.df$Visit!="Screening"]=="<LLOQ")
sum(orig.df$LLOD[orig.df$Visit!="Screening"]=="<LLOD")

# lets do a histogram of all the red samples
pdf("figures/Subfig1B.pdf", width = 8, height = 6)
plt<-ggplot(unique(orig.df[orig.df$Red=="Yes",c("PID","Working days from Receipt to report")]),
       aes(x=`Working days from Receipt to report`)) + 
  geom_histogram(stat="count",fill="blue",alpha=0.3)+
  xlab("Days to return ctDNA result")+theme_bw(base_size=14)+
  ylab("Frequency")

orig.df<-arrange(orig.df,PID,`Sample Date`)
dev.off()
write.csv(plt$data, "data/data_SupFig1B.csv", row.names = FALSE)



# bring in table 1 data ---------------------------------------------------


# we will now bring in randomisation date and treatmnet arm
d1<-read.csv("../Data/CSV_export/CAcTUS.Table 1.csv")
colnames(d1)<-d1[1,]
d1<-d1[-1,]
colnames(d1)
d1<-d1[d1$`Eligible for study`=="01|Yes",c("PID","Dab and Tram immediately","Date of randomisation",
          "Investigator first line treatment choice (Arm A only)",
          "Treatment arm")]
d1$Treatment<-"Targeted"
d1$Treatment[d1$`Investigator first line treatment choice (Arm A only)`=="02|Nivolumab and Ipilimumab"]<-"IO"
d1<-d1[,c("PID","Date of randomisation",
          "Treatment arm","Treatment","Dab and Tram immediately")]
colnames(d1)[3]<-"ARM"
orig.df<-merge(orig.df,d1,by=c("PID"),all.x=F,all.y=T)

# calculate days form randomisation
orig.df$`Date of randomisation`<-as.Date(orig.df$`Date of randomisation`,
                                         "%d/%m/%Y")

orig.df$`Sample Date`<-as.Date(orig.df$`Sample Date`,
                                         "%Y-%m-%d")

orig.df$Days<-as.numeric(orig.df$`Sample Date`-orig.df$`Date of randomisation`)

# note, patients started on D+T straight away whereas those on IO have C1D1 value given
# we need to adjust the Time value accordingly, so go through each patient


# also the data-set does tell us when therapy changed so we can use that to say under what treatment
orig.df<-arrange(orig.df,PID,`Sample Date`)

colnames(orig.df)[19]<-"Firstlinetreatment"
orig.df$Treatment<-orig.df$Firstlinetreatment
unique(orig.df$Visit)
orig.df$Treatment[str_detect(orig.df$Visit,"N")]<-"IO"


# now we need baseline dates for each patient

d3<-read.csv("../Data/CSV_export/CAcTUS.Table 3.csv")
colnames(d3)<-d3[1,]
d3<-d3[-1,]
d3<-d3[,c("PID","Baseline visit date")]

orig.df<-merge(orig.df,d3,by=c("PID"),all.x=T,all.y=F)


d2<-read.csv("../Data/CSV_export/CAcTUS.Table 2.csv")
colnames(d2)<-d2[1,]
d2<-d2[-1,]
d2<-d2[d2$`Dab dispensed date`!="8888-88-88",c("PID","Dab dispensed date")]

orig.df<-merge(orig.df,d2,by=c("PID"),all.x=T,all.y=F)

d11<-read.csv("../Data/CSV_export/CAcTUS.Table 11.csv")
colnames(d11)<-d11[1,]
d11<-d11[-1,]
d11<-d11[d11$`Treatment cycle`=="01|Cycle 01",c("PID","Visit date (Immune)")]

orig.df<-merge(orig.df,d11,by=c("PID"),all.x=T,all.y=F)
orig.df$`Baseline visit date`[orig.df$Firstlinetreatment=="Targeted" & is.na(orig.df$`Baseline visit date`)]<-
  orig.df$`Dab dispensed date`[orig.df$Firstlinetreatment=="Targeted" & is.na(orig.df$`Baseline visit date`)]

orig.df$`Baseline visit date`[orig.df$Firstlinetreatment=="IO"]<-
  orig.df$`Visit date (Immune)`[orig.df$Firstlinetreatment=="IO"]

orig.df<-orig.df[,c(1:23)]
orig.df$`Baseline visit date`<-as.Date(orig.df$`Baseline visit date`,"%d/%m/%Y")
orig.df$Days2<-as.numeric(orig.df$`Sample Date`-orig.df$`Baseline visit date`)
orig.df<-arrange(orig.df,PID,Days)

# patients 21, 5 and 8 don't have a 0 time
seq(1,21,by=1)[!(seq(1,21,by=1) %in% as.numeric(unique(orig.df$PID[orig.df$Days2==0])))]

orig.df[orig.df$PID==5,c("Visit","Days2")]
orig.df$Days2[orig.df$PID==5 & orig.df$Visit=="Screening"]<-0

orig.df[orig.df$PID==8,c("Visit","Days2")]
orig.df$Days2[orig.df$PID==8]<-orig.df$Days2[orig.df$PID==8]+1

orig.df[orig.df$PID==21,c("Visit","Days2")]
orig.df$Days2[orig.df$PID==21 & orig.df$Visit=="Screening"]<-0

# patient 9 was teh odd one with screening/baseline issue
orig.df[orig.df$PID==9,c("Visit","Percentage Mutant VAF","Days2")]
# we will remove the baseline and use screening values
orig.df$Days2[orig.df$PID==9 & orig.df$Visit=="Baseline"]<- -15
orig.df$Days2[orig.df$PID==9 & orig.df$Visit=="Screening"]<- 0


# quick plot...
ggplot(orig.df,aes(Days2,`Percentage Mutant VAF`,group=PID,col=Firstlinetreatment))+
  geom_line(alpha=0.5)+xlim(0,63)+facet_wrap(~Firstlinetreatment)
# looks good we now need to put the mean value in on days when two samples were recorded...

for (ii in 2:length(orig.df$PID)){
  if (orig.df$PID[ii]==orig.df$PID[ii-1] & orig.df$`Sample Date`[ii]==orig.df$`Sample Date`[ii-1]){
    orig.df$`Percentage Mutant VAF`[ii]<-(orig.df$`Percentage Mutant VAF`[ii]+
      orig.df$`Percentage Mutant VAF`[ii-1])/2
    orig.df$`Percentage Mutant VAF`[ii-1]<-orig.df$`Percentage Mutant VAF`[ii]
  }
}

colnames(orig.df)

orig.df<-unique(orig.df[,c(1,3,4,7,17,18,19,21:24)])

# quick plot...
ggplot(orig.df,aes(Days2,`Percentage Mutant VAF`,group=PID,col=Firstlinetreatment))+
  geom_line(alpha=0.5)+xlim(0,63)+facet_wrap(~Firstlinetreatment)

# now we have a ctDNA data-set we can work with
dummy<-orig.df[orig.df$Days2==0,c("PID","Percentage Mutant VAF")]
colnames(dummy)[2]<-"BSL"

orig.df<-merge(orig.df,dummy,by=c("PID"),all.x=T,all.y=T)
orig.df$CHG<-orig.df$`Percentage Mutant VAF`-orig.df$BSL
orig.df$PCHG<-100*orig.df$CHG/orig.df$BSL

orig.df<-orig.df[orig.df$Days2>=0,]
orig.df$ARM[orig.df$ARM=="01|ARM A"]<-"ARM A"
orig.df$ARM[orig.df$ARM=="02|ARM B"]<-"ARM B"


# ctDNA plots -------------------------------------------------------------
# 
# Arm A #FFA741
# Arm B #3684BC

# #087383 Targeted
# #F37B6A IO


# CR #34C59A
# PR #FFA741
# SD #81B4DA
# NE #7F93A6



ggplot(orig.df,aes(Days2,`Percentage Mutant VAF`,group=PID,col=ARM))+
  geom_line(alpha=0.7)+facet_wrap(~ARM)+
  ylab("Mutant VAF [%]")+
  theme_bw(base_size = 14)+
  xlab("Days from Baseline")+
  scale_color_manual(values = c("ARM A" = "#FFA741", "ARM B" = "#3684BC"))

ggplot(orig.df,aes(Days2,PCHG,group=PID,col=ARM))+
  geom_line(alpha=0.7)+facet_wrap(~ARM)+
  ylab("% Change from Baseline Mutant VAF")+
  theme_bw(base_size = 14)+
  xlab("Days from Baseline")+
  geom_hline(yintercept = -80,lty=2,col=1)+
  scale_color_manual(values = c("ARM A" = "#FFA741", "ARM B" = "#3684BC"))


dummy<-orig.df[orig.df$Firstlinetreatment==orig.df$Treatment,]
dummy <- dummy %>%
  group_by(PID) %>%
  filter(row_number() <= which(PCHG < -80)[1]) %>%
  mutate(Firstlinetreatment = case_when(
    Firstlinetreatment == "IO" ~ "CPI",
    Firstlinetreatment == "Targeted" ~ "TT"))

# figure 1C
dummy_last_point <- dummy %>%
  group_by(PID, Firstlinetreatment) %>%
  filter(Days2 == max(Days2)) %>% # Filter to keep only the row with the maximum Days2 for each group
  ungroup() # Don't forget to ungroup if you're doing further operations

pdf("figures/fig1C.pdf", width = 8, height = 6)
plt<-ggplot(dummy,
       aes(Days2,PCHG,group=PID,col=Firstlinetreatment))+
  geom_line(alpha=0.5)+facet_wrap(~Firstlinetreatment)+
  geom_point(data = dummy_last_point, aes(Days2, PCHG, col = Firstlinetreatment), size = 1, shape = 19) +
  geom_hline(yintercept = -80,lty=2)+
  ylab("% Change from Baseline Mutant VAF")+
  scale_y_continuous(breaks=seq(-100,0,by=20),lim=c(-100,0))+
  scale_x_continuous(breaks=seq(0,63,by=7),lim=c(0,63))+
  theme_bw(base_size = 14)+
  xlab("Days from Baseline")+
  scale_color_manual(values = c("TT" = "#087383", "CPI" = "#F37B6A"))+
  theme(legend.position = "none")
dev.off()
write.csv(plt$data, "data/data_Fig1C.csv", row.names = FALSE)


orig.df$PID<-as.numeric(orig.df$PID)
orig.df<-arrange(orig.df,PID,Days2)

# lets assess time to reach <-80
result <- orig.df[orig.df$PCHG <= -80, ]
result <- result[!duplicated(result$PID), c("PID", "Days2","ARM","Firstlinetreatment")]
summary(result$Days2[result$ARM=="ARM B"])
summary(result$Days2[result$ARM=="ARM A" & result$Firstlinetreatment=="Targeted"])
summary(result$Days2[result$ARM=="ARM A" & result$Firstlinetreatment=="IO"])

# we now need to create visits - we already have them though...
orig.df$Visit
# "Baseline" are D+T Week 2

orig.df$VisitNew<-NA
orig.df$VisitNew[orig.df$Days2==0]<-1

orig.df$VisitNew[orig.df$Visit=="D+T Week 2" & orig.df$Firstlinetreatment=="Targeted"|
                   orig.df$Visit=="N+I C2D1" & orig.df$Firstlinetreatment=="IO"]<-2

orig.df$VisitNew[orig.df$Visit=="D+T Week 4" & orig.df$Firstlinetreatment=="Targeted"|
                   orig.df$Visit=="N+I C3D1" & orig.df$Firstlinetreatment=="IO"]<-3

orig.df$VisitNew[orig.df$Visit=="D+T Week 8" & orig.df$Firstlinetreatment=="Targeted"|
                   orig.df$Visit=="N+I C4D1" & orig.df$Firstlinetreatment=="IO"]<-4

sum(orig.df$VisitNew==2,na.rm=T)

require(EnvStats)
ggplot(orig.df[orig.df$Firstlinetreatment==orig.df$Treatment & is.na(orig.df$VisitNew)==F,])+
  geom_line(aes(VisitNew,PCHG,group=PID,col=Firstlinetreatment),alpha=0.1)+facet_wrap(~Firstlinetreatment)+
  geom_hline(yintercept = -80,lty=2)+
  ylab("% Change from Baseline Mutant VAF")+
  scale_y_continuous(breaks=seq(-100,50,by=20),lim=c(-110,50))+
  scale_x_continuous(breaks=seq(1,4,by=1),lim=c(1,4))+
  theme_bw(base_size = 14)+
  xlab("Visit")+
  scale_color_manual(values = c("Targeted" = "#087383", "IO" = "#F37B6A"))+
  stat_summary(fun = median,  aes(x=VisitNew,y=PCHG,group = VisitNew,col=Firstlinetreatment)) +  # Median line
  stat_summary(fun.data = function(y) {  # IQR bars
    return(data.frame(y = median(y), ymin = quantile(y, 0.25), ymax = quantile(y, 0.75)))
  },  aes(x=VisitNew,y=PCHG,group = VisitNew,col=Firstlinetreatment),
  geom = "errorbar", width = 0.2, size = 1)+
  stat_summary(fun = median,  aes(x=VisitNew,y=PCHG,col=Firstlinetreatment),geom="line")+
  stat_n_text(aes(x=VisitNew,y=PID),y.pos = -110)

ggplot(orig.df[orig.df$Firstlinetreatment==orig.df$Treatment & is.na(orig.df$VisitNew)==F,])+
  geom_line(aes(VisitNew,`Percentage Mutant VAF`,group=PID,col=Firstlinetreatment),alpha=0.1)+facet_wrap(~Firstlinetreatment)+
  ylab("Mutant VAF [%]")+
  scale_x_continuous(breaks=seq(1,4,by=1),lim=c(1,4))+
  theme_bw(base_size = 14)+
  xlab("Visit")+
  scale_color_manual(values = c("Targeted" = "#087383", "IO" = "#F37B6A"))+
  stat_summary(fun = median,  aes(x=VisitNew,y=`Percentage Mutant VAF`,group = VisitNew,col=Firstlinetreatment)) +  # Median line
  stat_summary(fun.data = function(y) {  # IQR bars
    return(data.frame(y = median(y), ymin = quantile(y, 0.25), ymax = quantile(y, 0.75)))
  },  aes(x=VisitNew,y=`Percentage Mutant VAF`,group = VisitNew,col=Firstlinetreatment),
  geom = "errorbar", width = 0.2, size = 1)+
  stat_summary(fun = median,  aes(x=VisitNew,y=`Percentage Mutant VAF`,col=Firstlinetreatment),geom="line")+
  stat_n_text(aes(x=VisitNew,y=PID),y.pos = -10)

# we will now add the treatment start and end dates - we know the start as being 
# the baseline date and so we can go into table 7 and 11, merge them to figure out
# when each one started and stopped


# add treatment periods ---------------------------------------------------



d7<-read.csv("../Data/CSV_export/CAcTUS.Table 7.csv")
colnames(d7)<-d7[1,]
d7<-d7[-1,]
colnames(d7)

d11<-read.csv("../Data/CSV_export/CAcTUS.Table 11.csv")
colnames(d11)<-d11[1,]
d11<-d11[-1,]
colnames(d11)


d7<-d7[,c(1,3,4)]
d7
d11<-d11[,c(1,2,32)]
d11<-d11[d11$`Niv Infusion date`!="8888-88-88",]
colnames(d11)[c(2,3)]<-c("VISIT","Date")
colnames(d7)[c(2,3)]<-c("VISIT","Date")
d11$Treatment<-"IO"
d7$Treatment<-"Targeted"



trt<-rbind(d7,d11)
trt$PID<-as.numeric(trt$PID)
trt$Date<-as.Date(trt$Date,"%d/%m/%Y")
trt<-arrange(trt,PID,Date)

trt<-trt[,c(1,3,4)]
colnames(orig.df)[2]<-"Date"
orig.df<-merge(orig.df[,c(1,2,4,5,6,10,11)],trt,by=c("PID","Date"),
               all.x=T,all.y=T)

# we've now put in the right treatment dates. what we need to do is fill in the
# blanks
# we now want to add the radiological data to this so column called SUM for diameters
# and a column called LDH for LDH. Now we hope the dates work well...


# add RECIST --------------------------------------------------------------



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
d6<-d6[,c("PID","Schedule","TL1 mm","TL2 mm","TL3 mm","TL4 mm","TL5 mm","Overall response","Treatment","Date")]
d6$SUM<-NA
d6$`TL1 mm`<-as.numeric(d6$`TL1 mm`)
d6$`TL2 mm`<-as.numeric(d6$`TL2 mm`)
d6$`TL3 mm`<-as.numeric(d6$`TL3 mm`)
d6$`TL4 mm`<-as.numeric(d6$`TL4 mm`)
d6$`TL5 mm`<-as.numeric(d6$`TL5 mm`)
d6$SUM<-rowSums(d6[,c("TL1 mm","TL2 mm","TL3 mm","TL4 mm","TL5 mm")],na.rm=T)

dummy<-d6[d6$Schedule=="00|Baseline",c("PID","SUM")]


colnames(dummy)[2]<-"BSL"
d6<-merge(d6,dummy[,c("PID","BSL")],by=c("PID"),all.x=T,all.y=T)
d6$CHG<-d6$SUM-d6$BSL
d6$PCHG<-100*d6$CHG/d6$BSL
unique(d6$Treatment)
d6$Treatment[d6$Treatment=="88|Not applicable"]<-"00|Baseline"

colnames(d6)
d6$Date<-as.Date(d6$Date,"%d/%m/%Y")
colnames(d6)[2]<-"Visit"
colnames(d6)[12]<-"BSLmm"
# need to add in the progression dates....

d6$`Overall response`[d6$PID=="4" & d6$Date=="2019-12-23"]<-"04|PD"

# patient 9 had PD on 10/11/2021
d6[d6$PID=="9",]
dummy<-d6[d6$PID=="9" & d6$Visit=="09|09",]
dummy$Visit<-"10|10"
dummy$`TL1 mm`<-NA
dummy$`TL2 mm`<-NA
dummy$`Overall response`<-"04|PD"
dummy$Treatment<-"01|D + T 1"
dummy$SUM<-NA
dummy$CHG<-NA
dummy$PCHG<-NA
dummy$Date<-"2021-11-10"
d6<-rbind(d6,dummy)


d6[d6$PID=="14",]
dummy<-d6[d6$PID=="14" & d6$Visit=="10|10",]
dummy$`TL1 mm`<-NA
dummy$`TL2 mm`<-NA
dummy$`TL3 mm`<-NA
dummy$`TL4 mm`<-NA
dummy$`TL5 mm`<-NA
dummy$`Overall response`<-"04|PD"
dummy$Treatment<-"02|Immune"
dummy$SUM<-NA
dummy$CHG<-NA
dummy$PCHG<-NA
dummy$Date<-"2022-02-21"
d6<-rbind(d6,dummy)


orig.df<-merge(orig.df,d6[,c("PID","Date","Overall response","SUM")],
            by=c("PID","Date"),all.x=T,all.y=T)



# add LDH -----------------------------------------------------------------


# we now bring in LDH...
d25<-read.csv("../Data/CSV_export/CAcTUS.Table 25.csv")

colnames(d25)<-d25[1,]
d25<-d25[-1,]
colnames(d25)
d25<-d25[,c("PID","Sample date","LDH result")]
colnames(d25)[c(2)]<-c("Date")
d25$Date<-as.Date(d25$Date,"%d/%m/%Y")
d25$`LDH result`[d25$`LDH result`=="9998"]<-NA
orig.df<-merge(orig.df,d25,
               by=c("PID","Date"),all.x=T,all.y=T)


colnames(orig.df)

# lets fill in the missing data...
id<-unique(orig.df$PID)
for(ii in 1:length(id)){
  dummy<-orig.df[orig.df$PID==id[ii],]
  orig.df$ARM[orig.df$PID==id[ii]]<-unique(dummy$ARM[is.na(dummy$ARM)==F])
  orig.df$`Baseline visit date`[orig.df$PID==id[ii]]<-unique(dummy$`Baseline visit date`[is.na(dummy$`Baseline visit date`)==F])
}

# now fill in treatment
for (ii in 2:length(orig.df$PID)){
  if(is.na(orig.df$Treatment[ii]) & orig.df$PID[ii]==orig.df$PID[ii-1]){
    orig.df$Treatment[ii]<-orig.df$Treatment[ii-1]
  }
}


for (ii in length(orig.df$PID):2){
  if(is.na(orig.df$Treatment[ii-1]) & orig.df$PID[ii]==orig.df$PID[ii-1]){
    orig.df$Treatment[ii-1]<-orig.df$Treatment[ii]
  }
}


orig.df$Days2[is.na(orig.df$Days2)]<-orig.df$Date[is.na(orig.df$Days2)]-orig.df$`Baseline visit date`[is.na(orig.df$Days2)]
orig.df$Days2[orig.df$Days2<0]<-0
orig.df$Days2[orig.df$Days2<=10]<-0 # these can be set to 0

# note when we do a switch from IO to targeted there should be a 21 day gap
# so if there is a "targeted" label within 21 days then that needs changing to IO

# for (ii in 2:length(orig.df$PID)){
#   if(orig.df$PID[ii]==orig.df$PID[ii-1] & orig.df$Treatment[ii]=="Targeted" & orig.df$Treatment[ii-1]=="IO" &
#      (orig.df$Days2[ii]-orig.df$Days2[ii-1])<21){
#     orig.df$Treatment[ii]<-"IO"
#   }
# }

# when there is a PD and the previous treatment was IO


orig.df$`Overall response`[orig.df$`Overall response`=="88|Not applicable"]<-NA


# # note, patient 19 had triple combo
# # we will call them IO
# orig.df$Treatment[orig.df$PID==19][9]<-"IO"
orig.df$`Overall response`[orig.df$PID=="6" & orig.df$Days2==512]<-"02|PR"



dummy<-orig.df[is.na(orig.df$`Overall response`)==F & 
                 orig.df$`Overall response`=="01|CR"|
                 orig.df$`Overall response`=="02|PR",c("Days2","PID","ARM")]
result<-dummy[is.na(dummy$PID)==F,]
# result <- dummy %>%
#   group_by(PID) %>%
#   slice(1) %>%
#   ungroup()

dummy<-orig.df[is.na(orig.df$`Overall response`)==F & 
                 orig.df$`Overall response`=="02|PR",c("Days2","PID","ARM")]
result2<-dummy[is.na(dummy$PID)==F,]
# result2 <- dummy %>%
#   group_by(PID) %>%
#   slice(1) %>%
#   ungroup()

dummy<-orig.df[is.na(orig.df$`Overall response`)==F & 
                 orig.df$`Overall response`=="03|SD",c("Days2","PID","ARM")]
result3<-dummy[is.na(dummy$PID)==F,]

dummy<-orig.df[is.na(orig.df$`Overall response`)==F & 
                 orig.df$`Overall response`=="04|PD",c("Days2","PID","ARM")]
result4<-dummy[is.na(dummy$PID)==F,]

dummy<-unique(orig.df[,c("PID","Treatment","Days2","ARM")])
dummy <- dummy %>%
  group_by(PID) %>%
  mutate(Start = lag(Days2, default = 0), End = Days2) %>%
  ungroup()
dummy<-dummy[dummy$Days2!=0,]
#IO will be #FF0000 and Targeted will be #0000FF 


orig.df$PID<-as.numeric(orig.df$PID)
dummy$PID<-as.numeric(dummy$PID)
result$PID<-as.numeric(result$PID)
result2$PID<-as.numeric(result2$PID)
result3$PID<-as.numeric(result3$PID)
result4$PID<-as.numeric(result4$PID)

# patient 18 days 42 shoudl eb IO as it takes 3 weeks before you can start targeted
# drugs were dispensed but they are still IO

d7<-read.csv("../Data/CSV_export/CAcTUS.Table 7.csv")
colnames(d7)<-d7[1,]
d7<-d7[-1,]
colnames(d7)

# patient 1 is odd as thetd ispensed targeted even though teh patient had progressed radiologically

dummy<-d7[d7$PID==18,]
orig.df[orig.df$PID==18,]
# they have PD but that should be associated with IO as targeted had just beend dispensed
orig.df$Treatment[orig.df$Days2<=44 & orig.df$PID==18]<-"IO"



dummy<-unique(orig.df[,c("PID","Treatment","Days2","ARM")])
dummy <- dummy %>%
  group_by(PID) %>%
  mutate(Start = lag(Days2, default = 0), End = Days2) %>%
  ungroup()
dummy<-dummy[dummy$Days2!=0,]

# this is now correct

orig.df$`LDH result`<-as.numeric(orig.df$`LDH result`)




# individual patient plots ------------------------------------------------
# Arm A #FFA741
# Arm B #3684BC

# #087383 Targeted
# #F37B6A IO


# CR #34C59A
# PR #FFA741
# SD #81B4DA
# NE #7F93A6

# we want to do side-by-side plots for each patients - best to do this by arm as
# a landscape plot
# edit dummy, for when pateint off treatment - using excel files

dummy$Treatment[dummy$PID==1][c(8,9)]<-"Off"
dummy$Treatment[dummy$PID==19][c(13)]<-"Off"

pdf("figures/Subfig3A1.pdf", width = 20, height = 2.1)
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM A" ,],
                   aes(Days2,`Percentage Mutant VAF`,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM A" ,],
             aes(Days2,`Percentage Mutant VAF`,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM A" ,],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
  geom_vline(data=result[result$ARM=="ARM A" ,],
             aes(xintercept=Days2),col="#34C59A")+
  geom_vline(data=result2[result2$ARM=="ARM A" ,],
             aes(xintercept=Days2),col="#FFA741")+
  geom_vline(data=result3[result3$ARM=="ARM A",],
             aes(xintercept=Days2),col="#81B4DA")+theme_bw(base_size=14)+
  geom_vline(data=result4[result4$ARM=="ARM A",],
             aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
  ylab("Mutant VAF [%]")+
  xlab("Days from Treatment Start")+
  theme(legend.position = "none")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383", Off=NA))
dev.off()
write.csv(orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM A" ,], "data/data_SupFig3A1.csv", row.names = FALSE)


pdf("figures/Fig3A1.pdf", width = 20, height = 2.1)
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" ,],
                   aes(Days2,`Percentage Mutant VAF`,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" ,],
             aes(Days2,`Percentage Mutant VAF`,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM B" ,],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
  geom_vline(data=result[result$ARM=="ARM B" ,],
             aes(xintercept=Days2),col="#34C59A")+
  geom_vline(data=result2[result2$ARM=="ARM B" ,],
             aes(xintercept=Days2),col="#FFA741")+
  geom_vline(data=result3[result3$ARM=="ARM B",],
             aes(xintercept=Days2),col="#81B4DA")+theme_bw(base_size=14)+
  geom_vline(data=result4[result4$ARM=="ARM B",],
             aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
  ylab("Mutant VAF [%]")+
  xlab("Days from Treatment Start")+
  theme(legend.position = "none")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383", Off=NA))
dev.off()
write.csv(orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" ,], "data/data_Fig3A1.csv", row.names = FALSE)

# for Supp table 10
df_filtered_dplyr <- orig.df %>%
  filter(PID == "14") %>% # Keep only rows where PID is "P12"
  select(PID, Date, `Percentage Mutant VAF`, `Overall response`)
# Table Sup 10
# 2	No*					
#   4	Yes	Yes	63	PD	2019-10-21	2019-12-23
# 5	No**					
#   7	Yes	No				
# 11	Yes	Yes	40	PD	2021-05-05	2021-06-14
# 12	Yes	No		***		
#   13	Yes	Yes	54	PD	2021-07-14	2021-09-06
# 14	Yes	Yes	22	PR	2022-01-19	2022-02-10
# 19	Yes	Yes	35	PD	2022-07-04	2022-08-08
# 20	Yes	No		***		
#   21	Yes	Yes	7	PD	2022-08-03	2022-08-10
# as.numeric(as.Date("2019-07-25") - as.Date("2019-09-23"), units = "days")


# task 5: Arm B of the 9 patients that switched from targeted therapy to immunotherapy what was the median time it took
df_first_date <- orig.df %>%
  filter(ARM == "ARM B") %>%
  distinct(PID, Treatment, .keep_all = TRUE) %>%
  select(PID, Date, Treatment)%>%
  filter(!PID %in% c("2", "5")) 

date_diffs <- df_first_date%>%
  group_by(PID) %>%                   # Group by PID
  arrange(Date) %>%                  # Sort by Date within each group (ascending)
  summarize(
    first_date = first(Date),        # Get the first date
    last_date = last(Date),          # Get the last date
    date_difference = last(Date) - first(Date) # Calculate the difference
  ) %>%
  ungroup()        
mean(date_diffs$date_difference)
median(date_diffs$date_difference)

# Task 6: Would you be able to provide the median time from start of TT until the ctDNA rise in Arm A?  the days and range? 
orig.df%>%
  filter(ARM == "ARM A") %>%
  select(PID, Date, `Percentage Mutant VAF`, Treatment)
times<-c(115,60,250,131,181,315)
median(times)
mean(times)
sd(times)
# PID 1: 2019-09-16 - 2019-05-24
# PID 3: 2019-07-25 -2019-09-23
# PID 6: 2021-03-22 - 2020-07-15
# PID 8: 2020-12-03 - 2021-04-13 
# PID 9: 2021-02-02 - 2021-08-02
# PID 10: 2022-01-10 - 2021-03-01

orig.df %>%
  filter(PID %in% c("5","12")) %>%
  select(c("PID","Treatment",`Percentage Mutant VAF`))

# task 8: 
df_final_output<-orig.df%>%
  filter(!is.na(`Percentage Mutant VAF`))%>%
  select(PID, Days2, `Percentage Mutant VAF`, Treatment)%>%
  group_by(PID) %>% # Regroup by PID to perform calculations per PID
  mutate(
    First_VAF = first(`Percentage Mutant VAF`), # Get the first VAF for the current PID
    `Percentage Mutant VAF (Normalized)` = (`Percentage Mutant VAF` / First_VAF) * 100 # Calculate percentage
  ) %>%
  ungroup()%>%
  select(PID, "Days2", `Percentage Mutant VAF (Normalized)`) %>%
  filter(`Percentage Mutant VAF (Normalized)`<20) %>%
  distinct(PID,  .keep_all = TRUE)
write.csv(df_final_output, file = "data/values_days.csv", row.names = FALSE)

# 
# ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM A" & orig.df$PID<10,],
#                    aes(Days2,`Percentage Mutant VAF`,group=PID))+
#   geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM A" & orig.df$PID<10,],
#              aes(Days2,`Percentage Mutant VAF`,group=PID))+
#   facet_wrap(~PID,scales="free",nrow=1)+
#   geom_rect(data=dummy[dummy$ARM=="ARM A" & dummy$PID<10,],
#             aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
#   geom_vline(data=result[result$ARM=="ARM A" & result$PID<10,],
#              aes(xintercept=Days2),col="green")+
#   geom_vline(data=result2[result2$ARM=="ARM A" & result2$PID<10,],
#              aes(xintercept=Days2),col="blue")+
#   geom_vline(data=result3[result3$ARM=="ARM A" & result3$PID<10,],
#              aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
#   ylab("Mutant VAF [%]")+
#   xlab("Days from Treatment Start")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c(IO="red",Targeted="blue"))
# 
# 
# ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM A" & orig.df$PID>=10,],
#                    aes(Days2,`Percentage Mutant VAF`,group=PID))+
#   geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM A" & orig.df$PID>=10,],
#              aes(Days2,`Percentage Mutant VAF`,group=PID))+
#   facet_wrap(~PID,scales="free",nrow=1)+
#   geom_rect(data=dummy[dummy$ARM=="ARM A" & dummy$PID>=10,],
#             aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
#   geom_vline(data=result[result$ARM=="ARM A" & result$PID>=10,],
#              aes(xintercept=Days2),col="green")+
#   geom_vline(data=result2[result2$ARM=="ARM A" & result2$PID>=10,],
#              aes(xintercept=Days2),col="blue")+
#   geom_vline(data=result3[result3$ARM=="ARM A" & result3$PID>=10,],
#              aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
#   ylab("Mutant VAF [%]")+
#   xlab("Days from Treatment Start")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c(IO="red",Targeted="blue"))

pdf("figures/Subfig3A2.pdf", width = 20, height = 2.1)
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM A" ,],
                   aes(Days2,SUM,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM A" ,],
             aes(Days2,SUM,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM A" ,],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
  geom_vline(data=result[result$ARM=="ARM A" ,],
             aes(xintercept=Days2),col="#34C59A")+
  geom_vline(data=result2[result2$ARM=="ARM A" ,],
             aes(xintercept=Days2),col="#FFA741")+
  geom_vline(data=result3[result3$ARM=="ARM A" ,],
             aes(xintercept=Days2),col="#81B4DA")+
  geom_vline(data=result4[result4$ARM=="ARM A" ,],
             aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
  ylab("SLD [mm]")+
  xlab("Days from Treatment Start")+
  theme(legend.position = "none")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383"))
dev.off()
write.csv(orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM A" ,], "data/data_SupFig3A2.csv", row.names = FALSE)


pdf("figures/Fig3A2.pdf", width = 20, height = 2.1)
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM B" ,],
                   aes(Days2,SUM,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM B" ,],
             aes(Days2,SUM,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM B" ,],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
  geom_vline(data=result[result$ARM=="ARM B" ,],
             aes(xintercept=Days2),col="#34C59A")+
  geom_vline(data=result2[result2$ARM=="ARM B" ,],
             aes(xintercept=Days2),col="#FFA741")+
  geom_vline(data=result3[result3$ARM=="ARM B" ,],
             aes(xintercept=Days2),col="#81B4DA")+
  geom_vline(data=result4[result4$ARM=="ARM B" ,],
             aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
  ylab("SLD [mm]")+
  xlab("Days from Treatment Start")+
  theme(legend.position = "none")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383"))
dev.off()

write.csv(orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM B" ,], "data/data_Fig3A2.csv", row.names = FALSE)
# ggplot()+geom_line(data=orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM A" & orig.df$PID<10,],
#                    aes(Days2,SUM,group=PID))+
#   geom_point(data=orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM A" & orig.df$PID<10,],
#              aes(Days2,SUM,group=PID))+
#   facet_wrap(~PID,scales="free",nrow=1)+
#   geom_rect(data=dummy[dummy$ARM=="ARM A" & dummy$PID<10,],
#             aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
#   geom_vline(data=result[result$ARM=="ARM A" & result$PID<10,],
#              aes(xintercept=Days2),col="green")+
#   geom_vline(data=result2[result2$ARM=="ARM A" & result2$PID<10,],
#              aes(xintercept=Days2),col="blue")+
#   geom_vline(data=result3[result3$ARM=="ARM A" & result3$PID<10,],
#              aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
#   ylab("SD [mm]")+
#   xlab("Days from Treatment Start")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c(IO="red",Targeted="blue"))
# 
# 
# ggplot()+geom_line(data=orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM A" & orig.df$PID>=10,],
#                    aes(Days2,SUM,group=PID))+
#   geom_point(data=orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM A" & orig.df$PID>=10,],
#              aes(Days2,SUM,group=PID))+
#   facet_wrap(~PID,scales="free",nrow=1)+
#   geom_rect(data=dummy[dummy$ARM=="ARM A" & dummy$PID>=10,],
#             aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
#   geom_vline(data=result[result$ARM=="ARM A" & result$PID>=10,],
#              aes(xintercept=Days2),col="green")+
#   geom_vline(data=result2[result2$ARM=="ARM A" & result2$PID>=10,],
#              aes(xintercept=Days2),col="blue")+
#   geom_vline(data=result3[result3$ARM=="ARM A" & result3$PID>=10,],
#              aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
#   ylab("SD [mm]")+
#   xlab("Days from Treatment Start")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c(IO="red",Targeted="blue"))

orig.df$`LDH result`<-as.numeric(orig.df$`LDH result`)
pdf("figures/Subfig3A3.pdf", width = 20, height = 2.1)
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM A",],
                   aes(Days2,`LDH result`,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM A",],
             aes(Days2,`LDH result`,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM A",],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
  geom_vline(data=result[result$ARM=="ARM A",],
             aes(xintercept=Days2),col="#34C59A")+
  geom_vline(data=result2[result2$ARM=="ARM A",],
             aes(xintercept=Days2),col="#FFA741")+
  geom_vline(data=result3[result3$ARM=="ARM A",],
             aes(xintercept=Days2),col="#81B4DA")+
  geom_vline(data=result4[result4$ARM=="ARM A",],
             aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
  ylab("LDH [IU/L]")+
  xlab("Days from Treatment Start")+
  theme(legend.position = "none")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383"))
dev.off()

write.csv(orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM A",], "data/data_SupFig3A3.csv", row.names = FALSE)

pdf("figures/fig3A3.pdf", width = 20, height = 2.1)
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM B",],
                   aes(Days2,`LDH result`,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM B",],
             aes(Days2,`LDH result`,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM B",],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
  geom_vline(data=result[result$ARM=="ARM B",],
             aes(xintercept=Days2),col="#34C59A")+
  geom_vline(data=result2[result2$ARM=="ARM B",],
             aes(xintercept=Days2),col="#FFA741")+
  geom_vline(data=result3[result3$ARM=="ARM B",],
             aes(xintercept=Days2),col="#81B4DA")+
  geom_vline(data=result4[result4$ARM=="ARM B",],
             aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
  ylab("LDH [IU/L]")+
  xlab("Days from Treatment Start")+
  theme(legend.position = "none")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383"))
dev.off()

write.csv(orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM B",], "data/data_Fig3A3.csv", row.names = FALSE)
# ggplot()+geom_line(data=orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM A" & orig.df$PID<10,],
#                    aes(Days2,`LDH result`,group=PID))+
#   geom_point(data=orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM A" & orig.df$PID<10,],
#              aes(Days2,`LDH result`,group=PID))+
#   facet_wrap(~PID,scales="free",nrow=1)+
#   geom_rect(data=dummy[dummy$ARM=="ARM A" & dummy$PID<10,],
#             aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
#   geom_vline(data=result[result$ARM=="ARM A" & result$PID<10,],
#              aes(xintercept=Days2),col="green")+
#   geom_vline(data=result2[result2$ARM=="ARM A" & result2$PID<10,],
#              aes(xintercept=Days2),col="blue")+
#   geom_vline(data=result3[result3$ARM=="ARM A" & result3$PID<10,],
#              aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
#   ylab("LDH [IU/L]")+
#   xlab("Days from Treatment Start")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c(IO="red",Targeted="blue"))
# ggplot()+geom_line(data=orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM A" & orig.df$PID>=10,],
#                    aes(Days2,`LDH result`,group=PID))+
#   geom_point(data=orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM A" & orig.df$PID>=10,],
#              aes(Days2,`LDH result`,group=PID))+
#   facet_wrap(~PID,scales="free",nrow=1)+
#   geom_rect(data=dummy[dummy$ARM=="ARM A" & dummy$PID>=10,],
#             aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
#   geom_vline(data=result[result$ARM=="ARM A" & result$PID>=10,],
#              aes(xintercept=Days2),col="green")+
#   geom_vline(data=result2[result2$ARM=="ARM A" & result2$PID>=10,],
#              aes(xintercept=Days2),col="blue")+
#   geom_vline(data=result3[result3$ARM=="ARM A" & result3$PID>=10,],
#              aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
#   ylab("LDH [IU/L]")+
#   xlab("Days from Treatment Start")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c(IO="red",Targeted="blue"))

# Arm B
# 
# 
# 
# ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID<12,],
#                    aes(Days2,`Percentage Mutant VAF`,group=PID))+
#   geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID<12,],
#              aes(Days2,`Percentage Mutant VAF`,group=PID))+
#   facet_wrap(~PID,scales="free",nrow=1)+
#   geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID<12,],
#             aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
#   geom_vline(data=result[result$ARM=="ARM B" & result$PID<12,],
#              aes(xintercept=Days2),col="green")+
#   geom_vline(data=result2[result2$ARM=="ARM B" & result2$PID<12,],
#              aes(xintercept=Days2),col="blue")+
#   geom_vline(data=result3[result3$ARM=="ARM B" & result3$PID<12,],
#              aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
#   ylab("Mutant VAF [%]")+
#   xlab("Days from Treatment Start")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c(IO="red",Targeted="blue"))
# ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID>=12,],
#                    aes(Days2,`Percentage Mutant VAF`,group=PID))+
#   geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID>=12,],
#              aes(Days2,`Percentage Mutant VAF`,group=PID))+
#   facet_wrap(~PID,scales="free",nrow=1)+
#   geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID>=12,],
#             aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
#   geom_vline(data=result[result$ARM=="ARM B" & result$PID>=12,],
#              aes(xintercept=Days2),col="green")+
#   geom_vline(data=result2[result2$ARM=="ARM B" & result2$PID>=12,],
#              aes(xintercept=Days2),col="blue")+
#   geom_vline(data=result3[result3$ARM=="ARM B" & result3$PID>=12,],
#              aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
#   ylab("Mutant VAF [%]")+
#   xlab("Days from Treatment Start")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c(IO="red",Targeted="blue"))
# 
# ggplot()+geom_line(data=orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM B" & orig.df$PID<12,],
#                    aes(Days2,SUM,group=PID))+
#   geom_point(data=orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM B" & orig.df$PID<12,],
#              aes(Days2,SUM,group=PID))+
#   facet_wrap(~PID,scales="free",nrow=1)+
#   geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID<12,],
#             aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
#   geom_vline(data=result[result$ARM=="ARM B" & result$PID<12,],
#              aes(xintercept=Days2),col="green")+
#   geom_vline(data=result2[result2$ARM=="ARM B" & result2$PID<12,],
#              aes(xintercept=Days2),col="blue")+
#   geom_vline(data=result3[result3$ARM=="ARM B" & result3$PID<12,],
#              aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
#   ylab("SD [mm]")+
#   xlab("Days from Treatment Start")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c(IO="red",Targeted="blue"))
# ggplot()+geom_line(data=orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM B" & orig.df$PID>=12,],
#                    aes(Days2,SUM,group=PID))+
#   geom_point(data=orig.df[is.na(orig.df$SUM)==F & orig.df$ARM=="ARM B" & orig.df$PID>=12,],
#              aes(Days2,SUM,group=PID))+
#   facet_wrap(~PID,scales="free",nrow=1)+
#   geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID>=12,],
#             aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
#   geom_vline(data=result[result$ARM=="ARM B" & result$PID>=12,],
#              aes(xintercept=Days2),col="green")+
#   geom_vline(data=result2[result2$ARM=="ARM B" & result2$PID>=12,],
#              aes(xintercept=Days2),col="blue")+
#   geom_vline(data=result3[result3$ARM=="ARM B" & result3$PID>=12,],
#              aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
#   ylab("SD [mm]")+
#   xlab("Days from Treatment Start")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c(IO="red",Targeted="blue"))
# 
# orig.df$`LDH result`<-as.numeric(orig.df$`LDH result`)
# ggplot()+geom_line(data=orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM B" & orig.df$PID<12,],
#                    aes(Days2,`LDH result`,group=PID))+
#   geom_point(data=orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM B" & orig.df$PID<12,],
#              aes(Days2,`LDH result`,group=PID))+
#   facet_wrap(~PID,scales="free",nrow=1)+
#   geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID<12,],
#             aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
#   geom_vline(data=result[result$ARM=="ARM B" & result$PID<12,],
#              aes(xintercept=Days2),col="green")+
#   geom_vline(data=result2[result2$ARM=="ARM B" & result2$PID<12,],
#              aes(xintercept=Days2),col="blue")+
#   geom_vline(data=result3[result3$ARM=="ARM B" & result3$PID<12,],
#              aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
#   ylab("LDH [IU/L]")+
#   xlab("Days from Treatment Start")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c(IO="red",Targeted="blue"))
# ggplot()+geom_line(data=orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM B" & orig.df$PID>=12,],
#                    aes(Days2,`LDH result`,group=PID))+
#   geom_point(data=orig.df[is.na(orig.df$`LDH result`)==F & orig.df$ARM=="ARM B" & orig.df$PID>=12,],
#              aes(Days2,`LDH result`,group=PID))+
#   facet_wrap(~PID,scales="free",nrow=1)+
#   geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID>=12,],
#             aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
#   geom_vline(data=result[result$ARM=="ARM B" & result$PID>=12,],
#              aes(xintercept=Days2),col="green")+
#   geom_vline(data=result2[result2$ARM=="ARM B" & result2$PID>=12,],
#              aes(xintercept=Days2),col="blue")+
#   geom_vline(data=result3[result3$ARM=="ARM B" & result3$PID>=12,],
#              aes(xintercept=Days2),col="red")+theme_bw(base_size=14)+
#   ylab("LDH [IU/L]")+
#   xlab("Days from Treatment Start")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c(IO="red",Targeted="blue"))




# zoomed in arm B plots ---------------------------------------------------
# we want in zoomed in plots leading up to treatment switch for IO
# 4,7,11,13,14,19,20,21 from ARM B for ctDNA only

id<-4
pdf("figures/Subfig3B_n4.pdf", width = 4, height = 4)
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
                   aes(Days2,`Percentage Mutant VAF`,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
             aes(Days2,`Percentage Mutant VAF`,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID==id,],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+theme_bw(base_size=14)+
  ylab("Mutant VAF [%]")+
  xlab("Days from Treatment Start")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383",Off="blue"))+scale_x_continuous(lim=c(0,35))+
  theme(legend.position = "none")
dev.off()
write.csv(orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,], "data/data_SupFig3B_n4.csv", row.names = FALSE)

pdf("figures/Subfig3B_n7.pdf", width = 4, height = 4)
id<-7
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
                   aes(Days2,`Percentage Mutant VAF`,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
             aes(Days2,`Percentage Mutant VAF`,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID==id,],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+theme_bw(base_size=14)+
  ylab("Mutant VAF [%]")+
  xlab("Days from Treatment Start")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383"))+scale_x_continuous(lim=c(0,35))+
  theme(legend.position = "none")
dev.off()
write.csv(orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,], "data/data_SupFig3B_n7.csv", row.names = FALSE)

pdf("figures/Subfig3B_n11.pdf", width = 4, height = 4)
id<-11
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
                   aes(Days2,`Percentage Mutant VAF`,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
             aes(Days2,`Percentage Mutant VAF`,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID==id,],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+theme_bw(base_size=14)+
  ylab("Mutant VAF [%]")+
  xlab("Days from Treatment Start")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383"))+scale_x_continuous(lim=c(0,60))+
  theme(legend.position = "none")
dev.off()
write.csv(orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,], "data/data_SupFig3B_n11.csv", row.names = FALSE)

pdf("figures/Subfig3B_n13.pdf", width = 4, height = 4)
id<-13
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
                   aes(Days2,`Percentage Mutant VAF`,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
             aes(Days2,`Percentage Mutant VAF`,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID==id,],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+theme_bw(base_size=14)+
  ylab("Mutant VAF [%]")+
  xlab("Days from Treatment Start")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383"))+scale_x_continuous(lim=c(0,40))+
  theme(legend.position = "none")
dev.off()
write.csv(orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,], "data/data_SupFig3B_n13.csv", row.names = FALSE)

pdf("figures/Subfig3B_n14.pdf", width = 4, height = 4)
id<-14
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
                   aes(Days2,`Percentage Mutant VAF`,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
             aes(Days2,`Percentage Mutant VAF`,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID==id,],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+
  ylab("Mutant VAF [%]")+theme_bw(base_size=14)+
  xlab("Days from Treatment Start")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383"))+scale_x_continuous(lim=c(0,65))+
  theme(legend.position = "none")
dev.off()
write.csv(orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,], "data/data_SupFig3B_n14.csv", row.names = FALSE)

pdf("figures/Subfig3B_n19.pdf", width = 4, height = 4)
id<-19
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
                   aes(Days2,`Percentage Mutant VAF`,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
             aes(Days2,`Percentage Mutant VAF`,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID==id,],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+theme_bw(base_size=14)+
  ylab("Mutant VAF [%]")+
  xlab("Days from Treatment Start")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383",Targeted="blue"))+scale_x_continuous(lim=c(0,35))+
  theme(legend.position = "none")
dev.off()
write.csv(orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,], "data/data_SupFig3B_n19.csv", row.names = FALSE)

pdf("figures/Subfig3B_n20.pdf", width = 4, height = 4)
id<-20
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
                   aes(Days2,`Percentage Mutant VAF`,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
             aes(Days2,`Percentage Mutant VAF`,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID==id,],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+theme_bw(base_size=14)+
  ylab("Mutant VAF [%]")+
  xlab("Days from Treatment Start")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383"))+scale_x_continuous(lim=c(0,35))+
  theme(legend.position = "none")
dev.off()
write.csv(orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,], "data/data_SupFig3B_n20.csv", row.names = FALSE)

pdf("figures/Subfig3B_n21.pdf", width = 4, height = 4)
id<-21
plt<-ggplot()+geom_line(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
                   aes(Days2,`Percentage Mutant VAF`,group=PID))+
  geom_point(data=orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,],
             aes(Days2,`Percentage Mutant VAF`,group=PID))+
  facet_wrap(~PID,scales="free",nrow=1)+
  geom_rect(data=dummy[dummy$ARM=="ARM B" & dummy$PID==id,],
            aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, fill = Treatment), alpha = 0.2)+theme_bw(base_size=14)+
  ylab("Mutant VAF [%]")+
  xlab("Days from Treatment Start")+
  scale_fill_manual(values=c(IO="#F37B6A",Targeted="#087383"))+scale_x_continuous(lim=c(0,80))+
  theme(legend.position = "none")
dev.off()
write.csv(orig.df[is.na(orig.df$`Percentage Mutant VAF`)==F & orig.df$ARM=="ARM B" & orig.df$PID==id,], "data/data_SupFig3B_n21.csv", row.names = FALSE)


# add toxicity and death to create swimmer plot data-set ------------------



# we now want to add toxicity and death dates to the data frame which can then
# eb used tio create a swimmer plot
recist<-orig.df[is.na(orig.df$`Overall response`)==F,c("Days2","PID","ARM",
                                                       "Overall response")]
colnames(recist)[4]<-"Response"
recist$Response[recist$Response=="01|CR"]<-"CR"
recist$Response[recist$Response=="02|PR"]<-"PR"
recist$Response[recist$Response=="03|SD"]<-"SD"
recist$Response[recist$Response=="04|PD"]<-"PD"
recist$Response[recist$Response=="05|NE"]<-"NE"

# lets collect death dates and add them
d19<-read.csv("../Data/CSV_export/CAcTUS.Table 19.csv")
colnames(d19)<-d19[1,]
d19<-d19[-1,]
colnames(d19)

# its the 1st 21 patients we want data for
d19<-d19[c(1:21),]
# 4 column are most useful
d19<-d19[,c(1:3,7)]

# merge with d1
orig.df<-merge(orig.df,d19,by=c("PID"),all.x=T,all.y=F)

# calculate OS times and outcomes
orig.df$`Date of withdrawal`<-as.Date(orig.df$`Date of withdrawal`,"%d/%m/%Y")

orig.df$OS<-as.numeric(orig.df$`Date of withdrawal`-orig.df$`Baseline visit date`)
orig.df$DEATH<-0
orig.df$DEATH[orig.df$`Main Reason`=="03|Participant died"]<-1
orig.df$OS[orig.df$OS<=orig.df$Days2]<-orig.df$Days2[orig.df$OS<=orig.df$Days2]


os.df<-unique(orig.df[,c("PID","OS","DEATH","ARM")])
os.df$DEATH[os.df$DEATH==1]<-"Died"
os.df$DEATH[os.df$DEATH==0]<-"Alive"
os.df$PID<-factor(os.df$PID,levels=c(3,1,8,18,16,17,15,9,10,6,
                                     5,2,7,13,19,20,21,4,14,12,11))

os.df<-os.df %>%
  group_by(PID) %>%
  filter(OS == max(OS))%>%
  ungroup()


dummy<-unique(orig.df[,c("PID","Treatment","Days2","ARM")])
dummy <- dummy %>%
  group_by(PID) %>%
  mutate(Start = lag(Days2, default = 0), End = Days2) %>%
  ungroup()
dummy<-dummy[dummy$Days2!=0,]


# HERE WE CHOOSE WHETHER TO ALLOW TREATMENT COLOURS TO GO UP TO ALIVE/DEATH OR NOT
# dummy<-merge(dummy,os.df,by=c("PID","ARM"),all.x=T)
# dummy <- dummy %>%
#   group_by(PID) %>%
#   mutate(End = ifelse(row_number() == n(), OS, End)) %>%
#   ungroup()
# ENDS HERE

dummy$PID<-factor(dummy$PID,levels=c(3,1,18,16,17,15,8,9,10,6,
                                     2,7,13,21,20,5,19,14,4,12,11))
recist$PID<-factor(recist$PID,levels=c(3,1,18,16,17,15,8,9,10,6,
                                       2,7,13,21,20,5,19,14,4,12,11))


# now we add in toxicity....

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
d13$PID<-factor(d13$PID,levels=c(3,1,18,16,17,15,8,9,10,6,
                                   2,7,13,21,20,5,19,14,4,12,11))
d13<-merge(d13,unique(orig.df[,c("PID","Baseline visit date","ARM")]),
           by=c("PID"),all.x=T)
d13$`Start date`<-as.Date(d13$`Start date`,"%d/%m/%Y")
d13$Days<-as.numeric(d13$`Start date`-d13$`Baseline visit date`)
d13$Days[d13$Serious=="01|Yes"]
d13<-d13[is.na(d13$Days)==F,]
d13<-d13[d13$Days>=0,]

# swimmer plots -----------------------------------------------------------

os.df$DEATH[os.df$PID==7]<-"COVID-19-related death"
os.df$DEATH[os.df$DEATH=="Died"]<-"Melanoma-related death"
# for all those that are alive we will allow the treatment colour continue till the end
dummy <- dummy %>% mutate(Treatment = case_when(
  Treatment == "IO" ~ "CPI",
  Treatment == "Targeted" ~ "TT"))

pdf("figures/Fig1D.pdf", width = 14, height = 6.5)
plt<-ggplot() +
  geom_segment(dat=dummy, aes(xend = End, yend = PID,
                              x = Start, y = PID, color = Treatment), 
               size = 5,alpha=0.2) +
  labs(
    x = "Days from Treatment Start",
    y = "Patient ID",
    color="Treatment"
  ) +
  theme_minimal(base_size=14) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    #legend.position = "bottom",   # Position the legend below the x-axis
    legend.title = element_blank()  # Remove the legend title
  )+
  scale_color_manual(name="Treatment",values = c("TT" = "#087383", "CPI" = "#F37B6A",
                                "CR"="#34C59A","PR"="#FFA741","SD"="#81B4DA","PD"="red",
                                "NE"="#7F93A6"),
                     breaks = c("TT", "CPI", "CR","PR", "SD", "PD", "NE"))+
  scale_shape_manual(values = c("Melanoma-related death" = 16, "Alive" = 17, "COVID-19-related death" = 15))+
  geom_point(data=recist,aes(x=Days2,y=PID,col=Response),size=3)+
  facet_wrap(~ARM,scales="free_y")+
  geom_point(data=os.df,aes(x=OS,y=PID,shape=DEATH),size=3,alpha=0.3)+
  guides(
    color = guide_legend(title = "Treatment"),
    shape = guide_legend(order = 4)
  )
dev.off()
write.csv(dummy, "data/data_Fig1D.csv", row.names = FALSE)



# Arm A #FFA741
# Arm B #3684BC

# #087383 Targeted
# #F37B6A IO


# CR #34C59A
# PR #FFA741
# SD #81B4DA
# NE #7F93A6



# additional analyses -----------------------------------------------------

df<-read.csv("temp_data/prognostic_variables.csv",header=T)
df$CRP[df$CRP=="<5.0"]<-2.5
df$LDH
df$CRP<-as.numeric(df$CRP)
df$Platelet
df$Lymphocyte
df$ANC
df$ECOG
# We want to know how various aspects associated with prognosis change during the 
# targeted therapy run in in Arm B so to compare the baseline before starting immune 
# therapy (ie C1D1) vs the baseline/screening value of ECOG, LDH, ANC/lymphocyte 
# ratio, platelets, CRP. 
df$Visit<-factor(df$Visit,levels=c("Baseline/Screening - Targeted","C1D1 - IO"))
df$LDH<-as.numeric(df$LDH)
pdf("figures/Fig2B.pdf", width = 4, height = 4)
plt<-ggplot(df[df$ARM=="ARM B",])+geom_point(aes(Visit,LDH,group=PID),alpha=0.3)+
  geom_line(aes(Visit,LDH,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=LDH,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("LDH [IU/L]")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 
dev.off()
write.csv(plt$data, "data/data_Fig2B.csv", row.names = FALSE)

dfnow <- df[df$ARM=="ARM B",]
p<-wilcox.test(dfnow$LDH[dfnow$Period == 1],
               dfnow$LDH[dfnow$Period == 2], alternative = "two.sided", exact = FALSE, correct = TRUE)
print(paste0('Fig2B', ' pvalue=', p$p.value))

pdf("figures/Fig2F.pdf", width = 4, height = 4)
plt<-ggplot(df[df$ARM=="ARM B",])+geom_point(aes(Visit,CRP,group=PID),alpha=0.3)+
  geom_line(aes(Visit,CRP,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=CRP,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("CRP [mg/L]")+scale_y_log10()+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 
dev.off()
write.csv(plt$data, "data/data_Fig2F.csv", row.names = FALSE)
p<-wilcox.test(dfnow$CRP[dfnow$Period == 1],
               dfnow$CRP[dfnow$Period == 2], alternative = "two.sided", exact = FALSE, correct = TRUE)
print(paste0('Fig2F', ' pvalue=', p$p.value))
pdf("figures/Fig2G.pdf", width = 4, height = 4)
plt<-ggplot(df[df$ARM=="ARM B",])+geom_point(aes(Visit,Platelet,group=PID),alpha=0.3)+
  geom_line(aes(Visit,Platelet,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=Platelet,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("Platelet [x10^9/L]")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 
dev.off()
write.csv(plt$data, "data/data_Fig2G.csv", row.names = FALSE)
p<-wilcox.test(dfnow$Platelet[dfnow$Period == 1],
               dfnow$Platelet[dfnow$Period == 2], alternative = "two.sided", exact = FALSE, correct = TRUE)
print(paste0('Fig2G', ' pvalue=', p$p.value))
pdf("figures/Fig2E.pdf", width = 4, height = 4)
plt<-ggplot(df[df$ARM=="ARM B",])+geom_point(aes(Visit,Lymphocyte,group=PID),alpha=0.3)+
  geom_line(aes(Visit,Lymphocyte,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=Lymphocyte,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("Lymphocyte [x10^9/L]")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 
dev.off()
write.csv(plt$data, "data/data_Fig2E.csv", row.names = FALSE)
p<-wilcox.test(dfnow$Lymphocyte[dfnow$Period == 1],
               dfnow$Lymphocyte[dfnow$Period == 2], alternative = "two.sided", exact = FALSE, correct = TRUE)
print(paste0('Fig2E', ' pvalue=', p$p.value))
pdf("figures/Fig2C.pdf", width = 4, height = 4)
plt<-ggplot(df[df$ARM=="ARM B",])+geom_point(aes(Visit,ANC,group=PID),alpha=0.3)+
  geom_line(aes(Visit,ANC,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=ANC,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("ANC [x10^9/L]")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 
dev.off()
write.csv(plt$data, "data/data_Fig2C.csv", row.names = FALSE)
p<-wilcox.test(dfnow$ANC[dfnow$Period == 1],
               dfnow$ANC[dfnow$Period == 2], alternative = "two.sided", exact = FALSE, correct = TRUE)
print(paste0('Fig2C', ' pvalue=', p$p.value))
pdf("figures/Fig2D.pdf", width = 4, height = 4)
plt<-ggplot(df[df$ARM=="ARM B",])+geom_point(aes(Visit,ANC/Lymphocyte,group=PID),alpha=0.3)+
  geom_line(aes(Visit,ANC/Lymphocyte,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=ANC/Lymphocyte,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("ANC/Lymphocyte")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 
dev.off()
write.csv(plt$data, "data/data_Fig2D.csv", row.names = FALSE)
dfnow$ANC_Lympocyte<-dfnow$ANC/dfnow$Lymphocyte
p<-wilcox.test(dfnow$ANC_Lympocyte[dfnow$Period == 1],
               dfnow$ANC_Lympocyte[dfnow$Period == 2], alternative = "two.sided", exact = FALSE, correct = TRUE)
print(paste0('Fig2D', ' pvalue=', p$p.value))

ggplot(df[df$ARM=="ARM B",])+geom_point(aes(Visit,ANC/Lymphocyte,group=PID),alpha=0.3)+
  geom_line(aes(Visit,ANC/Lymphocyte,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=ANC/Lymphocyte,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("ANC/Lymphocyte")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 

df.ctdna<-data.frame(Visit=c("B","B","B","B","B","B","B","B","B",
                             "I","I","I","I","I","I","I","I","I"),
                     ctDNA=c(7.88,26.03,24,1.8,10,31,42,4.5,10,
                             0,0,0.23,0,1.3,0,0,0,0),
                     PID=c(4,7,11,12,13,14,19,20,21,
                           4,7,11,12,13,14,19,20,21))
pdf("figures/Fig2A.pdf", width = 4, height = 4)
plt<-ggplot(df.ctdna)+geom_point(aes(Visit,ctDNA,group=PID),alpha=0.3)+
  geom_line(aes(Visit,ctDNA,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=ctDNA,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("Mutant VAF [%]")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 
dev.off()
write.csv(plt$data, "data/data_Fig2A.csv", row.names = FALSE)
p<-wilcox.test(df.ctdna$ctDNA[df.ctdna$Visit == "B"],
            df.ctdna$ctDNA[df.ctdna$Visit == "I"], alternative = "two.sided", exact = FALSE, correct = TRUE)
print(paste0('Fig2A', ' pvalue=', p$p.value))

# ECOG we do as a table
mean(df$ECOG[df$ARM=="ARM B" & df$Visit=="Baseline/Screening - Targeted"]-
  df$ECOG[df$ARM=="ARM B" & df$Visit=="C1D1 - IO"])
sqrt(var(df$ECOG[df$ARM=="ARM B" & df$Visit=="Baseline/Screening - Targeted"]-
           df$ECOG[df$ARM=="ARM B" & df$Visit=="C1D1 - IO"]))
ggplot(df[df$ARM=="ARM B",])+geom_jitter(aes(Visit,ECOG,group=PID),alpha=0.3,height=0)+
  geom_line(aes(Visit,ECOG,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=ECOG,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("ECOG")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 

# Could you also look at the patients on Arm A starting 
# targeted therapy first so baseline/screening with the C1 D1 of immune therapy. 


ggplot(df[df$ARM=="ARM A" & df$First=="Targeted",])+geom_point(aes(Visit,LDH,group=PID),alpha=0.3)+
  geom_line(aes(Visit,LDH,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=LDH,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("LDH [IU/L]")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 

ggplot(df[df$ARM=="ARM A" & df$First=="Targeted",])+geom_point(aes(Visit,CRP,group=PID),alpha=0.3)+
  geom_line(aes(Visit,CRP,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=CRP,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("CRP [mg/L]")+scale_y_log10()+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 

ggplot(df[df$ARM=="ARM A" & df$First=="Targeted",])+geom_point(aes(Visit,Platelet,group=PID),alpha=0.3)+
  geom_line(aes(Visit,Platelet,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=Platelet,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("Platelet [x10^9/L]")

ggplot(df[df$ARM=="ARM A" & df$First=="Targeted",])+geom_point(aes(Visit,Lymphocyte,group=PID),alpha=0.3)+
  geom_line(aes(Visit,Lymphocyte,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=Lymphocyte,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("Lymphocyte [x10^9/L]")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 

ggplot(df[df$ARM=="ARM A" & df$First=="Targeted",])+geom_point(aes(Visit,ANC,group=PID),alpha=0.3)+
  geom_line(aes(Visit,ANC,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=ANC,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("ANC [x10^9/L]")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 

ggplot(df[df$ARM=="ARM A" & df$First=="Targeted",])+geom_point(aes(Visit,ANC/Lymphocyte,group=PID),alpha=0.3)+
  geom_line(aes(Visit,ANC/Lymphocyte,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=ANC/Lymphocyte,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("ANC/Lymphocyte")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 

df.ctdna<-data.frame(Visit=c("B","B","B",
                             "I","I","I"),
                     ctDNA=c(15,70,6,
                             20,20,25),
                     PID=c(1,8,9,
                           1,8,9))

ggplot(df.ctdna)+geom_point(aes(Visit,ctDNA,group=PID),alpha=0.3)+
  geom_line(aes(Visit,ctDNA,group=PID),alpha=0.5)+
  theme_bw(base_size=14)+stat_summary(aes(x=Visit,y=ctDNA,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("Mtant VAF [%]")+
  scale_x_discrete(labels = c("Baseline (targeted)", "Pre-Immune therapy")) 

# ECOG Table
mean(df$ECOG[df$ARM=="ARM A" & df$First=="Targeted" & df$Visit=="Baseline/Screening - Targeted"]-
       df$ECOG[df$ARM=="ARM A" & df$First=="Targeted" & df$Visit=="C1D1 - IO"])
sqrt(var(df$ECOG[df$ARM=="ARM A" & df$First=="Targeted" & df$Visit=="Baseline/Screening - Targeted"]-
           df$ECOG[df$ARM=="ARM A" & df$First=="Targeted" & df$Visit=="C1D1 - IO"]))


# In addition the patients starting IO first on arm A ie baseline/screening with 
# the C1D1 immune therapy of Arm B?

pdf("figures/Subfig2A.pdf", width = 4, height = 4)
plt<-ggplot(df[(df$First=="IO" & df$ARM=="ARM A")|(df$Visit=="C1D1 - IO" & df$ARM=="ARM B"),])+geom_point(aes(ARM,LDH),alpha=0.3)+
  theme_bw(base_size=14)+stat_summary(aes(x=ARM,y=LDH,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("LDH [IU/L]")
dev.off()
write.csv(plt$data, "data/data_SupFig2A.csv", row.names = FALSE)

dfnow<-df[(df$First=="IO" & df$ARM=="ARM A")|(df$Visit=="C1D1 - IO" & df$ARM=="ARM B"),]
p<-wilcox.test(dfnow$LDH[dfnow$ARM == 'ARM A'],
               dfnow$LDH[dfnow$ARM == 'ARM B'],  exact = FALSE, correct = TRUE)
print(paste0('SubFig2A', ' pvalue=', p$p.value))

pdf("figures/Subfig2E.pdf", width = 4, height = 4)
plt<-ggplot(df[(df$First=="IO" & df$ARM=="ARM A")|(df$Visit=="C1D1 - IO" & df$ARM=="ARM B"),])+geom_point(aes(ARM,CRP),alpha=0.3)+
  theme_bw(base_size=14)+stat_summary(aes(x=ARM,y=CRP,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("CRP [mg/L]")
dev.off()
write.csv(plt$data, "data/data_SupFig2E.csv", row.names = FALSE)
p<-wilcox.test(dfnow$CRP[dfnow$ARM == 'ARM A'],
               dfnow$CRP[dfnow$ARM == 'ARM B'],  exact = FALSE, correct = TRUE)
print(paste0('SubFig2E', ' pvalue=', p$p.value))

pdf("figures/Subfig2D.pdf", width = 4, height = 4)
plt<-ggplot(df[(df$First=="IO" & df$ARM=="ARM A")|(df$Visit=="C1D1 - IO" & df$ARM=="ARM B"),])+geom_point(aes(ARM,Platelet),alpha=0.3)+
  theme_bw(base_size=14)+stat_summary(aes(x=ARM,y=Platelet,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("Platelet [x10^9/L]")
dev.off()
write.csv(plt$data, "data/data_SupFig2D.csv", row.names = FALSE)
p<-wilcox.test(dfnow$Platelet[dfnow$ARM == 'ARM A'],
               dfnow$Platelet[dfnow$ARM == 'ARM B'],  exact = FALSE, correct = TRUE)
print(paste0('SubFig2D', ' pvalue=', p$p.value))

pdf("figures/Subfig2F.pdf", width = 4, height = 4)
plt<-ggplot(df[(df$First=="IO" & df$ARM=="ARM A")|(df$Visit=="C1D1 - IO" & df$ARM=="ARM B"),])+geom_point(aes(ARM,Lymphocyte),alpha=0.3)+
  theme_bw(base_size=14)+stat_summary(aes(x=ARM,y=Lymphocyte,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("Lymphocyte [x10^9/L]")
dev.off()
write.csv(plt$data, "data/data_SupFig2F.csv", row.names = FALSE)
p<-wilcox.test(dfnow$Lymphocyte[dfnow$ARM == 'ARM A'],
               dfnow$Lymphocyte[dfnow$ARM == 'ARM B'],  exact = FALSE, correct = TRUE)
print(paste0('SubFig2F', ' pvalue=', p$p.value))

pdf("figures/Subfig2B.pdf", width = 4, height = 4)
plt<-ggplot(df[(df$First=="IO" & df$ARM=="ARM A")|(df$Visit=="C1D1 - IO" & df$ARM=="ARM B"),])+geom_point(aes(ARM,ANC),alpha=0.3)+
  theme_bw(base_size=14)+stat_summary(aes(x=ARM,y=ANC,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("ANC [x10^9/L]")
dev.off()
write.csv(plt$data, "data/data_SupFig2B.csv", row.names = FALSE)
p<-wilcox.test(dfnow$ANC[dfnow$ARM == 'ARM A'],
               dfnow$ANC[dfnow$ARM == 'ARM B'],  exact = FALSE, correct = TRUE)
print(paste0('SubFig2B', ' pvalue=', p$p.value))
pdf("figures/Subfig2C.pdf", width = 4, height = 4)
plt<-ggplot(df[(df$First=="IO" & df$ARM=="ARM A")|(df$Visit=="C1D1 - IO" & df$ARM=="ARM B"),])+geom_point(aes(ARM,ANC/Lymphocyte),alpha=0.3)+
  theme_bw(base_size=14)+stat_summary(aes(x=ARM,y=ANC/Lymphocyte,group=Visit),col=2,size=1,alpha=0.5)+
  ylab("ANC/Lymphocyte")
dev.off()
write.csv(plt$data, "data/data_SupFig2C.csv", row.names = FALSE)
dfnow$ANC_Lymphocyte<-dfnow$ANC/dfnow$Lymphocyte
p<-wilcox.test(dfnow$ANC_Lymphocyte[dfnow$ARM == 'ARM A'],
               dfnow$ANC_Lymphocyte[dfnow$ARM == 'ARM B'],  exact = FALSE, correct = TRUE)
print(paste0('SubFig2C', ' pvalue=', p$p.value))

# ECOG Table
mean(df$ECOG[(df$First=="IO" & df$ARM=="ARM A")])
mean(df$ECOG[(df$Visit=="C1D1 - IO" & df$ARM=="ARM B")])


# KM plots ctDNA ---------------------------------------------

bsl<-orig.df[orig.df$Days2==0 & is.na(orig.df$`Percentage Mutant VAF`)==F,]

# set options for plotting survival curves
font.lab     <- 14
font.title   <- 16
font.tick    <- 12
conf.int     <- T
risk.table   <- T
ncensor.plot <- T
palette      <- c('red','blue')
ggtheme      <- theme_bw()
risk.table.height   <- 0.2
ncensor.plot.height <- 0.15

# bsl<-bsl%>%
#   filter(response_code != "77|Other specify")
bsl$Group<-"VAF<5%"
# censoring patient 7 due to non-cancer related death
# bsl$DEATH[bsl$PID==7] <- 0
bsl$Group[bsl$`Percentage Mutant VAF`>=5]<-"VAF>=5%"

# Figure Sup 3A
extra_patients<-read_csv("temp_data/CAcTUS screening patients data overall.csv") # I cannot convert it
extra_patients<- extra_patients %>%
  select(c(`screening number`,"ctDNA level (VAF)","Died no=0  yes=1","Date of death or last follow up","date starting treatment" ))
names(extra_patients)<-c('PID', 'VAF_pct','DEATH', 'date_of_death_censor', 'baseline_date')
extra_patients$baseline_date <- as.Date(extra_patients$baseline_date, format = "%d/%m/%Y")
extra_patients$date_of_death_censor <- as.Date(extra_patients$date_of_death_censor, format = "%d/%m/%Y")
extra_patients$OS<-extra_patients$date_of_death_censor - extra_patients$baseline_date 
# extra_patients$DEATH[extra_patients$PID=="01-909"] <- 0
# extra_patients$DEATH[extra_patients$PID=="01-901"] <- 0
extra_patients$OS <- as.integer(extra_patients$OS)
extra_patients$Group<-"VAF<5%"
extra_patients$Group[extra_patients$`VAF_pct`>=5]<-"VAF>=5%"
write_csv(extra_patients,'data/extra_patients_KM.csv')

# concatenate both
# bsl<-bsl %>%
#   select(names(extra_patients))
bsl$PID <- as.character(bsl$PID)
bsl2 <- bind_rows(bsl, extra_patients)
p.value <- survminer::surv_pvalue(fit = survfit(Surv(OS, DEATH) ~ Group, data=bsl2), data = bsl2)

# plot OS first
ggsurv<-survminer::ggsurvplot(
  survfit(Surv(OS, DEATH) ~ Group, data=bsl2),
  conf.int = conf.int,
  pval = round(p.value$pval,3),
  risk.table = risk.table,
  ncensor.plot = ncensor.plot,
  palette = palette,
  risk.table.col = "strata",
  legend.labs =
    c("VAF<5%","VAF>=5%"),
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
    plot.margin = unit(c(10, 10, 10, 64), "points"))
ggsurv$table <- ggsurv$table + 
  theme(
    plot.margin = unit(c(10, 10, 10, 45), "points"))
ggsurv$ncensor.plot <- ggsurv$ncensor.plot + 
  theme(
    plot.margin = unit(c(10, 10, 10, 90), "points"))
ggsurv


