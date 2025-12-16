# install.packages(c("survival", "survminer"))
# https://www.r-bloggers.com/2016/12/cox-proportional-hazards-model/
library("survival")
library("survminer")

#### Longitudinal COX PH
df1 <- read.csv('coxph_df1.csv')
df2 <- read.csv('important_biomarkers_all_interp2.csv')

# transform the data
df1$dead_or_PD_date <- as.Date(df1$death_date, format = "%Y-%m-%d")
df1$treatment_start_date <- as.Date(df1$treatment_start_date, format = "%Y-%m-%d")

# re-arrange dataset
tdata <- with(df1, data.frame(PID = PID,
                               pfstime = pmax(.5, dead_or_PD_date - treatment_start_date),
                               pfs = dead_or_PD
))

tdata2 <- tmerge(tdata, tdata, id=PID, 
                pfs = event(pfstime, pfs)) #set range

# this data is for the whole dataset
sdata <- tmerge(tdata2, df2, id=PID, 
                mm_SLD = tdc(days_from_baseline, mm_SLD_interp),
                ldh = tdc(days_from_baseline, LDH_interp), 
                Percentage_Mutant_VAF = tdc(days_from_baseline, Percentage_Mutant_VAF))

# this fits the Cox PH overall
# we probably want to log-transform these variables as they are skewed
# why is that important? well the cox PH model is a log-linear model and it assumes a
# unit increase in risk is asscoaited witha unit increase in teh predcitor

# we can't have 0s when we log transform 

fit<-coxph(Surv(tstart, tstop, pfs) ~ log(mm_SLD+5) + log(ldh) + log(Percentage_Mutant_VAF+0.01),
      data= sdata)

coef(fit)
summary(fit)

# 

fit<-coxph(Surv(tstart, tstop, pfs) ~ log(mm_SLD+5) ,
           data= sdata)
coef(fit)
summary(fit)

fit<-coxph(Surv(tstart, tstop, pfs) ~ log(ldh) ,
           data= sdata)
coef(fit)
summary(fit)

fit<-coxph(Surv(tstart, tstop, pfs) ~ log(Percentage_Mutant_VAF+0.01) ,
           data= sdata)
coef(fit)
summary(fit)

# so whats interesting is that we know the study is under-powered
# SLD should be a predictor and we signs of it
# If we look at the conodrance that tells us there is a signal in both SLD
# and VAF

fit<-coxph(Surv(tstart, tstop, pfs) ~ log(Percentage_Mutant_VAF+0.01) +
             log(mm_SLD+5),
           data= sdata)
coef(fit)
summary(fit)

# whats interesting is that combining the two doesn't improve anything
# if you look at the regression coefficients they both drop in value
# suggesting they both capture similar information





# this looks at each covariate independently
covariates <- c("Percentage_Mutant_VAF", "ldh","mm_SLD")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(tstart, tstop, pfs)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = sdata)})
# Extract results
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

# now only of the first part
df2$first_descent %<>% as.logical(df2$first_descent)
df2_red <- subset(df2, first_descent)
sdata <- tmerge(tdata2, df2_red, id=PID, 
                mm_SLD = tdc(days_from_baseline, mm_SLD_interp),
                ldh = tdc(days_from_baseline, LDH_interp), 
                Percentage_Mutant_VAF = tdc(days_from_baseline, Percentage_Mutant_VAF))

fit2<-coxph(Surv(tstart, tstop, pfs) ~ mm_SLD + ldh + Percentage_Mutant_VAF,
           data= sdata)

# this fits the Cox PH for the first part of the curve only
coef(fit2)
summary(fit2)

# this looks at each covariate independently
covariates <- c("Percentage_Mutant_VAF", "ldh")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(tstart, tstop, pfs)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = sdata)})
# Extract results
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)


##### cross-sectional COX PH
df <- read.csv('/Users/juandelgado/Desktop/Juan/code/systems_forecasting/cactus/data/processed_data/coxph.csv')

# multivariate analysis
covariates <- c("Percentage_Mutant_VAF", "ldh")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(days_from_eligibility, PFS)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = df)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
