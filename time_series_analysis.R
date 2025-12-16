require(ggplot2)
require(survival)
require(nlme)
library(tidyverse)


df<-read.csv("important_biomarkers_all2.csv")
data<-df

data$drdata$drdata$drug <- as.character(data$drug) # Ensure 'drug' is character type for manipulation
data$drug[data$drug == ""] <- NA

# Use group_by and fill to carry forward the last non-NA value within each patient group
data <- data %>%
  group_by(PID) %>%
  fill(drug, .direction = "down") %>%
  ungroup() 

# For the first entry of each patient where it should remain NA if it was NA,
# we reset it back to NA if it was the first occurrence and originally NA.
data <- data %>%
  group_by(PID) %>%
  mutate(drug = if_else(row_number() == 1 & is.na(drug[1]), NA_character_, drug)) %>%
  ungroup()

# Ensure the 'date' column is in Date format
data$date <- as.Date(data$date, format = "%Y-%m-%d")

# Calculate the baseline date for each patient
# The baseline date is the first non-NA date for each patient
data <- data %>%
  group_by(PID) %>%
  mutate(baseline_date = min(date[!is.na(drug)])) %>%
  ungroup()

# Calculate the difference in days from the baseline date for each entry
data <- data %>%
  mutate(days_from_baseline = as.integer(date - baseline_date))

# Convert variables to numeric (if it's not already)
data$Percentage_Mutant_VAF <- as.numeric(as.character(data$Percentage_Mutant_VAF))
data$ldh <- as.numeric(as.character(data$ldh))
data$mm_SLD <- as.numeric(as.character(data$mm_SLD))

# lets change drug to targeted or IO

data$drug[is.na(data$drug)==F & c(data$drug=="Dab+Tram"|data$drug=="Dab"|
                                    data$drug=="Tram"|data$drug=="Dab+Niv+Tram")]<-"Targeted"
data$drug[is.na(data$drug)==F & c(data$drug=="Ipi+Niv"|data$drug=="Niv"|
                                    data$drug=="Ipi")]<-"IO"



# Function to replace NA values at time 0 according to specified rules
replace_na_at_baseline <- function(df, var_name) {
  df %>%
    arrange(PID, date) %>%
    group_by(PID) %>%
    mutate(
      # Identify rows with NA at time 0
      is_na_at_0 = is.na(.data[[var_name]]) & days_from_baseline == 0,
      # Find the value before time 0 (if any)
      value_before_0 = ifelse(days_from_baseline < 0, .data[[var_name]], NA_real_),
      # Carry the last value forward, stopping at time 0
      value_before_0_filled = zoo::na.locf(value_before_0, na.rm = FALSE),
      # Find the first on-treatment value within 7 days after time 0 (if any)
      value_within_7_days = ifelse(days_from_baseline >= 0 & days_from_baseline <= 7, .data[[var_name]], NA_real_),
      # Carry the first value backward, starting from time 0
      value_within_7_days_filled = zoo::na.locf(value_within_7_days, fromLast = TRUE, na.rm = FALSE),
      # Determine the replacement value according to the rules
      replacement = coalesce(value_before_0_filled, value_within_7_days_filled),
      # Replace NA at time 0 with the determined replacement value
      !!sym(var_name) := ifelse(is_na_at_0, replacement, .data[[var_name]])
    ) %>%
    dplyr::select(-is_na_at_0, -value_before_0, -value_before_0_filled, -value_within_7_days, -value_within_7_days_filled, -replacement) %>%
    ungroup()
}
# Apply the function to replace NA values at time 0 for 'Percentage_Mutant_VAF'
data <- replace_na_at_baseline(data, "Percentage_Mutant_VAF")
data <- replace_na_at_baseline(data, "ldh")
data <- replace_na_at_baseline(data, "mm_SLD")

data<-data[data$days_from_baseline>=0,]


# Data frame for line plotting (exclude NA Percentage_Mutant_VAF)
line_data <- data %>% 
  filter(!is.na(Percentage_Mutant_VAF))


# Step 1: Directly track treatment transitions and assign unique identifiers to each period
data <- data %>%
  arrange(PID, date) %>%
  group_by(PID) %>%
  mutate(
    treatment_change = drug != lag(drug, default = drug[1]) | row_number() == 1,
    period_id = cumsum(treatment_change)
  ) %>%
  ungroup()

# Generate a unique identifier for each treatment period, considering the treatment name
data <- data %>%
  group_by(PID, period_id) %>%
  mutate(period_label = paste(drug, period_id, sep = "_")) %>%
  ungroup()

# Step 2: Calculate xmin and xmax for each unique treatment period
adjusted_data_for_shading <- data %>%
  group_by(PID, period_label) %>%
  summarise(
    drug = first(drug),
    xmin = min(days_from_baseline),
    xmax = max(days_from_baseline),
    .groups = 'drop'
  ) %>%
  arrange(PID, xmin)

# Adjust xmax to ensure non-overlapping periods
adjusted_data_for_shading <- adjusted_data_for_shading %>%
  group_by(PID) %>%
  mutate(xmax = if_else(row_number() < n(), lead(xmin) - 1, xmax)) %>%
  ungroup()

# Remove NA regions
adjusted_data_for_shading <- adjusted_data_for_shading %>%
  filter(!is.na(drug))

# Step 3: Plot with accurate shading for each treatment period
ggplot() +
  geom_line(data = line_data[line_data$days_from_baseline>=0,], aes(x = days_from_baseline, y = Percentage_Mutant_VAF, group = PID), na.rm = TRUE) +
  geom_point(data = data[data$days_from_baseline>=0,], aes(x = days_from_baseline, y = Percentage_Mutant_VAF), na.rm = TRUE) +
  facet_wrap(~ PID, scales = "free") +
  geom_rect(data = adjusted_data_for_shading, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = drug), alpha = 0.2) +
  scale_fill_manual(values = c("Targeted" = "blue", "IO" = "red")) +
  labs(title = "Percentage Mutant VAF Over Time by Patient",
       x = "Days from Treatment Start",
       y = "Mutant VAF [%]")


line_data <- data %>% 
  filter(!is.na(ldh))

ggplot() +
  geom_line(data = line_data, aes(x = days_from_baseline, y = ldh, group = PID), na.rm = TRUE) +
  geom_point(data = data, aes(x = days_from_baseline, y = ldh), na.rm = TRUE) +
  facet_wrap(~ PID, scales = "free") +
  geom_rect(data = adjusted_data_for_shading, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = drug), alpha = 0.2) +
  scale_fill_manual(values = c("Targeted" = "blue", "IO" = "red")) +
  labs(title = "LDH Over Time by Patient",
       x = "Days from Treatment Start",
       y = "LDH [IU/L]")


line_data <- data %>% 
  filter(!is.na(mm_SLD))

ggplot() +
  geom_line(data = line_data, aes(x = days_from_baseline, y = mm_SLD, group = PID), na.rm = TRUE) +
  geom_point(data = data, aes(x = days_from_baseline, y = mm_SLD), na.rm = TRUE) +
  facet_wrap(~ PID, scales = "free_x") +
  geom_rect(data = adjusted_data_for_shading, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = drug), alpha = 0.2) +
  scale_fill_manual(values = c("Targeted" = "blue", "IO" = "red")) +
  labs(title = "SLD Over Time by Patient",
       x = "Days from Treatment Start",
       y = "SLD [mm]")


# now we want to select ARM B
d1<-read.csv("../Data/CSV_export/CAcTUS.Table 1.csv")
colnames(d1)<-d1[1,]
d1<-d1[-1,]
colnames(d1)

data<-merge(data,d1[,c("PID","Treatment arm")],all.x=T,all.y=F)
dummy<-data[data$`Treatment arm`=="02|ARM B",]
dummy<-dummy[dummy$drug=="IO" & dummy$period_id==2,]

dummy <- dummy %>%
  group_by(PID) %>%
  mutate(baseline_date = min(date[!is.na(drug)])) %>%
  ungroup()
dummy <- dummy %>%
  mutate(days_from_baseline = as.integer(date - baseline_date))

d6<-read.csv("../Data/CSV_export/CAcTUS.Table 6_RECIST.csv")
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

dummy$date<-as.character(dummy$date)
d6$Date<-as.character(d6$Date)
colnames(d6)[2]<-"date"

dummy<-merge(dummy,d6,by=c("PID","date"),all.x=T,all.y=F)

dummy2<-unique(dummy[,c("PID","Overall response","period_label","Treatment arm")])
dummy[,c("PID","Overall response")]

dummy <- dummy %>%
  group_by(PID) %>%
  fill(`Overall response`, .direction = "updown") %>%
  ungroup()

ggplot(dummy[is.na(dummy$Percentage_Mutant_VAF)==F,],aes(days_from_baseline,
                                                         Percentage_Mutant_VAF,group=PID,
                 col=`Overall response`))+geom_line()+facet_wrap(~PID,scales="free_y")+
  xlab("Days from end of Targeted Therapy")+
  ylab("Mutant VAF [%]")+theme_bw(base_size=14)+geom_point()


library(dplyr)

data2<-data
# Ensure mm_SLD is numeric (if not already)
data2$mm_SLD <- as.numeric(data2$mm_SLD)

# Calculate % change from baseline for each patient
data2 <- data2 %>%
  group_by(PID, `Treatment arm`) %>%
  mutate(baseline_SLD = first(mm_SLD[!is.na(mm_SLD)]), # Get the first non-NA mm_SLD as baseline
         percent_change_SLD = (mm_SLD - baseline_SLD) / baseline_SLD * 100) %>%
  ungroup()
# Identify the best % change for each patient in each treatment arm
best_changes <- data2 %>%
  group_by(PID, `Treatment arm`) %>%
  summarize(best_percent_change = min(percent_change_SLD, na.rm = TRUE), .groups = 'drop') # Use min() for best reduction

library(ggplot2)

# Waterfall plot for best % change in mm_SLD, stratified by treatment arm
ggplot(best_changes, aes(x = reorder(PID, best_percent_change), y = best_percent_change, fill = `Treatment arm`)) +
  geom_col() +
  labs(x = "Patient ID",
       y = "Best % Change in SLD") +
  facet_wrap(~`Treatment arm`, scales = "free_x") +
  theme_bw(base_size=14)+
  scale_fill_manual(values = c("01|ARM A" = "orange", "02|ARM B" = "#558DAE"))
library(dplyr)
data3<-data
data3 <- data3 %>%
  mutate(mm_SLD = as.numeric(mm_SLD),  # Ensuring mm_SLD is numeric
         baseline_SLD = ifelse(days_from_baseline == 0, mm_SLD, NA),  # Identifying baseline SLD
         baseline_SLD = zoo::na.locf(baseline_SLD),  # Forward-filling baseline SLD values
         percent_change_SLD = (mm_SLD - baseline_SLD) / baseline_SLD * 100)  # Calculating percent change
library(lme4)

# Assuming 'percent_change_SLD' has been correctly calculated and exists in 'data3'
# Fit the linear mixed-effects model
data3$`Treatment arm`
model0 <- lmer(percent_change_SLD ~ days_from_baseline + (1|PID), data = data3, na.action = na.exclude)
model1 <- lmer(percent_change_SLD ~ `Treatment arm`+days_from_baseline + (1|PID), data = data3, na.action = na.exclude)
anova(model1,model0)#0.329

# Summary of the model
summary(model)
confint(model)

###############################################################################

require(lme4)
require(splines)
require(performance)
library(lme4)
library(merTools)

# for the mixed-effecst model with VAF and LDH as predictors we will
# use last value carried forward

data <- data %>%
  group_by(PID) %>%
  arrange(PID, date) %>%
  fill(Percentage_Mutant_VAF, ldh, .direction = "down") %>%
  ungroup()

summary(data$ldh)
# Model with time as a fixed effect and random intercepts and slopes for time
model <- lmer(mm_SLD ~ log(ldh) + days_from_baseline + 
                            (1 | PID), data = data, REML = FALSE)

# Display the summary of the model
summary(model)
confint.merMod(model)
# Calculate R-squared values
r_squared <- r2(model)
print(r_squared)

# lets do the prediction
newdata <- expand.grid(
  ldh = seq(min(data$ldh, na.rm = TRUE), max(data$ldh, na.rm = TRUE), length.out = 100),
  days_from_baseline = median(data$days_from_baseline, na.rm = TRUE), # Example: Using median days
  PID = factor(1) # Dummy PID since we're focusing on fixed effects
)

# Generate prediction intervals
pred_intervals <- predictInterval(model, newdata = newdata, level = 0.95, n.sims = 10000, include.resid = TRUE)

# Add ldh to the prediction intervals for plotting
pred_intervals$ldh <- newdata$ldh

ggplot(pred_intervals, aes(x = ldh)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) + # Prediction interval
  geom_line(aes(y = fit), color = "blue") + # Predicted values
  labs(x = "LDH", y = "Predicted mm_SLD", title = "Predicted mm_SLD with Prediction Intervals as a Function of LDH") +
  theme_minimal()




# Model with time as a fixed effect and random intercepts and slopes for time
model <- lmer(mm_SLD ~ log(Percentage_Mutant_VAF) + days_from_baseline + 
                (1 | PID), data = data, REML = FALSE)
# Display the summary of the model
summary(model)
confint.merMod(model)
# Calculate R-squared values
r_squared <- r2(model)
print(r_squared)

# Example: Generating a data frame for prediction
# Example: Predicting mm_SLD for a range of ldh values at a specific time point
newdata <- expand.grid(
  Percentage_Mutant_VAF = seq(min(data$Percentage_Mutant_VAF, na.rm = TRUE), max(data$Percentage_Mutant_VAF, na.rm = TRUE), length.out = 100),
  days_from_baseline = 42, # Example: Using median days
  PID = factor(1) # Dummy PID since we're focusing on fixed effects
)

# Generate predictions
# Generate prediction intervals
pred_intervals <- predictInterval(model, newdata = newdata, level = 0.95, n.sims = 10000, include.resid = F)

# Add ldh to the prediction intervals for plotting
pred_intervals$Percentage_Mutant_VAF <- newdata$Percentage_Mutant_VAF

ggplot(pred_intervals, aes(x = Percentage_Mutant_VAF)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) + # Prediction interval
  geom_line(aes(y = fit), color = "blue") + # Predicted values
  labs(x = "% VAF", y = "Predicted mm_SLD", title = "Predicted mm_SLD with Prediction Intervals as a Function of LDH") +
  theme_minimal()+scale_x_log10()



dummy<-data[is.na(data$mm_SLD)==F,]
model0 <- lmer(mm_SLD ~ days_from_baseline + 
                 (1 | PID), data = dummy, REML = FALSE)
model1 <- lmer(mm_SLD ~ log(ldh) + days_from_baseline + 
                (1 | PID), data = dummy, REML = FALSE)
model2 <- lmer(mm_SLD ~ log(Percentage_Mutant_VAF) + days_from_baseline + 
                 (1 | PID), data = dummy, REML = FALSE)
model3 <- lmer(mm_SLD ~ log(ldh) + log(Percentage_Mutant_VAF) + days_from_baseline + 
                 (1 | PID), data = dummy, REML = FALSE)
model4 <- lmer(mm_SLD ~ log(ldh)*log(Percentage_Mutant_VAF)*log(days_from_baseline+1) + 
                 (1 | PID), data = dummy, REML = FALSE)
summary(model4)
anova(model0,model1,model2,model3,model4)
r2(model0)
r2(model1)
r2(model2)
r2(model3)
summary(model0)
summary(model1)
summary(model2)
summary(model3)

residuals <- resid(model0)

# Plotting residuals vs. fitted values
ggplot(dummy, aes(x = predict(model3), y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values", y = "Residuals", title = "Residuals vs. Fitted Values")

# Checking for normality of residuals
ggplot(dummy, aes(x = residuals)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.5, color = "black", fill = "white") +
  geom_density(alpha = .2, fill = "#FF6666") +
  labs(x = "Residuals", y = "Density", title = "Histogram of Residuals")

qqnorm(residuals)
qqline(residuals)


# Calculate the standard deviation of the residuals
std_dev_residuals <- sd(residuals)

# Standardize the residuals
dummy$std_residuals <- residuals / std_dev_residuals

# Identify potential outliers
potential_outliers <- subset(dummy, abs(std_residuals) > 2)

# Plot to visualize potential outliers
ggplot(dummy, aes(x = predict(model3), y = std_residuals)) +
  geom_point() +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Standardized Residuals", title = "Standardized Residuals vs. Fitted Values")

