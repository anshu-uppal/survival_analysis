# Install and load relevant packages 

install.packages("survival") 

## Installing package into  
## (as 'lib' is unspecified) 

## package 'survival' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 
##  install.packages("ggplot2") 

## Installing package into  
## (as 'lib' is unspecified) 

## package 'ggplot2' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 
##  install.packages("survminer") 

## Installing package into  
## (as 'lib' is unspecified) 

## package 'survminer' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 
##   

install.packages("ggfortify") 

## Installing package into  
## (as 'lib' is unspecified) 

## package 'ggfortify' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 


require(survival)  

## Loading required package: survival 

## Warning: package 'survival' was built under R version 3.5.1 

require(ggplot2)  

## Loading required package: ggplot2 

## Warning: package 'ggplot2' was built under R version 3.5.1 

require(survminer) 

## Loading required package: survminer 

## Warning: package 'survminer' was built under R version 3.5.1 

## Loading required package: ggpubr 

## Warning: package 'ggpubr' was built under R version 3.5.1 

## Loading required package: magrittr 

require(ggfortify) 

## Loading required package: ggfortify 

## Warning: package 'ggfortify' was built under R version 3.5.1 

# Load dataset 
g <- read.csv(here("data", "heart_failure.csv"), header = TRUE, sep = ",")


# Define variables 
gender <- factor(g[,"gender"]) 
fu_time <- g[,"fu_time"]  
death <-  g[,"death"] 
age <- g[,"age"] 
copd <- factor(g[,"copd"]) 
ethnicgroup <- factor(g[,"ethnicgroup"]) 
quintile <- factor(g[,"quintile"]) 
ihd <- factor(g[,'ihd']) 
valvular <- factor(g[,'valvular_disease']) 
pvd <- factor(g[,'pvd']) 
stroke <- factor(g[,'stroke']) 
pneumonia <- factor(g[,'pneumonia']) 
renal <- factor(g[,'renal_disease']) 
ca <- factor(g[,'cancer']) 
mets <- factor(g[,'metastatic_cancer']) 
mental_health <- factor(g[,'mental_health']) 
ht <- factor(g[,"hypertension"]) 
cog_imp <- factor(g[,"senile"]) 
prior_dnas <- g[,"prior_dnas"] 



# Plotting a Kaplan-Meier curve 
###################### 

# 1. Generate the survival curve 
km_fit <- survfit(Surv(fu_time, death) ~ 1) 
# 2b. Alternative plot with ggplot2 
autoplot(km_fit) + theme_bw() # theme_bw() is a predesigned "theme" which makes the plot prettier 
# 1. Run the full model with all of your predictors  
cox <- coxph(Surv(fu_time, death) ~ age + gender + ethnicgroup + ihd +                 
                     valvular + pvd + stroke + copd + pneumonia + ht + renal +                 
                     ca + mets + mental_health + cog_imp + los + prior_dnas, data = g)  
summary(cox)


###################### 

# Output the probability of survival at certain times after hospital admission 
summary(km_fit, times = c(1:7,30,60,90*(1:10))) 

# Plotting a Kaplan-Meier curve by gender 
###################### 

# 1. Generate the survival curve 
km_gender_fit <- survfit(Surv(fu_time, death) ~ gender) 

# 2. Plot the curve 
plot(km_gender_fit)
# 2b. Alternative plot with ggplot2 
autoplot(km_gender_fit) + theme_bw() 

###################### 

# Perform log rank test to see whether survival varies by gender 
survdiff(Surv(fu_time, death) ~ gender, rho = 0) 

# Testing whether those over the age of 65 have different survival to those under it 
###################### 

# 1. Dichotomise age into categorical (binary in this case) variable 
age_65plus <- ifelse(g[,'age']>=65, 1, 0) 

# 2. Perform log rank test 
survdiff(Surv(fu_time, death) ~ age_65plus, rho = 0) 

###################### 

# Plot survival curve by age above or below 65 
###################### 

# 1. Generate survival curve 
km_old_fit <- survfit(Surv(fu_time, death) ~ age_65plus) 

# 2. Plot 
plot(km_old_fit) 
# 2b. Alternative plot in ggplot2 
autoplot(km_old_fit) + theme_bw()

###################### 

# Run Cox regression model with age as predictor (continuous variable) 
###################### 

# 1. Generate model 
cox <- coxph(Surv(fu_time, death) ~ age, data = g) 

# 2. Summarise model 
summary(cox) 
