pacman::p_load(
        here,
        tidyverse,
        survival,
        survminer,
        janitor
)

heart_failure <- read_csv(here("data", "heart_failure.csv")) |> 
        mutate(
                # Convert NA to Unknown group
                ethnicgroup = factor(case_when(
                        is.na(ethnicgroup) ~ 8, 
                        .default = ethnicgroup)),
                ethnicgroup = fct_recode(ethnicgroup,
                                         "White" = "1", "Black" = "2",
                                         "Indian subcontinent" = "3",
                                         "Not known" = "8",
                                         "Other" = "9"
                ),
                gender = fct_recode(factor(gender),
                                    Male = "1", Female = "2"),
                cognitive = case_when(
                        senile == 1 | dementia == 1 ~ 1,
                        .default = 0
                ),
                quintile = fct_recode(factor(quintile),
                                      "1 (most affluent)" = "1",
                                      "5 (least affluent)" = "5"
                ),
                quintile = fct_relevel(quintile, "1 (most affluent)"),
                quintile_5groups = factor(fct_recode(quintile, "5 (least affluent)" = "0")),
                quintile_5groups = case_when(quintile == "0" ~ NA,
                                             .default = quintile)
        ) |> 
        droplevels()

heart_failure |> UEPtools::metadata_generator_any()

# The data are simulated based on real hospital administrative data for England called Hospital Episodes Statistics. Every public (National Health Service, NHS) hospital in the country must submit records for every admission; private hospitals also submit records for any NHS patients that they treat. The other UK countries and Ireland have similar databases. These can be linked to the national death registry in order to captures deaths that occur after discharge. 
# 
# Your simulated extract contains a random sample of emergency (unplanned) admissions for heart failure (ICD10 code I50). Here's a list of the fields and an explanation for some of them. Many of the fields are comorbidities coded as 0/1, where 1 indicates that the patient had it recorded. All comorbidities are recorded in HES's secondary diagnosis fields, of which there are currently 19. There are 24 fields to capture procedures and operations.
# 
# death (0/1)
# 
# los (hospital length of stay in nights)
# 
# age (in years)
# 
# gender (1=male, 2=female)
# 
# cancer
# 
# cabg (previous heart bypass)
# 
# crt (cardiac resynchronisation device - a treatment for heart failure)
# 
# defib (defibrillator implanted)
# 
# dementia
# 
# diabetes (any type)
# 
# hypertension
# 
# ihd (ischaemic heart disease)
# 
# mental_health (any mental illness)
# 
# arrhythmias
# 
# copd (chronic obstructive lung disease)
# 
# obesity
# 
# pvd (peripheral vascular disease)
# 
# renal_disease
# 
# valvular_disease (disease of the heart valves)
# 
# metastatic_cancer
# 
# pacemaker
# 
# pneumonia
# 
# prior_appts_attended (number of outpatient appointments attended in the previous year)
# 
# prior_dnas (number of outpatient appointments missed in the previous year)
# 
# pci (percutaneous coronary intervention)
# 
# stroke (history of stroke)
# 
# senile
# 
# quintile (socio-economic status for patient's neighbourhood, from 1 (most affluent) to 5 (poorest))
# 
# ethnicgroup (see below for categories)
# 
# fu_time (follow-up time, i.e. time in days since admission to hospital) 
# 
# Ethnic group has the following categories in this extract:
# 
# 1=white 
# 
# 2=black 
# 
# 3=Indian subcontinent 
# 
# 8=not known 
# 
# 9=other

cox <- coxph(Surv(fu_time, death)~ ethnicgroup, data = heart_failure)
summary(cox)
heart_failure |> summary()
table(heart_failure$gender, exclude=NULL)
table(heart_failure$prior_dnas, exclude=NULL)
table(heart_failure$ethnicgroup, exclude=NULL)
table(heart_failure$copd, exclude=NULL)
heart_failure |>  janitor::tabyl(gender)

cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + quintile_5groups + ethnicgroup, data = heart_failure)

summary(cox)

cox <- coxph(Surv(fu_time, death) ~ gender, data = heart_failure)
cox_test <- cox.zph(cox)
cox_test
plot(cox)
fit <- survfit(Surv(fu_time, death) ~ gender, data = heart_failure)
plot(fit)


cox <- coxph(Surv(fu_time, death) ~ copd, data = heart_failure)
cox_test <- cox.zph(cox)
cox_test
plot(cox_test)

names(heart_failure)

cox <- coxph(Surv(fu_time, death) ~ age + gender + ethnicgroup + ihd + stroke + arrhythmias + valvular_disease + pvd + copd + pneumonia + renal_disease + hypertension + cancer + mental_health + cognitive + prior_appts_attended + prior_dnas + los + defib + cabg, 
             data = heart_failure)
summary(cox)

cox2 <- coxph(Surv(fu_time, death) ~ age + gender + ethnicgroup + ihd + arrhythmias + valvular_disease + copd + pneumonia + renal_disease + cancer + mental_health + cognitive + prior_appts_attended + prior_dnas + los + defib + cabg, 
             data = heart_failure)
summary(cox2)
survminer::ggforest(cox, data = as.data.frame(heart_failure))
survminer::ggforest(cox2, data = as.data.frame(heart_failure))


cox_test <- cox.zph(cox)
cox_test
survminer::ggforest(cox, data = as.data.frame(heart_failure))

summary(coxph(Surv(fu_time, death) ~ gender, 
              data = heart_failure))
summary(coxph(Surv(fu_time, death) ~ gender + tt(gender), 
                   data = heart_failure))

survminer::ggsurvplot(
        fit = survfit(Surv(fu_time, death) ~ ihd,
                      data = heart_failure),
        conf.int = TRUE,
        risk.table = TRUE,
        break.time.by = 100            # present the time axis with an increment of 10 days
)

survminer::ggsurvplot(
        fit = survfit(Surv(fu_time, death) ~ hypertension,
                      data = heart_failure),
        risk.table = TRUE,
        break.time.by = 100            # present the time axis with an increment of 10 days
)

survminer::ggsurvplot(
        fit = survfit(Surv(fu_time, death) ~ cancer,
                      data = heart_failure),
        risk.table = TRUE,
        break.time.by = 100            # present the time axis with an increment of 10 days
)

survminer::ggsurvplot(
        fit = survfit(Surv(fu_time, death) ~ ethnicgroup,
                      data = heart_failure),
        risk.table = TRUE,
        break.time.by = 100            # present the time axis with an increment of 10 days
)

survminer::ggsurvplot(
        fit = survfit(Surv(fu_time, death) ~ stroke,
                      data = heart_failure),
        risk.table = TRUE,
        conf.int = TRUE,
        break.time.by = 100            # present the time axis with an increment of 10 days
)
