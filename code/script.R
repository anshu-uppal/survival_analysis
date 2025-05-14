pacman::p_load(
        here,
        tidyverse,
        survival,
        janitor
)

linelist_case_data <- readRDS(here("data", "linelist_cleaned.rds"))


linelist_surv <-  linelist_case_data %>% 
        
        dplyr::filter(
                # remove observations with wrong or missing dates of onset or date of outcome
                date_outcome > date_onset) %>% 
        
        dplyr::mutate(
                # create the event var which is 1 if the patient died and 0 if he was right censored
                event = ifelse(is.na(outcome) | outcome == "Recover", 0, 1), 
                
                # create the var on the follow-up time in days
                futime = as.double(date_outcome - date_onset), 
                
                # create a new age category variable with only 3 strata levels
                age_cat_small = dplyr::case_when( 
                        age_years < 5  ~ "0-4",
                        age_years >= 5 & age_years < 20 ~ "5-19",
                        age_years >= 20   ~ "20+"),
                
                # previous step created age_cat_small var as character.
                # now convert it to factor and specify the levels.
                # Note that the NA values remain NA's and are not put in a level "unknown" for example,
                # since in the next analyses they have to be removed.
                age_cat_small = fct_relevel(age_cat_small, "0-4", "5-19", "20+")
        )

summary(linelist_surv$futime)
linelist_surv %>% janitor::tabyl(outcome, event)
linelist_surv %>% janitor::tabyl(age_cat_small, age_cat)

linelist_surv %>% 
        tabyl(gender, age_cat_small, show_na = F) %>% 
        adorn_totals(where = "both") %>% 
        adorn_percentages() %>% 
        adorn_pct_formatting() %>% 
        adorn_ns(position = "front")

# Surv-type object ####

# Use Suv() syntax for right-censored data
survobj <- Surv(time = linelist_surv$futime,
                event = linelist_surv$event)

# fit the KM estimates using a formula where the Surv object "survobj" is the response variable.
# "~ 1" signifies that we run the model for the overall survival  
linelistsurv_fit <-  survival::survfit(survobj ~ 1)

#print its summary for more details
summary(linelistsurv_fit)

#print its summary at specific times
summary(linelistsurv_fit, times = c(5,10,20,30,60))

# We can also use the print() function. The print.rmean = TRUE argument is used to obtain the mean survival time and its standard error (se).

# NOTE: The restricted mean survival time (RMST) is a specific survival measure more and more used in cancer survival analysis and which is often defined as the area under the survival curve, given we observe patients up to restricted time T (more details in Resources section).

# print linelistsurv_fit object with mean survival time and its se. 
print(linelistsurv_fit, print.rmean = TRUE)


# Plotting Kaplan-Meir curves ####
plot(linelistsurv_fit, 
     xlab = "Days of follow-up",    # x-axis label
     ylab="Survival Probability",   # y-axis label
     main= "Overall survival curve" # figure title
)

# original plot
plot(
        linelistsurv_fit,
        xlab = "Days of follow-up",       
        ylab = "Survival Probability",       
        mark.time = TRUE,              # mark events on the curve: a "+" is printed at every event
        conf.int = FALSE,              # do not plot the confidence interval
        main = "Overall survival curve and cumulative mortality"
)

# draw an additional curve to the previous plot
lines(
        linelistsurv_fit,
        lty = 3,             # use different line type for clarity
        fun = "event",       # draw the cumulative events instead of the survival 
        mark.time = FALSE,
        conf.int = FALSE
)

# add a legend to the plot
legend(
        "topright",                               # position of legend
        legend = c("Survival", "Cum. Mortality"), # legend text 
        lty = c(1, 3),                            # line types to use in the legend
        cex = .85,                                # parametes that defines size of legend text
        bty = "n"                                 # no box type to be drawn for the legend
)

# Comparison of survival curves ####
## By gender ####
linelistsurv_fit_sex <-  survfit(Surv(futime, event) ~ gender, data = linelist_surv)

survminer::ggsurvplot(
        linelistsurv_fit_sex, 
        data = linelist_surv,          # again specify the data used to fit linelistsurv_fit_sex 
        conf.int = FALSE,              # do not show confidence interval of KM estimates
        surv.scale = "percent",        # present probabilities in the y axis in %
        break.time.by = 10,            # present the time axis with an increment of 10 days
        xlab = "Follow-up days",
        ylab = "Survival Probability",
        pval = T,                      # print p-value of Log-rank test 
        pval.coord = c(40,.91),        # print p-value at these plot coordinates
        risk.table = T,                # print the risk table at bottom 
        legend.title = "Gender",       # legend characteristics
        legend.labs = c("Female","Male"),
        font.legend = 10, 
        palette = "Dark2",             # specify color palette 
        surv.median.line = "hv",       # draw horizontal and vertical lines to the median survivals
        ggtheme = theme_light()        # simplify plot background
)

## By infection source ####
linelistsurv_fit_source <-  survfit(
        Surv(futime, event) ~ source,
        data = linelist_surv
)

# plot
survminer::ggsurvplot( 
        linelistsurv_fit_source,
        data = linelist_surv,
        size = 1, linetype = "strata",   # line types
        conf.int = T,
        surv.scale = "percent",  
        break.time.by = 10, 
        xlab = "Follow-up days",
        ylab= "Survival Probability",
        pval = T,
        pval.coord = c(40,.91),
        risk.table = T,
        legend.title = "Source of \ninfection",
        legend.labs = c("Funeral", "Other"),
        font.legend = 10,
        palette = c("#E7B800","#3E606F"),
        surv.median.line = "hv", 
        ggtheme = theme_light()
)

# Cox regression analysis ####
#fitting the cox model
linelistsurv_cox_sexage <-  survival::coxph(
        Surv(futime, event) ~ gender + age_cat_small, 
        data = linelist_surv
)


#printing the model fitted
linelistsurv_cox_sexage
#summary of the model
summary(linelistsurv_cox_sexage)

## Model assumptions ####
# It was interesting to run the model and look at the results but a first look to verify whether the proportional hazards assumptions is respected could help saving time.
test_ph_sexage <- survival::cox.zph(linelistsurv_cox_sexage)
test_ph_sexage

## Expanded model ####
#fit the model
linelistsurv_cox <-  coxph(
        Surv(futime, event) ~ gender + age_years+ source + days_onset_hosp,
        data = linelist_surv
)


#test the proportional hazard model
linelistsurv_ph_test <- cox.zph(linelistsurv_cox)
linelistsurv_ph_test
# graphical verification of this assumption:
survminer::ggcoxzph(linelistsurv_ph_test)

# Forest plots ####
survminer::ggforest(linelistsurv_cox, data = linelist_surv)

