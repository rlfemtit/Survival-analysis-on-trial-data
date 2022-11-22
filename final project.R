

library(survMisc)
library(survival)
library(flexsurv)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(survsim)
library(dplyr)

data1 = read.csv("phase2.csv")
data2 = read.csv("phase3.csv")

data1["time"] = as.numeric(difftime(data1$end, data1$start, units="days", tz="GMT")) # survived time
data1["medical_status"] = rep(0, nrow(data1))

for (i in 1:nrow(data1)){
    if (data1["karno"][i,] %in% c(10, 20, 30) == TRUE){
        data1["medical_status"][i,] = "completely hospitalized"
    }
    else if (data1["karno"][i,] %in% c(40, 50, 60) == TRUE){
        data1["medical_status"][i,] = "partial confinement to a hospital"    
    }
        else{ 
                data1["medical_status"][i,] = "able to care for self"
        }
    
}

# Assumption for KM and log-rank : censoring is unrelated to time-to-event.

### Part1
# 1
# descriptive statistics for celltype/karnoage/age

data1_t1 = data1[data1["trt"] == 1,] # standard
data1_t2 = data1[data1["trt"] == 2,] # new


# celltype
other_table = table(data1$trt, data1$celltype)
barplot(other_table,
        xlab = "Cell type", ylab = "Frequency",
        col = c("darkgrey", "black"),
        legend.text = c("standard", "new"),
        beside = TRUE) # Grouped bars

# karno
plot(density(data1_t1$karno, xlim=c(0,90)), col="blue", xlab="Medical Status", main="Density funciton of medical status")
lines(density(data1_t2$karno), col='red', add=T)
legend(x = "topleft",          # Position
       legend = c("standard", "new"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c(2, 4),           # Line colors
       lwd = 2)                 # Line width

# age
summary(data1_t1$age)
summary(data1_t2$age)

plot(density(data1_t1$age), col="blue", main="Age Distribution", xlab="Age")
lines(density(data1_t2$age), col='red', add=T)
legend(x = "topleft",          # Position
       legend = c("standard", "new"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c(2, 4),           # Line colors
       lwd = 2)                 # Line width


# 2
data1$trt = ifelse(data1$trt==1,"standard","new")
data1
# statistics for the survival time : mean survival time / median survival time / survival curve / CI
km_fit <- survfit(Surv(time, status) ~ trt, data=data1)
autoplot(km_fit, conf.int = TRUE, censor = FALSE, main="Survival curve per the treatment group")


# standard : tr1
hist(data1$time, breaks=20) # outliers exist -> median would be of great interest for many stake holders
summary(km_fit)
print(km_fit, print.rmean=TRUE) # Survival probability in KM plot approximately reached to 0. So we can get the mean.
103 # median since S(t=103) = 0.4863 or we can use the linear interpolation to get 
# confidence interval for median : 59 ~ 132

# new : tr2
mean(data1_t2$time)
52 # median since S(t=52)=0.5
# confidence interval for median : 44 ~ 95

# 3
# hypothesis test if the two survival functions are identical.
# Since there are many crossing points between functions, use the Renyi type test
comp(ten(km_fit)) # p-value when weight=1 is 0.26013 -> can not reject H0 that two hazard functions are same.

# 4
# parametric model : expo, weibull, generalized gamma, log-logistic, log-normal
# no difference between treatment -> no need to seperate the data

# exponential
fitE <- flexsurvreg(formula = Surv(time, status) ~ 1, data = data1, dist="exp")
plot(fitE, conf.int=FALSE, main="Parametric survival model", 
     ylab="Survival probability", xlab="Time", ci=FALSE, lwd=1)

fitE

fitW <- flexsurvreg(formula = Surv(time, status) ~ 1, data = data1, dist="weibull")
lines(fitW, col="green", ci=FALSE, lwd=1)

fitG <- flexsurvreg(formula = Surv(time, status) ~ 1, data = data1, dist="gengamma")
lines(fitG, col="blue", ci=FALSE, lwd=1)

fitL <- flexsurvreg(formula = Surv(time, status) ~ 1, data = data1, dist="lnorm")
lines(fitL, col="purple", col.ci="purple", ci=FALSE, lwd=1)

legend(x = "topright",          # Position
       legend = c("Kaplan-Meier", "Exponential", "Weibull", "Generalized gamma", "log-normal"),  # Legend texts
       lty = c(1, 1, 1, 1, 1),             # Line types
       col = c(1, 2, 3, 4, "blueviolet"),  # Line colors
       lwd = 3)                            # Line width
# Generalized gamma and Weibull seem the best fit.

# Exp
fitg <- flexsurvreg(formula = Surv(time, status) ~ 1, data = data1, dist="exp")
plot(fitg, conf.int=FALSE, main="Exponential distribution", ylab="Survival probability", xlab="Time")
legend(x = "topright",          # Position
       legend = c("Exponential", "Kaplan-Meier"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c(2, 1),           # Line colors
       lwd = 2)                 # Line width

# Weibull
fitg <- flexsurvreg(formula = Surv(time, status) ~ 1, data = data1, dist="Weibull")
plot(fitg, conf.int=FALSE, main="Weibull distribution", ylab="Survival probability", xlab="Time")
legend(x = "topright",          # Position
       legend = c("Weibull", "Kaplan-Meier"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c(2, 1),           # Line colors
       lwd = 2)                 # Line width

# generalized gamma
fitg <- flexsurvreg(formula = Surv(time, status) ~ 1, data = data1, dist="gengamma")
fitg
plot(fitg, conf.int=FALSE, main="Generalized gamma", ylab="Survival probability", xlab="Time")
legend(x = "topright",          # Position
       legend = c("Generalized Gamma", "Kaplan-Meier"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c(2, 1),           # Line colors
       lwd = 2)                 # Line width

# log normal
fitg <- flexsurvreg(formula = Surv(time, status) ~ 1, data = data1, dist="lnorm")
plot(fitg, conf.int=FALSE, main="log-normal distribution", ylab="Survival probability", xlab="Time")
legend(x = "topright",          # Position
       legend = c("log-normal", "Kaplan-Meier"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c(2, 1),           # Line colors
       lwd = 2)                 # Line width

# log logistic
fitLL <- flexsurvreg(formula = Surv(time, status) ~ 1, data = data1, dist="llogis")
plot(fitLL, conf.int=FALSE, main="log-logistic distribution", ylab="Survival probability", xlab="Time")
legend(x = "topright",          # Position
       legend = c("log-logistic", "Kaplan-Meier"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c(2, 1),           # Line colors
       lwd = 2)                 # Line width

# Gompertz
fitGO <- flexsurvreg(formula = Surv(time, status) ~ 1, data = data1, dist="gompertz")
plot(fitGO, conf.int=FALSE, main="Gompertz distribution", ylab="Survival probability", xlab="Time")
legend(x = "topright",          # Position
       legend = c("Gompertz", "Kaplan-Meier"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c(2, 1),           # Line colors
       lwd = 2)                 # Line width

# Compare the parametric models using likelihood or AIC

fitE
fitW
fitG # the best according to AIC
fitL
 
1 - pchisq(2 * (-748.0912 + 751.2212), df = 1) # p-value



# Bonus : treatment effect on the different cell type or karno or age
# cell type
data1_squa = data1[data1$celltype == "squamous",]
data1_adeno = data1[data1["celltype"] == "adeno",]
data1_large = data1[data1["celltype"] == "large",]
data1_small = data1[data1["celltype"] == "smallcell",]

# stratified by celltype : significant different in squamous cell type.
km_fit <- survfit(Surv(time, status) ~ trt, data=data1_squa)
summary(km_fit)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Group of squamous cell type")
comp(ten(km_fit),p=0, q=1) # FM test -> 0.03(idenical rather than crossed)

km_fit <- survfit(Surv(time, status) ~ trt, data=data1_adeno)
summary(km_fit)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Group of adeno cell type")
comp(ten(km_fit),p=0, q=1) # Renyi type test w=1/FH -> 0.46/0.24

km_fit <- survfit(Surv(time, status) ~ trt, data=data1_large)
summary(km_fit)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Group of large cell type")
comp(ten(km_fit),p=1, q=0) # Renyi type test w=1/FH -> 0.27/0.1134

km_fit <- survfit(Surv(time, status) ~ trt, data=data1_small)
summary(km_fit)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Group of small cell type")
comp(ten(km_fit),p=0, q=1) # log rank/FM -> 0.1325/0.0582

# karno
data1_completely =data1[data1["medical_status"] == "completely hospitalized",]
data1_partial =data1[data1["medical_status"] == "partial confinement to a hospital",]
data1_able =data1[data1["medical_status"] == "able to care for self",]

# stratified by treament
km_fit <- survfit(Surv(time, status) ~ medical_status, data=data1_t1)
summary(km_fit)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Group of standard treatment ")

km_fit <- survfit(Surv(time, status) ~ medical_status, data=data1_t2)
summary(km_fit)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Group of new treatment")

# stratified by medical status
km_fit <- survfit(Surv(time, status) ~ trt, data1_completely)
summary(km_fit)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Group of completely hospitalized ")
comp(ten(km_fit)) # Reyni type test -> p=0.53

km_fit <- survfit(Surv(time, status) ~ trt, data=data1_partial)
summary(km_fit)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Group of partial confinement to a hospital")
comp(ten(km_fit), p=0, q=1) # Reyni type test + FH-> p=0.13

km_fit <- survfit(Surv(time, status) ~ trt, data=data1_able)
summary(km_fit)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Group of able to care for self")
comp(ten(km_fit), p=0, q=1) # Reyni type test + FH-> p=0.1



### Part2
## 1
set.seed(1523856)

# data generation 

# baseline hazard: Weibull
# N = sample size    
# lambda = scale parameter in h0()
# rho = shape parameter in h0()
# beta = fixed effect parameter
# followup = follow-up time
# source : https://stats.stackexchange.com/questions/135124/how-to-create-a-toy-survival-time-to-event-data-with-right-censoring

simulWeib <- function(N, lambda, rho, beta, followup){
    # assign treatment variable
    x <- sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5))
    
    # Weibull event times
    v <- runif(n=N)
    Tlat <- (- log(v) / (lambda * exp(x * beta)))^(1 / rho)
    
    # censoring times
    nr_random_censored = round(N * 0.05)
    C = vector(length=nr_random_censored)
    
    # random censoring time
    j=1
    for (i in Tlat[1:nr_random_censored]){
        if (i >= followup){
            C[j] = runif(1, min=0, max=followup)
        }else{
            C[j] = runif(1, min=0, max=i)
        }
        j=j+1
    }
    
    # right censoring time
    C <- append(C, rep(followup, N-nr_random_censored)) # if not dropped, C is the follow-up time.
    
    # follow-up times and event indicators
    time <- pmin(Tlat, C)
    status <- as.numeric(Tlat <= C)
    
    # data set
    data.frame(id=1:N,
               time=time,
               status=status,
               x=x)
}
data_sim = simulWeib(N=300, lambda=0.01, rho=1, beta= 0, followup=120)
data_sim

km_fit <- survfit(Surv(time, status) ~ x, data=data_sim)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Survival function per treatment")
result = survdiff(formula = Surv(time, status) ~ x, data = data_sim) # log-rank test
p = 1 - pchisq(result$chisq, length(result$n) - 1)
p

# log-rank test : type1 error, power, coverage prob
# factors to impact the power : HR, a, sample proportion(determined), sample size(determined)
log_rank_power <- function(N, followup, a, beta){
    # N : sample size 
    # followup : follow-up time 
    # a : confidence level  
    # beta : effect of treatment
    
    p_vector = vector(length=1000)
    
    for (i in 1:1000){
        data_sim = simulWeib(N, lambda=0.01, rho=1, beta= beta, followup) 
        result = survdiff(formula = Surv(time, status) ~ x, data = data_sim) # log-rank test
        p = 1 - pchisq(result$chisq, length(result$n) - 1) # p-value
        p_vector[i] = p
    }
    
    power_log_rank = mean(p_vector < a) # power when the confidence level is 0.05
    return(power_log_rank)
}
log_rank_power(N=300, followup=120, a=0.05, beta=-0.5)

# power per different sample size
x = seq(200, 400, 50)
y = vector(length=length(x))
j=1
for (i in x){
    y[j] = log_rank_power(N=i, followup=120, a=0.05, beta=-0.5)
    j = j+1
}
plot(x,y, type='o', xlab="Nr of patients", ylab="Power", main="Power of log-rank test"
     , ylim=c(0.7, 1))

# power per different followup
x = seq(90, 240, 30)
y = vector(length=length(x))
j=1
for (i in x){
    y[j] = log_rank_power(N=300, followup=i, a=0.05, beta=-0.5)
    j = j+1
}
plot(x,y, type='o', xlab="follow-up time", ylab="Power", main="Power of log-rank test"
     , ylim=c(0.8,1))

# power per beta
x = seq(0.1, 1, 0.2)
y = vector(length=length(x))
j=1
for (i in x){
    y[j] = log_rank_power(N=300, followup=120, a=0.05, beta=i)
    j=j+1
}
plot(x,y, type='o', xlab="Coefficient", ylab="Power", main="Power of log-rank test"
     , ylim=c(0,1))

# type 1 error can be calculated by setting beta to 0 : generate data under H0
log_rank_power(N=300, followup=120, a=0.05, beta= 0) # 0.048




# median survival time : coverage probability of 95% 
coverage_median <- function(N, followup, a, beta){
    # N : sample size 
    # followup : follow-up time 
    # a : confidence level  
    # beta : effect of treatment
    
    # True median value using the population data
    data_sim = simulWeib(N=100000, lambda=0.01, rho=1, beta= 0, followup=followup)
    km = survfit(Surv(time, status) ~ 1, data=data_sim)
    true_median = summary(km)$table["median"] 
    
    coverage_vector = vector(length=1000)
    
    for (i in 1:1000){
        data_sim = simulWeib(N, lambda=0.01, rho=1, beta= beta, followup)
        km_fit = survfit(Surv(time, status) ~ 1, data=data_sim)
        result = summary(km_fit)$table
        
        lower_bound = result[8]
        upper_bound = result[9]

        coverage_vector[i] = between(true_median, lower_bound, upper_bound)
    }
    
    coverage_rate = mean(coverage_vector) # power when the confidence level is 0.05
    return(coverage_rate)
}
coverage_median(N=300, followup=120, a=0.05, beta=0) # 0.944
coverage_median(N=300, followup=180, a=0.05, beta=0) # 0.955

# coverage per followup
x = seq(120, 210, 30)
y = vector(length=length(x))
j=1
for (i in x){
    y[j] = coverage_median(N=300, followup=i, a=0.05, beta=0)
    j=j+1
}
plot(x,y, type='o', xlab="Follow-up time", ylab="Coverage rate", 
     main="Coverage rate of median survival time", ylim=c(0.8,1), xlim=c(120,210))

### Part3
data2 = read.csv("phase3.csv")
km_fit <- survfit(Surv(time, death) ~ 1, data=data2)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Survival function")

hist(data2$time)
data2$age_categories <- cut(data2$age, breaks = 3, labels = c("young","middle","senior"))
data2$trt = ifelse(data2$treatment == 1, "new", "standard")
data2$cdp_categories <- cut(data2$cdp, breaks = 3, labels = c("low","medium","high"))

# Before constructing the model, do tests if the effects of covariates are significant.
# Non parametric test is more robust method than the parametric way.
# plot the survival curve per covariates

#age : not effective
km_fit <- survfit(Surv(time, death) ~ age_categories, data=data2)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Survival function per age group")
comp(ten(km_fit))
comp(ten(km_fit))$tests$trendTests # p=0.3154 weight=1 -> there is no trend

#trt : effective 
km_fit <- survfit(Surv(time, death) ~ trt, data=data2)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Survival function per treatment")
comp(ten(km_fit))# log-rank p=0.0002  

#gender : not effective
km_fit <- survfit(Surv(time, death) ~ gender, data=data2)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Survival function per gender")
comp(ten(km_fit)) # log-rank/renyi(weight=1) p=0.3107/0.5963

#smoke history : not effective
km_fit <- survfit(Surv(time, death) ~ smoke_history, data=data2)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Survival function per smoke history")
comp(ten(km_fit)) # log-rank/renyi(weight=1) p=0.5144/0.7238

#cell type : effective
km_fit <- survfit(Surv(time, death) ~ cell_type, data=data2)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Survival function per cell type")
comp(ten(km_fit)) # test for four samples : p-value(w=1) = close to 0

#cdp : effective
km_fit <- survfit(Surv(time, death) ~ cdp_categories, data=data2) # for visualization
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Survival function per cdp group")
km_fit <- survfit(Surv(time, death) ~ cdp, data=data2) # for test
comp(ten(km_fit))$tests$trendTests # p=0.0259 weight=1 -> there is trend

# Check the multi-collilnearity 


## cox model
# reference : female, non-smoker, small cell
data2$cell_type = relevel(data2$cell_type, ref = "smallcell")
data2$gender = relevel(data2$gender, ref = "female")
data2$smoke_history = relevel(data2$smoke_history, ref = "non-smoker")

# Assumptions
# 1. ln(hazard) is a linear function of x -> residual plot
# 2. proportional hazards  -> c log log or schoenfeld test
# 3. x is not time-series data

cox_model = coxph(Surv(time, death) ~ age + treatment  + 
                      cell_type + cdp, data=data2, ties=c("efron")) # Efron is better than Breslow

# 1. ln(hazard) is a linear function of x -> residual plot
plot(predict(cox_model), residuals(cox_model, type='martingale'), xlab="Fitted value"
     , ylab="Residual", main="Martingale residual")
abline(h=0, col="red")
lines(smooth.spline(predict(cox_model),
                    residuals(cox_model, type='martingale')), col="red")

# 2. proportional hazards  -> c log log or schoenfeld test
cox.zph(cox_model, terms=FALSE) 


# 3. Check time independence : every covariates are time-independent

data2[data2$death == 0,]
0.5651-1 
#  If the black solid line is fairly flat and straight, then proportionality is supported.
# Red line is the time-fixed effect estimated by a Cox model 
plot(cox.zph(cox_model, terms=FALSE)[1], main="Age") # treatment's effect changes over time
abline(h=summary(cox_model, terms=FALSE)$coefficient["age",][1], col="red")

plot(cox.zph(cox_model, terms=FALSE)[2], main="Treatment") # treatment's effect changes over time
abline(h=summary(cox_model, terms=FALSE)$coefficient["treatment",][1], col="red")

plot(cox.zph(cox_model, terms=FALSE)[3], main="Adeno", ylab="Beta(t) for Adeno") # Cell_type : time dependent
abline(h=summary(cox_model, terms=FALSE)$coefficient["cell_typeadeno",][1], col="red")

plot(cox.zph(cox_model, terms=FALSE)[4], main="Large", ylab="Beta(t) for Large") # Cell_type : time dependent
abline(h=summary(cox_model, terms=FALSE)$coefficient["cell_typelarge",][1], col="red")

plot(cox.zph(cox_model, terms=FALSE)[5], main="Squamous", ylab="Beta(t) for Squamous") # Cell_type : time dependent
abline(h=summary(cox_model, terms=FALSE)$coefficient["cell_typesquamous",][1], col="red")

plot(cox.zph(cox_model, terms=FALSE)[6], main="CDP") # CDP : not time dependent
abline(h=summary(cox_model, terms=FALSE)$coefficient["cdp",][1], col="red")


# LRT test
cox_model1 = coxph(Surv(time, death) ~ age + treatment  + 
                      cell_type + cdp, data=data2, ties=c("efron")) # Efron is better than Breslow
cox_model2 = coxph(Surv(time, death) ~ treatment  + 
                      cell_type + cdp, data=data2, ties=c("efron")) # Efron is better than Breslow
anova(cox_model1, cox_model2, test="LRT") # p=0.1362, 




# Goodness-of-fit test using cox-snell residual -> quite good
mart_residual = residuals(cox_model, type = "martingale")
cox_snell_residual = -(mart_residual - data2$death)
fit_coxsnell = coxph(Surv(cox_snell_residual, death) ~ 1, data=data2, ties=c("efron"))
df_base_haz <- basehaz(fit_coxsnell, centered = FALSE)

p <- ggplot(data = df_base_haz, mapping = aes(x = time, y = hazard)) +
    geom_point() +
    scale_x_continuous(limit = c(0,3.5)) +
    scale_y_continuous(limit = c(0,3.5)) +
    labs(x = "Cox-Snell residuals ",
         y = "H(r)") + 
    theme_bw() + theme(legend.key = element_blank())
p + geom_abline(slope=1, col="red", lwd=1)
abline(h=1)


# prediction
cox_model = coxph(Surv(time, death) ~ treatment  + 
                      cell_type + cdp, data=data2, ties=c("efron")) # Efron is better than Breslow
summary(cox_model)

# predicted survival for a 60 year old (3 month and 6 month)
newdata=data.frame(treatment=0,cell_type="squamous", cdp=4)
newdata2=data.frame(treatment=1,cell_type="squamous", cdp=4)
plot(survfit(cox_model, newdata=newdata, conf.int=F), xlab="Time(days)", ylab="Survival", 
     col="blue", main="Predicted survival probability")
lines(survfit(cox_model, newdata=newdata2, conf.int=F), xlab="Time(days)", ylab="Survival", col="Red")
legend(x = "topright",          # Position
       legend = c("new", "standard"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c(2, 4),           # Line colors
       lwd = 2)                 # Line width

summary(survfit(cox_model, newdata=newdata, conf.int=F)) # standard 0.0478
summary(survfit(cox_model, newdata=newdata2, conf.int=F)) # new 0.1793

data1

## Further discussion
# Check the different treatment effect per stata using stratified test
data2 = read.csv("phase3.csv")
data2$age_categories <- cut(data2$age, breaks = 2, labels = c("young","old")) # age<61 -> young
data2$trt = ifelse(data2$treatment == 1, "treatment", "untreated")
data2$cdp_categories <- cut(data2$cdp, breaks = 2, labels = c("low","high"))

# Again, when we create more specific strata, randomization might be gone.

# gender - treatment -> bigger effect in female
km_fit <- survfit(Surv(time, death) ~ gender + trt, data=data2)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Effect of treatment per gender")
comp(ten(km_fit),p=0, q=1) 

# smoke history - treatment -> bigger effect in smoker
km_fit <- survfit(Surv(time, death) ~ smoke_history + trt, data=data2)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Effect of treatment per gender")
comp(ten(km_fit),p=0, q=1) # later effect -> FM test(weight=1) -> 0.0364

# age - treatment -> bigger effect in the old
km_fit <- survfit(Surv(time, death) ~ age_categories + trt, data=data2)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Effect of treatment per gender")
comp(ten(km_fit),p=0, q=1) # later effect -> FM test(weight=1) -> 0.0364

# CDP - treatment
km_fit <- survfit(Surv(time, death) ~ cdp_categories + trt, data=data2)
autoplot(km_fit, conf.int = FALSE, censor=FALSE, main="Effect of treatment per gender")
comp(ten(km_fit),p=0, q=1) # later effect -> FM test(weight=1) -> 0.0364


# AFT model -> no need since the hazards are proportional.
















