# This R code produces Supplementary Table 1 results in the manuscript Gilbert, Huang et al. 

library(survival)
library(cmprsk)
library(dplyr)

# input data
surv_nab <- read.csv("../data/amp_survival_wk80_wk104.csv")
surv_nab$cov <- factor(surv_nab$cov, levels=c(1,0))

### Fine-Gray

# model 1: pooled VRC01 vs. placebo
dat1 <- surv_nab
model1 <- crr(dat1$fudays_new, dat1$status_new80ls, cov1= dat1$cov, failcode=1, cencode=0) 
print(summ1 <- summary(model1)$conf.int)

hr = summ1[, "exp(coef)"]
lo.hr = summ1[, "2.5%"]
up.hr = summ1[, "97.5%"]

print(pe <- 1 - hr)
print(lo.pe <- 1 - up.hr)
print(up.pe <- 1 - lo.hr)

# model 2: VRC01 10mg/kg vs. placebo
dat2 <- filter(surv_nab, rx_code %in% c("T1", "C3"))
model2 <- crr(dat2$fudays_new, dat2$status_new80ls, cov1= dat2$cov, failcode=1, cencode=0) 
print(summ2 <- summary(model2)$conf.int)

hr = summ2[, "exp(coef)"]
lo.hr = summ2[, "2.5%"]
up.hr = summ2[, "97.5%"]

print(pe <- 1 - hr)
print(lo.pe <- 1 - up.hr)
print(up.pe <- 1 - lo.hr)

# model 3: VRC01 30mg/kg vs. placebo
dat3 <- filter(surv_nab, rx_code %in% c("T2", "C3"))
model3 <- crr(dat3$fudays_new, dat3$status_new80ls, cov1= dat3$cov, failcode=1, cencode=0) 
print(summ3 <- summary(model3)$conf.int)

hr = summ3[, "exp(coef)"]
lo.hr = summ3[, "2.5%"]
up.hr = summ3[, "97.5%"]

print(pe <- 1 - hr)
print(lo.pe <- 1 - up.hr)
print(up.pe <- 1 - lo.hr)

### Cause-specific Cox

# model 1: pooled VRC01 vs. placebo
model1cox <- coxph(Surv(fudays_new, status_new80ls_cox) ~ cov, data=dat1)
print(summ1cox <- summary(model1cox)$conf.int)

hr = summ1cox[, "exp(coef)"]
lo.hr = summ1cox[, "lower .95"]
up.hr = summ1cox[, "upper .95"]

print(pe <- 1 - hr)
print(lo.pe <- 1 - up.hr)
print(up.pe <- 1 - lo.hr)
  
# model 2: VRC01 10mg/kg vs. placebo
model2cox <- coxph(Surv(fudays_new, status_new80ls_cox) ~ cov, data=dat2)
print(summ2cox <- summary(model2cox)$conf.int)

hr = summ2cox[, "exp(coef)"]
lo.hr = summ2cox[, "lower .95"]
up.hr = summ2cox[, "upper .95"]

print(pe <- 1 - hr)
print(lo.pe <- 1 - up.hr)
print(up.pe <- 1 - lo.hr)

# model 3: VRC01 30mg/kg vs. placebo
model3cox <- coxph(Surv(fudays_new, status_new80ls_cox) ~ cov, data=dat3)
print(summ3cox <- summary(model3cox)$conf.int)

hr = summ3cox[, "exp(coef)"]
lo.hr = summ3cox[, "lower .95"]
up.hr = summ3cox[, "upper .95"]

print(pe <- 1 - hr)
print(lo.pe <- 1 - up.hr)
print(up.pe <- 1 - lo.hr)

q(save="no")
