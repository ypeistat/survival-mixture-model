rm(list=ls())
library(survival)
library(glmnet)
library(foreach)
library(ggplot2)
library(caret)
library(survival)
library(survminer)
library(survivalROC)
#set.seed(123456)

source("fcts.R")
source("cox_EM.R")
source("cox_est.R")


########################
dat = read.csv('data_y.csv',header = T,stringsAsFactors=F)
age = dat$age
sex <- ifelse(dat$sex=="female",1,0)
X <- cbind(age,sex)
time <- dat$event_censoring_days
status <- ifelse(dat$status=="TRUE",1,0)
data <- as.data.frame(cbind(time, status,X))
###
result <- coxph(Surv(time, status) ~ age + sex, data = data)
summary(result)

###
mmax <- 5
n <- dim(data)[1]
p <- 2
#lambda0 <- seq(sqrt(log(n)/n)*0.1, sqrt(log(n)/n)*0.6,length.out = 10)
#############################
lambda <- sqrt(log(n)/n)/14
res <- cox_EM(data, lambda, mmax)
# require(PFMC)
# fit1 = pmixcox(Time = data$time,
#                Delta = data$status,
#                X = data[,-1:-2],
#                K = 1,
#                tparm = 0,
#                abstol = 1e-4,
#                reltol = 1e-6,
#                seed = 1)
# 
# fit2 = pmixcox(Time = data$time,
#                         Delta = data$status,
#                         X = data[,-1:-2],
#                         K = 2,
#                         tparm = 0,
#                         abstol = 1e-4,
#                         reltol = 1e-6,
#                         seed = 1)
# 
# 
# 
# fit3 = pmixcox(Time = data$time,
#                Delta = data$status,
#                X = data[,-1:-2],
#                K = 3,
#                tparm = 0,
#                abstol = 1e-4,
#                reltol = 1e-6,
#                seed = 1)
# 
# 
# fit4 = pmixcox(Time = data$time,
#                Delta = data$status,
#                X = data[,-1:-2],
#                K = 4,
#                tparm = 0,
#                abstol = 1e-4,
#                reltol = 1e-6,
#                seed = 1)
# 
# 
# fit5 = pmixcox(Time = data$time,
#                Delta = data$status,
#                X = data[,-1:-2],
#                K = 5,
#                tparm = 0,
#                abstol = 1e-4,
#                reltol = 1e-6,
#                seed = 1)
# 
# 
# fit6 = pmixcox(Time = data$time,
#                Delta = data$status,
#                X = data[,-1:-2],
#                K = 6,
#                tparm = 0,
#                abstol = 1e-4,
#                reltol = 1e-6,
#                seed = 1)
# 
# fit7 = pmixcox(Time = data$time,
#                Delta = data$status,
#                X = data[,-1:-2],
#                K = 7,
#                tparm = 0,
#                abstol = 1e-4,
#                reltol = 1e-6,
#                seed = 1)
# 
# k = 1:6
# n = 1159
# p = 2
# dof = k-1+k*p
# mixloglik = c(-560,-525,-516,-479,-466,-459)
# bic = -2*mixloglik + dof*log(n)
# plot(k,mixloglik)





