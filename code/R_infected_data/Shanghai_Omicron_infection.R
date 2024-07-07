########################################################
#注: 运行本程序前先运行backcalculate.xf脚本.
#2022上海Omicron
#由报道病例估计新发病例
########################################################
#rm(list=ls())
library(readxl)
rc=read_excel('../reported_case.xlsx')#根据表格reported_case.xlsx（文件夹中有该表格）所处位置进行导入
attach(rc)

###感染-报道间隔分布，已知感染-报道间隔分布的间隔时长c(0.5,1,2,3,4,5,6)，对应频数c(9,34,29,3,2,7,9)，拟合其分布
library(MASS)
library(survival)
library(fitdistrplus)
#每个元素对应真实间隔时长/天：c(0.5,1,2,3,4,5,6)，这里这样做是为了避免频率直方图把0.5和1归为一组
infect_report_interval=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5)
infect_report_interval_frequency=c(9,34,29,3,2,7,9)#频数
data_infect_report_interval=rep(infect_report_interval,infect_report_interval_frequency)
#write.csv(data_infect_report_interval,file='../data_infect_report_interval.csv')

fitg <- fitdist(data_infect_report_interval, "gamma")
summary(fitg)
plot(fitg)
fitW <- fitdist(data_infect_report_interval, "weibull")
fitln <- fitdist(data_infect_report_interval, "lnorm")
gofstat(list(fitW, fitg, fitln), fitnames=c("Weibull", "gamma", "lognormal"))
# Weibull    gamma lognormal
#Akaike's Information Criterion 340.9162 336.3352  335.9116
#对数正态分布，根据AIC最小准则，选择对数正态分布来拟合感染-报道间隔分布
par(cex.lab=1.1)
denscomp(list(fitW, fitg, fitln), xlab="time interval",ylab='Density',fitlwd=2,main = "(b)",legendtext=c("Weibull AIC=340.9162", "gamma AIC=336.3352", "lognormal AIC=335.9116"))
bc=seq(0.5,6.5,0.1)
gam=dgamma(bc,fitg$estimate[1],fitg$estimate[2])
Weib=dweibull(bc,fitW$estimate[1],fitW$estimate[2])
lognor=dlnorm(bc,fitln$estimate[1],fitln$estimate[2])
interval_fit=data.frame(time=bc,gam=gam,Weib=Weib,lognor=lognor)
write.csv(interval_fit,file='../interval_fit.csv')

summary(fitln)
#Fitting of the distribution ' lnorm ' by maximum likelihood 
#Parameters : 
#  estimate Std. Error
#meanlog 0.7490948 0.07066096
#sdlog   0.6814296 0.04996436
#Loglikelihood:  -165.9558   AIC:  335.9116   BIC:  340.9768 
#Correlation matrix:
#  meanlog         sdlog
#meanlog  1.000000e+00 -2.508592e-11
#sdlog   -2.508592e-11  1.000000e+00
plot(fitln)

######################################################################################3
a=fitln$estimate[1];b=fitln$estimate[2]
p = rep(0,1000);
p[1] = plnorm(1.5, meanlog = a,sdlog = b);
for(i in 2:1000){p[i] = plnorm(i+0.5, meanlog = a,sdlog = b) - plnorm(i-0.5, meanlog = a,sdlog = b)}
plot(p[1:50]);




######################################## 逆卷积估计新发病例infec ########################################
#黑色：新发数据，红色：报道数据，绿色：最终估计的报道数据
result = backcalculate.xf(lambda0 = rep(1,length(reported_case)+200),D = reported_case, p = p, n = 5);
result$time
# 检查结果(实心黑点为确诊数据data, 红线为感染者数估计值(逆卷积的结果), 黑线为感染者数实际值infec)
#par(mar = c(0.5, 1, 5, 1) ) 
par(cex.lab=1.1)
plot(reported_case,pch=1,col='red',xlim = c(result$time[1],length(reported_case)),main="(C)",xlab='time',ylab='cases')
lines(result$time[1]:result$time[2],result$estimation,col="black",lwd=2);
legend("topright", pch=c(1,NA),lty=c(0,1),lwd=c(1,2),col =   c('red',"black"),legend = c("newly reported cases","newly infected cases"))
time=seq(as.Date("2022/2/24"),as.Date("2022/7/2"),by="day")
infec=data.frame(time=time,infection=round(result$estimation))
write.csv(infec,file='../infected_data.csv')

######################################## 逆卷积估计新发病例infec ########################################
