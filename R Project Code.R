set.seed(23)
e=rnorm(100,0,2) #100 residuals for 1st series
plot.ts(e)
phi=c(0.3,0,0.4) #set phi values
y=c()
y[0]=0
y[1]=6+e[1]
y[2]=6+phi[1]*(y[1]-6)+e[2]
y[3]=6+phi[1]*(y[2]-6)+phi[2]*(y[1]-6)+e[3]
y[4]=6+phi[1]*(y[3]-6)+phi[2]*(y[2]-6)+phi[3]*(y[1]-6)+e[4]

for(i in 4:100)
  {y[i]=6+phi[1]*(y[i-1]-6)+phi[2]*(y[i-2]-6)+phi[3]*(y[i-3]-6)+e[i]} #generate series 1
  
plot.ts(y)

y1=y[51:100] #throw out first 50 obvs y1 is the 50 obvs series
plot.ts(y1)

set.seed(23)
e2=rnorm(100,0,5) #residuals for series 2
plot.ts(e2)

y2=c()
y2[51]=6+phi[1]*(y1[50]-6)+phi[2]*(y1[49]-6)+phi[3]*(y1[48]-6)+e2[51]
y2[52]=6+phi[1]*(y2[51]-6)+phi[2]*(y1[50]-6)+phi[3]*(y1[49]-6)+e2[52]
y2[53]=6+phi[1]*(y2[52]-6)+phi[2]*(y2[51]-6)+phi[3]*(y1[50]-6)+e2[53]
y2[54]=6+phi[1]*(y2[53]-6)+phi[2]*(y2[52]-6)+phi[3]*(y2[51]-6)+e2[54]

for(i in 54:100)
{y2[i]=6+phi[1]*(y2[i-1]-6)+phi[2]*(y2[i-2]-6)+phi[3]*(y2[i-3]-6)+e2[i]}
y20=y2[51:100] #throw out first 100 y20 is the second 50 obvs series
plot.ts(y20)

set.seed(23)
e3=rnorm(100,0,8) #residuals for series 3
plot.ts(e3)

y3=c()
y3[51]=6+phi[1]*(y20[50]-6)+phi[2]*(y20[49]-6)+phi[3]*(y20[48]-6)+e3[51]
y3[52]=6+phi[1]*(y3[51]-6)+phi[2]*(y20[50]-6)+phi[3]*(y20[49]-6)+e3[52]
y3[53]=6+phi[1]*(y3[52]-6)+phi[2]*(y3[51]-6)+phi[3]*(y20[50]-6)+e3[53]
y3[54]=6+phi[1]*(y3[53]-6)+phi[2]*(y3[52]-6)+phi[3]*(y3[51]-6)+e3[54]

for(i in 54:100)
{y3[i]=6+phi[1]*(y3[i-1]-6)+phi[2]*(y3[i-2]-6)+phi[3]*(y3[i-3]-6)+e3[i]}
y30=y3[51:100] #throw out first 50 and y30 is the 50 obvs series for series 3
plot.ts(y30)

yt=c()
yt=y1 #Add in the first 50 obvs
plot.ts(yt)
yt[51:100]=y20 #Then the second 50
plot.ts(yt)
yt[101:150]=y30 #Now the third 50
plot.ts(yt,main = "150 obvs Time Series") #AR(3) process and start trying to fit it
acf(yt, main="Sample ACF for Yt")
pacf(yt,main="Sample PACF for Yt")

diff1=diff(yt)
plot.ts(diff1,main="Differenced Series Y't")
acf(diff1,main="Sample ACF Y't")
pacf(diff1, main="Sample PACF Y't")

diff2=diff(diff1)
plot.ts(diff2,main="Second Difference Y''t")
acf(diff2,main="Sample ACF Y''t")
pacf(diff2,main="Sample PACF Y''t")

resid=c()
resid[1:50]=e[51:100]
resid[51:100]=e2[51:100]
resid[101:150]=e3[51:100]
plot.ts(resid)
resid_sq=resid*resid
plot.ts(resid_sq)
acf(resid_sq)
pacf(resid_sq)

qqnorm(yt); qqline(yt)
qqnorm(resid); qqline(resid)

arima(yt, order = c(3,0,0))
arima300=arima(yt,order=c(3,0,0))
t(confint(arima300))

arima(yt, order = c(3,0,1))
arima301=arima(yt,order=c(3,0,1))
t(confint(arima301))

arima(yt, order = c(3,0,2))
arima302=arima(yt,order=c(3,0,2))
t(confint(arima302))

arima(yt, order = c(3,0,3))
arima303=arima(yt,order=c(3,0,3))
t(confint(arima303))

arima(yt, order = c(4,0,2))
arima402=arima(yt,order=c(4,0,2))
t(confint(arima402))

arima302=arima(yt,order=c(3,0,2))
arimaresiduals=arima302$residual
t(confint(arima302))

plot.ts(arimaresiduals,main="Residuals of ARMA(3,2)")
shapiro.test(arimaresiduals)
ks.test(arimaresiduals,pnorm,alternative = c("two.sided", "less", "greater"),exact = NULL)
qqnorm(arimaresiduals, main="QQ Plot ARMA(3,2) Residuals"); qqline(arimaresiduals)
hist(arimaresiduals, freq=FALSE,main="Histogram ARMA(3,2) Residuals")
t <- -100:1000/10
lines(t, dnorm(t, 0, 5), col = 'darkgreen')
acf(arimaresiduals, main="ACF Plot of ARMA(3,2) Residuals")
pacf(arimaresiduals, main="PACF Plot of ARMA(3,2) Residuals")


residuals_squared=arimaresiduals*arimaresiduals
plot.ts(residuals_squared, main="Plot of Residuals Squared")

Box.test(arimaresiduals, lag = 1, type="Ljung")
Box.test(arimaresiduals, lag = 30, type="Ljung")
Box.test(arimaresiduals, lag = 60, type="Ljung")
Box.test(arimaresiduals, lag = 90, type="Ljung")
Box.test(arimaresiduals, lag = 120, type="Ljung")
Box.test(arimaresiduals, lag = 150, type="Ljung")


library(TSA)
McLeod.Li.test(arima302,main="McLeod-Li Test on ARIMA(3,0,2)")

arch01=garch(arimaresiduals,order=c(0,1),trace=F)
loglik01=logLik(arch01)
print(loglik01)
summary(arch01) 
arch01resid=arch01$residuals
aic01=-2*loglik01+2*4*(150/145)
print(aic01)
t(confint(arch01))

arch02=garch(arimaresiduals,order=c(0,2),trace=F)
loglik02=logLik(arch02)
print(loglik02)
summary(arch02) 
arch02resid=arch02$residuals
aic02=-2*loglik02+2*4*(150/145)
print(aic02)
t(confint(arch02))

arch04=garch(arimaresiduals,order=c(0,4),trace=F)
loglik04=logLik(arch04)
print(loglik04)
summary(arch04) 
arch04resid=arch04$residuals
aic04=-2*loglik04+2*4*(150/145)
print(aic04)
t(confint(arch04))


arch03=garch(arimaresiduals,order=c(0,3),trace=F)
loglik03=logLik(arch03)
print(loglik03)
summary(arch03) #ARCH(3) Model is most suitable with highest log-likelihood
arch03resid=arch03$residuals
aic03=-2*loglik03+2*4*(150/145)
print(aic03)
t(confint(arch03))

arch11=garch(arimaresiduals,order=c(1,1),trace=F)
loglik11=logLik(arch11)
print(loglik11)
summary(arch11)
arch11resid=arch11$residuals
aic11=-2*loglik11+2*4*(150/145)
print(aic11)
t(confint(arch11))

arch13=garch(arimaresiduals,order=c(1,3),trace=F)
loglik13=logLik(arch13)
print(loglik13)
summary(arch13)
arch13resid=arch13$residuals
aic13=-2*loglik13+2*4*(150/145)
print(aic13)
t(confint(arch13))

arch23=garch(arimaresiduals,order=c(2,3),trace=F)
loglik23=logLik(arch23)
print(loglik23)
summary(arch03) 
arch03resid=arch03$residuals
aic23=-2*loglik23+2*4*(150/145)
print(aic23)
t(confint(arch23))


Box.test(arch03resid, lag = 1,type="Ljung")
Box.test(arch03resid, lag = 30,type="Ljung")
Box.test(arch03resid, lag = 60,type="Ljung")
Box.test(arch03resid, lag = 90,type="Ljung")
Box.test(arch03resid, lag = 120,type="Ljung")


plot.ts(arch03resid,main="Residuals of ARCH(3)")
acf(arch03resid[4:150],main="Sample ACF for Residuals of ARCH(3)")
pacf(arch03resid[4:150],main="Sample PACF for Residuals of ARCH(3)")
qqnorm(arch03resid[4:150], main="QQ Plot ARCH(3) Residuals");qqline(arch03resid)
hist(arch03resid,freq = FALSE, main="Histogram ARCH(3) Residuals") 
t <- -100:1000/10
lines(t, dnorm(t, 0, 1), col = 'darkblue')
shapiro.test(arch03resid)
ks.test(arch03resid,pnorm,alternative = c("two.sided", "less", "greater"),exact = NULL)


archresid_sq=arch03resid*arch03resid
plot.ts(archresid_sq)
acf(archresid_sq[4:150])
pacf(archresid_sq[4:150])

