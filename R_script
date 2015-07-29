library(quantmod)
library(fBasics)
library(forecast)
library(fGarch)
library(MTS)
library(rugarch)
library(rmgarch)
library(ggplot2)
library(stringi)

###################################################################
### Univariate Time Series Analysis
### BP

# load BP
#(0), data preparation and cleaning
getSymbols("BP",from="2015-01-01", to="2015-01-31", src = "google")
head(BP)
dim(BP)
bp=BP[,4] 
head(bp)
chartSeries(BP,theme="white")
logbp=log(bp+1)  # take log bp
head(logbp)
chartSeries(logbp,theme="white")
logbp.dif=diff(logbp)  # take difference
head(logbp.dif)
par(mfrow=c(3,1))
plot.ts(logbp.dif1)
plot.ts(logchevron.dif1)
plot.ts(logex.dif1)

logbp.dif1 = logbp.dif[2:20,]  # get rid of NA

#(a), basic statistics
basicStats(logbp.dif1) 

#(b), test if mean is zero
t.test(logbp.dif1)
# analyze acf and pacf
acf(logbp.dif1,lag=30)
pacf(logbp.dif1, lag=30)]

#(c), perform skewness test
s3=skewness(logbp.dif1)
T=length(logbp.dif1)
tst=s3/sqrt(6/T)
tst
pv=2*pnorm(tst)
pv

#(d), perform kurtosis test.
k4=kurtosis(logbp.dif1)
tst=k4/sqrt(24/T)
tst
pv=2*(1-pnorm(tst))
pv
#(e), empirical density plot
par(mfrow=c(1,3))
den.bp=density(logbp.dif1)
plot(den.bp$x,den.bp$y, xlab="log return",ylab="density",main = "BP", type="l",col="red")
xfit<-seq(min(logbp.dif1),max(logbp.dif1),length=40) 
yfit<-dnorm(xfit,mean=mean(logbp.dif1),sd=sd(logbp.dif1)) 
lines(xfit, yfit, col="black", lwd=2)

#(f), ARIMA model
ar(logbp.dif1)
fit <- arima(logbp.dif1, order=c(4,0,0), include.mean=F)
fit
tsdiag(fit)

#(g), 1-st month forecast
forecast(fit,20)
plot(forecast(fit,20))
plot(fit)
plot(logbp.dif1)

#(h), test if there is ARCH effect
Box.test(fit$residuals^2, lag = 10, type = "Ljung")
par(mfrow=c(1,2))
acf(logbp.dif1^2, lag = 24)
pacf(logbp.dif1^2, lag = 24)

#(i), GARCH model
m1=garchFit(~arma(4,0)+garch(1,1),data=logbp.dif1,trace=F)
summary(m1)
plot(m1)
predict(m1,20)
#(j)
m11=garchFit(~arma(4,0)+garch(1,1),data=logbp.dif1,trace=F,cond.dist="std")
summary(m11)
plot(m11)

bp.Garch.fit=garchFit(~arma(4,0)+garch(1,1),data=logbp.dif1,trace=F,cond.dist="sstd")
summary(bp.Garch.fit)
plot(bp.Garch.fit)
predict(bp.Garch.fit)


###################################################################
### CHEVRON

# load Chevron
#(0)
getSymbols("CVX",from="2015-01-01", to="2015-01-31", src = "google")
head(CVX)
chartSeries(CVX,theme="white")
dim(CVX)
chevron=CVX[,4] 
head(chevron)
# take log
logchevron=log(chevron+1)
acf(logchevron)
pacf(logchevron)
head(logchevron)
chartSeries(logchevron,theme="white")
# take difference
logchevron.dif=diff(logchevron)
head(logchevron.dif)
chartSeries(logchevron.dif,theme="white")
# get rid of NA
logchevron.dif1 = logchevron.dif[2:20,]
#(a)
basicStats(logchevron.dif1) 
#(b)
t.test(logchevron.dif1)
# analyze acf and pacf
acf(logchevron.dif1,lag=100)
pacf(logchevron.dif1, lag=100)
#(c)
# Perform skewness test
s3=skewness(logchevron.dif1)
T=length(logchevron.dif1)
tst=s3/sqrt(6/T)
tst
pv=2*pnorm(tst)
pv
#(d)
# Perform kurtosis test.
k4=kurtosis(logchevron.dif1)
tst=k4/sqrt(24/T)
tst
pv=2*(1-pnorm(tst))
pv
#(e)
den.chevron=density(logchevron.dif1)
plot(den.chevron$x,den.chevron$y, xlab="log return",ylab="density",main = "Chevron", type="l",col="red")
xfit<-seq(min(logchevron.dif1),max(logchevron.dif1),length=40) 
yfit<-dnorm(xfit,mean=mean(logchevron.dif1),sd=sd(logchevron.dif1)) 
lines(xfit, yfit, col="black", lwd=2)
#(f)
fit2 <- arima(logchevron.dif1,order=c(5,0,2),include.mean=F)
fit2
#(g)
forecast(fit2, 20)
plot(forecast(fit2, 20))
chev=predict(fit2,4)
chev
tsdiag(fit2)
#(h)
Box.test(fit2$residuals^2, lag = 10, type = "Ljung")
par(mfrow=c(1,2))
acf(logchevron.dif1^2, lag = 24)
pacf(logchevron.dif1^2, lag = 24)
#(i)
m2=garchFit(~arma(5,2)+garch(1,1),data=logchevron.dif1,trace=F)
summary(m2)
plot(m2)
predict(m2,20)
#(j)
chevron.Garch.fit=garchFit(~arma(5,0)+garch(1,1),data=logchevron.dif1,trace=F,cond.dist="sstd")
summary(chevron.Garch.fit)
plot(chevron.Garch.fit)



###################################################################
### EXXON

# load Exxon
#(0)
getSymbols("XOM",from="2015-01-01", to="2015-01-31", src = "google")
head(XOM)
dim(XOM)
ex=XOM[,4] 
head(ex)
chartSeries(ex,theme="white")
# take log
logex=log(ex+1)
head(logex)
chartSeries(logex,theme="white")
# take difference
logex.dif=diff(logex)
head(logex.dif)
chartSeries(logex.dif,theme="white", main="Exxon")
# get rid of NA
logex.dif1 = logex.dif[2:20,]
#(a)
basicStats(logex.dif1) 
#(b)
t.test(logex.dif1)
#(c)
# Perform skewness test
s3=skewness(logex.dif1)
T=length(logex.dif1)
tst=s3/sqrt(6/T)
tst
pv=2*pnorm(tst)
pv
#(d)
# Perform kurtosis test.
k4=kurtosis(logex.dif1)
tst=k4/sqrt(24/T)
tst
pv=2*(1-pnorm(tst))
pv
#(e)
den.exxon=density(logex.dif1)
plot(den.exxon$x,den.exxon$y, xlab="log return",ylab="density",main = "Exxon", type="l", col="red")
xfit<-seq(min(logex.dif1),max(logex.dif1),length=40) 
yfit<-dnorm(xfit,mean=mean(logex.dif1),sd=sd(logex.dif1)) 
lines(xfit, yfit, col="black", lwd=2)
#(f)
ar(logex.dif1)
fit3 <- arima(logex.dif1, order=c(2,0,0), include.mean=F)
fit3
tsdiag(fit3)
#(g)
forecast(fit3,20)
plot(forecast(fit3,20))
#(h)
Box.test(fit3$residuals^2, lag = 10, type = "Ljung")
par(mfrow=c(1,2))
acf(logex.dif1^2, lag = 24)
pacf(logex.dif1^2, lag = 24)
#(i)
m3=garchFit(~arma(2,0)+garch(1,1),data=logex.dif1,trace=F)
summary(m3)
plot(m3)
predict(m3,20)
#(j)
exxon.Garch.fit=garchFit(~arma(2,0)+garch(2,1),data=logex.dif1,trace=F,cond.dist="sstd")
summary(exxon.Garch.fit)
plot(exxon.Garch.fit)

garchstd.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(4,0)),
                           distribution.model="sstd")
bp.Garch.fit = ugarchfit(spec=garchstd.spec, data=logbp.dif1,
                           solver.control=list(trace = 0))

garchstd.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(5,0)),
                           distribution.model="sstd")
chevron.Garch.fit = ugarchfit(spec=garchstd.spec, data=logchevron.dif1,
                         solver.control=list(trace = 0))

garchstd.spec = ugarchspec(variance.model = list(garchOrder=c(2,1)),
                           mean.model = list(armaOrder=c(2,0)),
                           distribution.model="sstd")
exxon.Garch.fit = ugarchfit(spec=garchstd.spec, data=logex.dif1,
                              solver.control=list(trace = 0))
# (k) Obtain 1-step ahead mean and volatility forecasts using the fitted ARMA-GARCH model with
# Student-t innovations with 95% confidence intervals for the first month in 2015.
bp.Garch.forecast = ugarchforecast(bp.Garch.fit, n.ahead=20,n.start=503)
chevron.Garch.forecast = ugarchforecast(chevron.Garch.fit, n.ahead=20,n.start=503)
exxon.Garch.forecast = ugarchforecast(exxon.Garch.fit, n.ahead=20,n.start=503)
bp.Garch.forecast
chevron.Garch.forecast
exxon.Garch.forecast
# Calculate 95% confidence intervals
bp.lcl<- as.data.frame(bp.Garch.forecast@forecast$seriesFor-
                            1.96*bp.Garch.forecast@forecast$sigmaFor)
names(bp.lcl) <- "BP.lcl"
bp.ucl<- as.data.frame(bp.Garch.forecast@forecast$seriesFor+1.96*bp.Garch.forecast@forecast$sigmaFor)
names(bp.ucl) <- "BP.ucl"
bp.pred <- as.data.frame(bp.Garch.forecast@forecast$seriesFor)
names(bp.pred) <- "BP.predicted"
bp.pred.df <-cbind(bp.lcl, bp.pred, bp.ucl)
bp.pred.df
# Calculate 95% confidence intervals
chevron.lcl<- as.data.frame(chevron.Garch.forecast@forecast$seriesFor-1.96*chevron.Garch.forecast@forecast$sigmaFor)
names(chevron.lcl) <- "Chevron.lcl"
chevron.ucl<- as.data.frame(chevron.Garch.forecast@forecast$seriesFor+1.96*chevron.Garch.forecast@forecast$sigmaFor)
names(chevron.ucl) <- "Chevron.ucl"
chevron.pred <- as.data.frame(chevron.Garch.forecast@forecast$seriesFor)
names(chevron.pred) <- "Chevron.predicted"
chevron.pred.df <-cbind(chevron.lcl, chevron.pred, chevron.ucl)
chevron.pred.df
# Calculate 95% confidence intervals
ex.lcl<- as.data.frame(exxon.Garch.forecast@forecast$seriesFor-1.96*exxon.Garch.forecast@forecast$sigmaFor)
names(ex.lcl) <- "Exxon.lcl"
ex.ucl<- as.data.frame(exxon.Garch.forecast@forecast$seriesFor+1.96*exxon.Garch.forecast@forecast$sigmaFor)
names(ex.ucl) <- "Exxon.ucl"
ex.pred <- as.data.frame(exxon.Garch.forecast@forecast$seriesFor)
names(ex.pred) <- "Exxon.predicted"
ex.pred.df <-cbind(ex.lcl, ex.pred, ex.ucl)
ex.pred.df

par(mfcol=c(1,2))
plot(bp.Garch.forecast, which=1)
mtext("BP",side=3,cex=0.8)
plot(bp.Garch.forecast, which=3)
mtext("BP",side=3,cex=0.8)
plot(chevron.Garch.forecast, which=1)
mtext("Chevron",side=3,cex=0.8)
plot(chevron.Garch.forecast, which=3)
mtext("Chevron",side=3,cex=0.8)
plot(exxon.Garch.forecast, which=1)
mtext("Exxon",side=3,cex=0.8)
plot(exxon.Garch.forecast, which=3)
mtext("Exxon",side=3,cex=0.8)

####################################################################
### Multivariate Time Series Analysis

#(l) cross correlation for all between all series
all = cbind(logbp.dif1,logchevron.dif1,logex.dif1)
ccm(all,lag=20)
#m1=Eccm(all)
#m1
MTSplot(all)

# scatterplot of returns
par(mfrow=c(3,1))
#Correlation plot between BP and Chevron
plot( coredata(logbp.dif1), coredata(logchevron.dif1), xlab="BP", ylab="Chevron", type="p", pch=16, lwd=2, col="blue")
abline(h=0,v=0)

#Correlation plot between BP and Exxon
plot( coredata(logbp.dif1), coredata(logex.dif1), xlab="BP", ylab="Exxon",
      type="p", pch=16, lwd=2, col="blue")
abline(h=0,v=0)

plot( coredata(logchevron.dif1), coredata(logex.dif1), xlab="Chevron", ylab="Exxon",
      type="p", pch=16, lwd=2, col="blue")
abline(h=0,v=0)

cor(logbp.dif1,logchevron.dif1)
cor(logbp.dif1,logex.dif1)
cor(logchevron.dif1,logex.dif1)

#(m) 
cor.fun = function(x){
  cor(x)[1,2]
}
cov.fun = function(x){
  cov(x)[1,2]
}

roll.cov = rollapply(as.zoo(all), FUN=cov.fun, width=30,
                     by.column=FALSE, align="right")
roll.cor = rollapply(as.zoo(all), FUN=cor.fun, width=30,
                     by.column=FALSE, align="right")
par(mfrow=c(2,1))
plot(roll.cov, main="30-day rolling covariances",
     ylab="covariance", lwd=2, col="blue")
grid()
abline(h=cov(all)[1,2], lwd=2, col="red")
plot(roll.cor, main="30-day rolling correlations",
     ylab="correlation", lwd=2, col="blue")
grid()
abline(h=cor(all)[1,2], lwd=2, col="red")

par(mfrow=c(1,1))
mean(roll.cor) 
min(roll.cor)
max(roll.cor) 
mean(roll.cov) 
min(roll.cov) 
max(roll.cov)

#(n), univariate normal GARCH(1,1) for each series

garch11.spec=ugarchspec(mean.model = list(armaOrder = c(0,0)),
                        variance.model=list(garchOrder = c(1,1), model="sGARCH"), distribution.model="norm")
dcc.garch11.spec = dccspec(uspec=multispec(replicate(3,garch11.spec)), dccOrder = c(1,1), distribution="mvnorm")
dcc.garch11.spec
dcc.fit = dccfit(dcc.garch11.spec, all)
class(dcc.fit)
slotNames(dcc.fit)
names(dcc.fit@mfit)
names(dcc.fit@model)
show(dcc.fit)
plot(dcc.fit)
coef(dcc.fit)
likelihood(dcc.fit)

#(o)
# Plot Conditional Covariance
par(mar=c(3,3,4,1))
plot(dcc.fit,which=3, series=c(1,2))
mtext("DCC with Normal Distribution",side=3,cex=0.9,,padj=-1.2)
par(mar=c(3,3,4,1))
plot(dcc.fit,which=3, series=c(1,3))
mtext("DCC with Normal Distribution",side=3,cex=0.9,padj=-1.2)
par(mar=c(3,3,4,1))
plot(dcc.fit,which=3, series=c(2,3))
mtext("DCC with Normal Distribution",side=3,cex=0.9,padj=-1.2)

#Plot Conditional Correlation
par(mar=c(3,3,4,1))
plot(dcc.fit,which=4, series=c(1,2))
mtext("DCC with Normal Distribution",side=3,cex=0.9,,padj=-1.2)
par(mar=c(3,3,4,1))
plot(dcc.fit,which=4, series=c(1,3))
mtext("DCC with Normal Distribution",side=3,cex=0.9,,padj=-1.2)
par(mar=c(3,3,4,1))
plot(dcc.fit,which=4, series=c(2,3))
mtext("DCC with Normal Distribution",side=3,cex=0.9,,padj=-1.2)

# Plot EW Portfolio with conditional density VaR limits
par(mar=c(3,3,4,1))
plot(dcc.fit,which=5, series=c(1,2))
mtext("DCC with Normal Distribution",side=3,cex=0.9,,padj=-1.2)
par(mar=c(3,3,4,1))
plot(dcc.fit,which=5, series=c(1,3))
mtext("DCC with Normal Distribution",side=3,cex=0.9,,padj=-1.2)
par(mar=c(3,3,4,1))
plot(dcc.fit,which=5, series=c(2,3))
mtext("DCC with Normal Distribution",side=3,cex=0.9,,padj=-1.2)
mtext("DCC with Normal Distribution",side=3,cex=0.9,padj=-1.2)

#(p)
# forecasting conditional volatility and correlations
dcc.fcst = dccforecast(dcc.fit, n.ahead=20)
dcc.fcst

# Plot the DCC Forecasted Conditional Covariances.
plot(dcc.fcst, which=3, series=c(1,2))
plot(dcc.fcst, which=3, series=c(1,3))
plot(dcc.fcst, which=3, series=c(2,3))

# Plot EW Portfolio with conditional density VaR limits
plot(dcc.fcst, which=5, series=c(1,2))

#plot the DCC forecasted conditional correlation
corr.dcc.fcst = as.data.frame(rcor(dcc.fcst))
names(corr.dcc.fcst)
row.names(corr.dcc.fcst)
corr.dcc.fcst.df <- as.data.frame(t(corr.dcc.fcst))
names(corr.dcc.fcst.df)
corr.dcc.fcst.df$X2014.12.31.CVX.Close <- as.vector(rep(c(0,1,0),times=20))
corr.dcc.fcst.df$X2014.12.31.XOM.Close <- as.vector(rep(c(0,0,1),times=20))

plot.BP.CVX.df <- subset(corr.dcc.fcst.df,corr.dcc.fcst.df$X2014.12.31.CVX.Close==1, select=BP.Close)
plot.BP.XOM.df <- subset(corr.dcc.fcst.df,corr.dcc.fcst.df$X2014.12.31.XOM.Close==1, select=BP.Close )
plot.CVX.XOM.df <- subset(corr.dcc.fcst.df,corr.dcc.fcst.df$X2014.12.31.XOM.Close==1, select=CVX.Close )

# Extract the fitted conditional correlations into a dataframe.
corr.dcc.fit = as.data.frame(rcor(dcc.fit))
corr.dcc.fit.df <- as.data.frame(t(corr.dcc.fit)) # transposte the rows and columns
corr.dcc.fit.df[1:10,] # verify the transposed columns
corr.dcc.fit.df$X2014.12.31.CVX.Close <- as.vector(rep(c(0,1,0),times=503)) 
corr.dcc.fit.df$X2014.12.31.XOM.Close <- as.vector(rep(c(0,0,1),times=503))

# subset the dataframe to extract the desired conditional correlation pairs from the fitted dataframe.
plot.BP.XOM.fit.df <- subset(corr.dcc.fit.df,corr.dcc.fit.df$X2014.12.31.XOM.Close==1, select=BP.Close )
plot.BP.CVX.fit.df <- subset(corr.dcc.fit.df,corr.dcc.fit.df$X2014.12.31.CVX.Close==1, select=BP.Close )
plot.CVX.XOM.fit.df <- subset(corr.dcc.fit.df,corr.dcc.fit.df$X2014.12.31.XOM.Close==1, select=CVX.Close )

# Combine the fitted and forecasted conditional correlations into a single dataframe.
plot.BP.CXV.fit.df <- rbind(plot.BP.XOM.fit.df,plot.BP.XOM.fit.df)
plot.BP.XOM.fit.df <- rbind(plot.BP.CVX.fit.df,plot.BP.CVX.fit.df)
plot.CVX.XOM.fit.df <- rbind(plot.CVX.XOM.fit.df,plot.CVX.XOM.fit.df)

# Plot BP-CVX Correlation
row.names(plot.BP.CXV.fit.df)[735:1006]
x.start <- 735 # begin plotting in the 3rd quarter of 2014 (October 1,2014)
x.end <- nrow(plot.BP.CXV.fit.df)
fit.length <- 1006 - x.start -19
forcast.length <- nrow(plot.BP.CVX.df)
line.color <- c(rep("darkgrey",fit.length), 
                rep("red", forcast.length)) 
df <- data.frame(
  x = x.start:x.end, 
  y = plot.BP.CXV.fit.df$BP.Close[x.start:x.end],
  col=line.color)

p1 <- ggplot(df, aes(x=x, y=y)) +
  geom_line(aes(colour=col, group=1)) +
  scale_colour_identity() +
  labs(y="Correlation", x = "Days") +
  ggtitle("BP-CVX Forecasted Conditional Correlation\n January 2013 through December 2014")
p1 + theme_bw() + theme(plot.title = element_text(color="black",size=20,face="bold"))


# Plot BP-XOM Correlation
df <- data.frame(x = x.start:x.end,y = plot.BP.XOM.fit.df$BP.Close[x.start:x.end],col = line.color)

p2 <- ggplot(df, aes(x=x, y=y)) +
  geom_line(aes(colour=col, group=1)) +
  scale_colour_identity() +
  labs(y="Correlation", x = "Days") +
  ggtitle("BP-XOM Forecasted Conditional Correlation\n January 2013 through December 2014")
p2 + theme_bw() + theme(plot.title = element_text(color="black",size=20,face="bold"))

# Plot CVX-XOM Correlation
df <- data.frame(
  x = x.start:x.end,
  y = plot.CVX.XOM.fit.df$CVX.Close[x.start:x.end],
  col =line.color)
p3 <- ggplot(df, aes(x=x, y=y)) +
  geom_line(aes(colour=col, group=1)) +
  scale_colour_identity() +
  labs(y="Correlation", x = "Days") +
  ggtitle("CVX-XOM Forecasted Conditional Correlation\n January 2013 through December 2014")
p3 + theme_bw() + theme(plot.title = element_text(color="black",size=20,face="bold"))

#q
garchstd.spec=ugarchspec(mean.model = list(armaOrder = c(0,0)),
                         variance.model=list(garchOrder = c(1,1), model="sGARCH"), distribution.model="std")
# dcc specification - GARCH(1,1) for conditional correlations
dcc.garchstd.spec = dccspec(uspec=multispec(replicate(3,garchstd.spec)), dccOrder = c(1,1), distribution="mvt")
dcc.garchstd.spec
dccstd.fit = dccfit(dcc.garchstd.spec, all)
show(dccstd.fit)
coef(dccstd.fit)
likelihood(dccstd.fit)
# Plot In Sample Conditional Covariance for the Student-t Distribution
par(mar=c(3,3,4,1))
plot(dccstd.fit,which=3, series=c(1,2))
mtext("DCC with Student-t Distribution",side=3,cex=0.9,,padj=-1.2)
par(mar=c(3,3,4,1))
plot(dccstd.fit,which=3, series=c(1,3))
mtext("DCC with Student-t Distribution",side=3,cex=0.9,padj=-1.2)
par(mar=c(3,3,4,1))
plot(dccstd.fit,which=3, series=c(2,3))
mtext("DCC with Student-t Distribution",side=3,cex=0.9,padj=-1.2) # add a title in the top margin of the plot
#Plot In Sample Conditional Correlation
par(mar=c(3,3,4,1))
plot(dccstd.fit,which=4, series=c(1,2))
mtext("DCC with Student-t Distribution",side=3,cex=0.9,,padj=-1.2)
par(mar=c(3,3,4,1))
plot(dccstd.fit,which=4, series=c(1,3))
mtext("DCC with Student-t Distribution",side=3,cex=0.9,,padj=-1.2)
par(mar=c(3,3,4,1))
plot(dccstd.fit,which=4, series=c(2,3))
mtext("DCC with Student-t Distribution",side=3,cex=0.9,,padj=-1.2)

###################################################################
#Fitting Multivariate model to our time series.
all = cbind(logbp.dif1,logchevron.dif1,logex.dif1)
all2 = as.matrix(all)    #Transforming data frame to matrix
#MTSplot(all)   #Plots time series
ccm(all2,5)    #calculates cross correlation among the time series equivalent to univariate ACF 
mq(all2,10)    #Computes multivariate Ljung-Box test statistics

#Fitting model and removing insignificant parameters
VARorder(all2,maxp=8)         # suggests an AR order
m1=VAR(all2,p=2)             #fitting AR model with suggested order
MTSdiag(m1)                  #Model diagnostics

#Removing insignificant parameters
m2=refVAR(m1, thres=1.0)
MTSdiag(m2)
VARpred(m2,4)       #Computes the forecasts of the VAR model, the associated errors of forecasts and the mean squared errors of forecasts. 


summary(m2)
