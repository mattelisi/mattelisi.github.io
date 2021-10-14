	# ------------------------------------------------------------------------
	#
	# linear and generalized linear mixed-effects model tutorial
	#
	# Matteo Lisi, 2015
	#
	# ------------------------------------------------------------------------
	
	# clear workspace
	rm(list=ls())
	
	# load libraries
	library(lme4) 	 # contains functions for mixed-effects model analysis
	library(lattice) # additional plotting functions
	
	# ------------------------------------------------------------------------
	# example 1: sleepstudy
	
	# plot data
	str(sleepstudy)
	xyplot(Reaction ~ Days | Subject, sleepstudy, type = c("g","p","r"),
	       index = function(x,y) coef(lm(y ~ x))[1],
	       xlab = "Days of sleep deprivation",
	       ylab = "Average reaction time (ms)", aspect = "xy")
	
	## fit model
	
	# fully parametrized variance-covariance matrix 
	# in addition to the two variances (random slope and intercept)
	# also the covariance is estimated
	sleep1.m <- lmer(Reaction ~ 1 + Days + (1 + Days|Subject), data = sleepstudy)
	summary(sleep1.m)
	VarCorr(sleep1.m) # note that here you have the extra field for the correlation between the two random effects.
	dotplot(ranef(sleep1.m, condVar=T),scale="free") # visualize the conditional modes
	
	# random slope and intercept model
	# slope and intercept are considered as independent (i.e., with 0 covariance)
	sleep2.m <- lmer(Reaction ~ Days + (1|Subject) + (0 + Days|Subject), data = sleepstudy)
	summary(sleep2.m)
	VarCorr(sleep2.m) # this extract the standard deviation the random effects
	
	# likelihood ratio test between the two model
	anova(sleep1.m, sleep2.m)
	
	# likelihood ratio test to obtain a p value for the fixed-effect predictor
	fixef(sleep2.m)
	anova(sleep2.m, update(sleep2.m, .~. - Days))
	
	# but the better solution is to use bootstrap, and to do that
	# the easiest way is to use the function confint().
	# this function will use the method added by the lme4 library
	# for objects of class "merMod" (that is models fitted with lmer)
	# (remember that R is an object-oriented language!)
	# by default it will be a parametric bootstrap, with each
	# new simulation generating new values from the appropriate distributions (see also ?bootMer)
	confint(sleep2.m, nsim=1000) # equivalent to confint.merMod(sleep2.m, nsim=1000)
	
	# to see the help of the method for merMod fits type:
	?confint.merMod
	
	# diagnostic plot: visualize the model residuals
	par(mfrow=c(1,2))
	plot(fitted(sleep2.m),resid(sleep2.m), xlab="fitted values", ylab="residuals") # against fitted values if more than one predictor
	plot(jitter(sleepstudy$Days),resid(sleep2.m), xlab="Days", ylab="residuals") # but here one could plot it agains the preditor Days in this way
	boxplot(resid(sleep2.m)~sleepstudy$Days, xlab="Days", ylab="residuals")
	abline(h=0,lty=2)
	# hist(resid(sleep2.m)) # this makes a simple histogram of the residuals
	qqnorm(resid(sleep2.m))
	qqline(resid(sleep2.m),lty=2)
	
	
	# ------------------------------------------------------------------------
	# example 2: ergoStool
	data(ergoStool,package="MEMSS")
	str(ergoStool)
	# data are from an ergonomic experiment, where 9 subjects  evaluated 
	# the difficulty to arise for each of 4 types of stool
	# the effor is rated in the Borg scale of perceived exertion (ranging from 6 to 20)
	
	# plot
	barchart(effort ~ Type, ergoStool,group=Subject,xlab = "Type of stool", ylab = "Effort to arise",stack=FALSE)
	
	# plot separately variability within and between subjects
	par(mfrow=c(1,2))
	xcoord<-barplot(with(ergoStool,tapply(effort,Subject,mean)),ylim=c(0,20), ylab = "Effort to arise", xlab = "Subject")
	require(Hmisc)	# this package has a handy function for errobars
	errbar(xcoord,with(ergoStool,tapply(effort,Subject,mean)),with(ergoStool,tapply(effort,Subject,mean))+with(ergoStool,tapply(effort,Subject,sd)),with(ergoStool,tapply(effort,Subject,mean))-with(ergoStool,tapply(effort,Subject,sd)),pch="",add=T)
	xcoord<-barplot(with(ergoStool,tapply(effort,Type,mean)),ylim=c(0,20), ylab = "Effort to arise", xlab = "Type of stool")
	errbar(xcoord,with(ergoStool,tapply(effort,Type,mean)),with(ergoStool,tapply(effort,Type,mean))+with(ergoStool,tapply(effort,Type,sd)),with(ergoStool,tapply(effort,Type,mean))-with(ergoStool,tapply(effort,Type,sd)),pch="",add=T)
	
	# fit model
	stool1.m <- lmer(effort ~ Type + (1|Subject), ergoStool)
	summary(stool1.m)
	
	# use model parameters to test contrasts of interests
	confint(stool1.m, parm=4:6) # vs TypeT1
	
	# you can switch the contrast matrix of the factor to test different contrasts
	contrasts(ergoStool$Type) # visualize the contrast matrix; T2, T3, T4 are tested against T1
	
	# you can also adjust the confidence level of the interval to correct for 
	# multiple comparisons
	confint(stool1.m, parm=4:6, level = 1 - 0.05/6) # T2, T3, T4 vs T1 
	
	stool2.m <- lmer(effort ~ Type + (1|Subject), within(ergoStool, Type <- relevel(Type, ref = "T2")))
	confint(stool2.m, parm=5:6, level = 1 - 0.05/6) # T3, T4 vs T2 
	
	stool2.m <- lmer(effort ~ Type + (1|Subject), within(ergoStool, Type <- relevel(Type, ref = "T2")))
	confint(stool2.m, parm=6, level = 1 - 0.05/6) # T4, vs TypeT3 
	
	# diagnostic plot
	par(mfrow=c(1,2))
	boxplot(resid(stool1.m)~ergoStool$Type, xlab="Stool type", ylab="residuals")
	abline(h=0,lty=2)
	qqnorm(resid(stool1.m))
	qqline(resid(stool1.m),lty=2)
	
	# Notes:
	#
	# Not a great advantage of using mixed models here (small and balanced data set)
	# but when there are umbalanced data or large samples, the flexibility of mixed models
	# becomes important in the estimation of parameters. In ANOVA (or generally fixed effects linear models) 
	# the estimation is by least-squares, and requires predictors (i.e., the column of model matrix X) to be linearly 
	# independent. In practice this is very difficult to determine.
	#
	# In addition inferences from classic fixed-effects only models like the ANOVA 
	# apply strictly speaking only to the sample studied and not the general population.
	# In practice is not a big problem to think of them as applying to the population,
	# UNLESS, e.g., we want to predict the score that a general member of the population would give to 
	# a particular stool type. We can make the prediction also with the ANOVA, but then it is 
	# difficult to assess the variability of that prediction.
	
	
	
	# ------------------------------------------------------------------------
	#
	# GENERALIZED linear mixed-effects models
	#
	# ------------------------------------------------------------------------
	
	# ------------------------------------------------------------------------
	# example 1: the blink study 
	
	# load dataset
	setwd("~/Dropbox/works 2015/tutorial mixed models")	# chenage it to set the working directory to where you have the data file
	bridge <- read.table("blinkStudy.txt",header=T,sep="\t")
	str(bridge)
	
	# durations are taken from a uniform distribution between 250 and 500
	# therefore we take the ratio of durations divided by the expected value
	# or first moment of the distribution, that is 1/2*(250+500)=375
	bridge$ratio <- bridge$DUR/375
	
	# plot 
	
	# average values with same durations for plotting
	bridgeMean <- with(bridge, aggregate(RESP,list(ratio,SUBJ),mean))
	colnames(bridgeMean) <- c("ratio","SUBJ","RESP")
	bridgeMean$n <- with(bridge, aggregate(RESP,list(ratio,SUBJ),length))$x
	
	# subfunctions to fit data in individual panels with glm()
	panel.psyfun <- function(x, y, n, lnk = "logit", ...) {
			xy.glm <- glm(cbind(n * y, n * (1 - y)) ~ x,  binomial(lnk)) 
			rr <- current.panel.limits()$xlim 
			xx <- seq(rr[1], rr[2], len = 100) 
			yy <- predict(xy.glm, data.frame(x = xx),  type = "response") 
			panel.lines(xx, yy, ...) 
			PSE <- -coef(xy.glm)[1]/coef(xy.glm)[2]
			panel.lines(c(PSE,PSE),c(0.5,-1),lty=1,...)
			} 
	
	xyplot(RESP ~ ratio | SUBJ, data = bridgeMean, 
		type = c("p"), col="black", cex=bridgeMean$n/5, subscripts = TRUE, id=bridgeMean$n,
		xlab = "ratio", 
		panel = function(x,y,id,subscripts,...){
			panel.abline(h=0.5,lty=2,col="dark grey")
			panel.abline(v=1,lty=2,col="dark grey")
			panel.xyplot(x,y,...)
			panel.psyfun(x,y,id[subscripts],lnk="probit", lwd=2, col="black")},			 
		ylab = expression(paste("P('longer' response | ratio)")), 
		auto.key = list(space = "right",title = "ratio ", cex = 0.75), 
		par.settings = list(strip.background = list(col = "lightgrey")))
	
	
	## fit model with random location and scale parameters
	bridge.m <- glmer(RESP ~ ratio + (ratio|SUBJ), data=bridge, family = binomial(link=probit))
	summary(bridge.m)
	
	## diagnostic plot
	require(ggplot2)
	ggplot()+geom_point(aes(fitted(bridge.m),y=residuals(bridge.m, type = "pearson")),pch=21)+stat_smooth(aes(fitted(bridge.m),y=residuals(bridge.m, type = "pearson")),method = "loess",se = T, n=5,col="black")+theme_bw()+labs(x="fitted", y = "residuals")
	
	ggplot()+geom_point(aes(x=bridge$ratio,y=residuals(bridge.m, type = "pearson"),color=bridge$SUBJ),size=2)+stat_smooth(aes(x=bridge$ratio,y=residuals(bridge.m, type = "pearson")),method = "loess",se = T, n=5,col="black")+theme_bw()+labs(x="ratio", y = "residuals")+guides(color=FALSE)
	
	## examine parameters
	fixef(bridge.m) # linear predictor parameters
	dotplot(ranef(bridge.m, condVar=T),scale="free") # visualize the conditional modes
	
	PSE <- unname(-fixef(bridge.m)[1]/fixef(bridge.m)[2]) # this give the PSE
	
	# the PSE is 1.15. Is this significantly greater than 1
	# or in other words, is there a significant underestimation of 
	# the duration during eye blinks?
	# We can test it with bootstrap
	
	# simple function to extract PSE
	myFUN <- function(.) {
		pse <- unname(-fixef(.)[1]/fixef(.)[2])
	}
	
	# perform bootstrap (warning: this is going to take some time...)
	bootResult <- bootMer(bridge.m, myFUN, 1000,.progress="txt")
	
	# compute CI from the quantiles of bootstrapped PSE distribution
	bootCI <- c(quantile(bootResult$t, probs = 0.025),quantile(bootResult$t, probs = 0.975))
	
	# make a nice plot of the bootstrapped distribution
	d<-density(bootResult$t)
	maxYplot <- max(d$y) 
	hd <- hist(bootResult$t, 
		main =" ", 
		xlab=expression(paste(-beta[italic("0")]/beta[italic("1")])), 
		ylab="", prob=T, freq =F, breaks = 30, 
		col= "grey", border ="white",
		yaxt="n",
		ylim=c(0,maxYplot+1),
		xlim=c(1,1.35))
	lines(d$x,d$y)
	require(Hmisc)
	lines(c(bootCI[1],bootCI[1]),c(-1,maxYplot+0.2), lty=2)
	lines(c(bootCI[2],bootCI[2]),c(-1,maxYplot+0.2), lty=2)
	boxed.labels(mean(bootResult$t),maxYplot+1,labels="95% confidence interval",cex=0.6,bg="white",col="black",border="white",xpad=2,ypad=1.5)
	text(quantile(bootResult$t, probs = 0.96),maxYplot*2/3,labels=expression(paste(10^3," rep.")),cex=0.6,pos=4)
	arrows(bootCI[1], maxYplot+0.2, bootCI[2], maxYplot+0.2, length = 0.1, code = 3)
	
	
	## compare glmm and individual fits 
	
	# compute location and scale parameters from individual fit with the function lmList
	bridge.single <- lmList(RESP ~ ratio |SUBJ, data=bridge, family = binomial(link=probit))
	PSE.single <- -coef(bridge.single)[,1]/coef(bridge.single)[,2]
	scale.single <- 1/coef(bridge.single)[,2]
	
	# compute individual parameters by adding the conditional modes to the fixed effect parameters of the glmm
	repmat <- function(X,m,n){
		# R equivalent of repmat (matlab)
		mx = dim(X)[1]
		nx = dim(X)[2]
		matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
	}
	
	bridge.mixed <- as.matrix(ranef(bridge.m)$SUBJ) + repmat(t(as.matrix(fixef(bridge.m))),11,1)
	PSE.mixed <- -bridge.mixed[,1]/bridge.mixed[,2]
	scale.mixed <- 1/bridge.mixed[,2]
	
	# plot parameters from mixed model and individual fit
	plot(PSE.single, scale.single,xlab=expression(paste(mu," (location parameter)")),ylab=expression(paste(sigma," (scale parameter)")),pch=19,cex.lab=1.2, xlim=c(0.85,1.5), ylim=c(0.1,0.5))
	grid()
	text(PSE.single, scale.single,labels=rownames(coef(bridge.single)),pos=2)
	points(PSE.mixed, scale.mixed)
	arrows(PSE.single, scale.single,PSE.mixed, scale.mixed, length=0.1)
	points(mean(PSE.mixed), mean(scale.mixed),pch=2)
	legend("bottomright",c("within-subjects","mixed-model","population"),pch=c(19,1,2),bty="n")
