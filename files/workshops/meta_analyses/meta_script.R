rm(list=ls())
hablar::set_wd_to_script_path()

# meta analysis workshop
library(tidyverse)
library(metafor)

# -------------------------------------------------------------- #
# example 1 meta-analysis
power_pose <- read.csv("./data/power_pose.csv")
str(power_pose)

# compute effect sizes

# function for computing pooled SD
pooled_SD <- function(sdA,sdB,nA,nB){
  sqrt((sdA^2*(nA-1) + sdB^2*(nB-1))/(nA+nB-1))
}

# function for Cohen's d
cohen_d_se <- function(n1,n2,d){
  sqrt((n1+n2)/(n1*n2) + (d^2)/(2*(n1+n2)))
}

# apply functions & create new columns in the dataset (using dplyr::mutate)
power_pose <- power_pose %>%
  mutate(pooled_sd =pooled_SD(sd_high_power, sd_low_power, n_high_power, n_low_power), 
         mean_diff =  mean_high_power - mean_low_power,
         cohen_d = mean_diff/pooled_sd,
         SE = cohen_d_se(n_high_power, n_low_power,cohen_d))

# plots
par(mfrow=c(1,2))
plot(with(power_pose, n_high_power + n_low_power),power_pose$SE, xlab="n", ylab="SE(Cohen's d)")
plot(with(power_pose, n_high_power + n_low_power),power_pose$cohen_d, xlab="n", ylab="Cohen's d")

# calculate also the t-test
power_pose$t <- with(power_pose, cohen_d / sqrt(1/n_high_power + 1/n_low_power))
power_pose$p <- 1 - pt(power_pose$t, df=with(power_pose, n_high_power + n_low_power - 2))


# run meta-analysis
res <- rma(cohen_d, SE^2, data = power_pose)
res


# methods for 'rma' objects
confint(res)
predict(res)

influence(res)
plot(influence(res))

funnel(res)
regtest(res)

# forest plot
forest(res)

# prettier forest plot
forest(res, slab=study, xlim=c(-5,2), xlab="Cohen's d",
       ilab=round(cbind(with(power_pose, n_high_power + n_low_power),
                        mean_low_power, sd_low_power, mean_high_power, sd_high_power),digits=2), 
       ilab.xpos=c(-3.5, -2.5,-2,-1.25,-0.75), cex=.75,
       header="Authors", mlab="")
par(cex=.75, font=2)
text(c(-2.5,-2,-1.25,-0.75), res$k+2, c("mean", "Std.", "mean", "Std."))
text(c(-2.25,-1),     res$k+3, c("Low power", "High power"))
text(-3.5, res$k+2, "N")

text(-5, -1, pos=4, cex=0.75, 
     bquote(paste("RE Model (Q = ",
                  .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
                  ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",
                  .(formatC(res$I2, digits=1, format="f")), "%)")))


# fixed-effect model
summary(rma(cohen_d, SE^2, data = power_pose, method="FE"))


# -------------------------------------------------------------- #
### Exercise RE meta-analysis

# 1
towels <- read.csv("./data/towels.csv")
str(towels)

# 2
load('./data/thirdwave.rda')
str(ThirdWave)


# -------------------------------------------------------------- #
# Meta-regression

# transform duration in a dummy variable
ThirdWave$duration_dummy <- ifelse(ThirdWave$InterventionDuration=="long", 1, 0)
res2 <- rma(yi=TE, vi=seTE^2, data = ThirdWave, mods = duration_dummy)
print(res2)

# LRT test
# refit models using maximum likelihood estimation (ML)
res0 <- rma(yi=TE, vi=seTE^2, data = ThirdWave, method="ML")
res2 <- rma(yi=TE, vi=seTE^2, data = ThirdWave, mods = duration_dummy, method="ML")

# run likelihood ratio test of nested models
anova(res0, res2)

# 
funnel(res0)
regtest(res0)

# -------------------------------------------------------------- #
# multilevel / three-level model

load('./data/Chernobyl.rda')
str(Chernobyl)

model_cherobyl <- rma.mv(yi = z, V = var.z, 
                     slab = author,
                     data = Chernobyl,
                     random = ~ 1 | author/es.id, # equivalent formulations below:
                     # random = list(~ 1 | author, ~ 1 | interaction(author,es.id)), 
                     # random = list(~ 1 | author, ~ 1 | es.id),
                     method = "REML")

summary(model_cherobyl)

funnel(model_cherobyl)
ranktest(model_cherobyl)

forest(model_cherobyl)

# -------------------------------------------------------------- #
### Exercise three-level model

library(metafor)
dat <- dat.konstantopoulos2011
head(dat)

# -------------------------------------------------------------- #
# funnel plot for ThirdWave dataset

res <- rma(yi=TE, vi=seTE^2, data = ThirdWave)
funnel(res, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), refline=0, legend=TRUE, xlab="Effect size (standardized mean difference)", refline2=res$b)

# "Egger’s regression test" for funnel plot asymmetry
regtest(res)

with(ThirdWave, summary(lm(I(TE/seTE) ~ I(1/seTE))))

ThirdWave %>% 
  mutate(y = TE/seTE, x = 1/seTE) %>% 
  lm(y ~ x, data = .) %>% 
  summary()

