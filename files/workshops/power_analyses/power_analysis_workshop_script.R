rm(list=ls())
hablar::set_wd_to_script_path()

library(tidyverse)
library(pwr)

# ---- 
# Power for t-test, analytical calculation

# Define parameters
d <-c(0.2,0.5,0.8)
N <- seq(15, 200, 1)
alpha <- 0.05

# Calculate power for each effect size
X <- expand_grid(d, N)
power <- with(X, pwr.t.test(d = d, n = N, 
                            sig.level = alpha, 
                            type = "paired")$power)

# Plot
X %>% mutate(
  d = ordered(d)) %>%
  ggplot(aes(x=N,y=power, 
             color=d))+
  geom_hline(yintercept = 0.8, linewidth=0.4, lty=2)+
  geom_line(linewidth=1.5)+
  scale_x_continuous(breaks = seq(0,300,25))+
  labs(y="Power (1 - β)",
       x="Sample size (N. participants)")


# ---- 
# Simulation one-way anova

# simulate 1 dataset and run ANOVA test
set.seed(123)
group1 <- rnorm(30, mean = 50, sd = 10)
group2 <- rnorm(30, mean = 55, sd = 10)
group3 <- rnorm(30, mean = 52, sd = 10)

data <- data.frame(
  value = c(group1, group2, group3),
  group = factor(rep(1:3, each = 30))
)

head(data)
summary(aov(value ~ group, data))

# compute standardized effect size for F test
M_i <- c(50, 55, 52) # means
f <- sqrt(sum((M_i - mean(M_i))^2)/3) / 10

# what would be the power for this N?
pwr.anova.test(k = 3, n = 30, f = f, sig.level = 0.05)


# write custom function to simulate data
simulate_data <- function(N_per_group){
  
  group1 <- rnorm(N_per_group, mean = 50, sd = 10)
  group2 <- rnorm(N_per_group, mean = 52, sd = 10)
  group3 <- rnorm(N_per_group, mean = 55, sd = 10)
  
  data <- data.frame(
    value = c(group1, group2, group3),
    group = factor(rep(1:3, each = N_per_group))
  )
  
  return(data)
}

# function that run test and return p-value
run_test <- function(data){
  
  anova_summary <- summary(aov(value ~ group, data))
  p_value <- anova_summary [[1]]$`Pr(>F)`[1]
  return(p_value)
  
}

# run simulations for 1 sample size
N_sim <- 10^3
p_values <- rep(NA,N_sim)
for(i in 1:N_sim){
  p_values[i] <- run_test(simulate_data(30))
}

# visualize p-values distribution
hist(p_values)

# statistical power
mean(p_values < 0.05)


# Run simulations for multiple effect size 
# (that is, compute power curve)
N <- seq(30, 100, 10)
N_sim <- 10^3
simres <- expand_grid(iteration=1:N_sim, N=N)
str(simres)
simres$p <- NA
for(i in 1:nrow(simres)){
  simres$p[i] <- run_test(simulate_data(simres$N[i]))
}

# custom function to compute binomial standard errors (for plotting)
binomSE <- function (v) {
  sqrt((mean(v) * (1 - mean(v)))/length(v))
}

# plot results and compare to analytical calculation
alpha <- 0.05
simres %>%
  group_by(N) %>%
  summarize(power = mean(p<alpha),
            se = 2*binomSE(p<alpha),
            exact = pwr.anova.test(k=3, n=N, f=f, sig.level=alpha)$power) %>%
  ggplot(aes(x=N, y=power))+
  geom_errorbar(aes(ymin=power-se, 
                    ymax=power+se),
                width=0)+
  geom_point()+
  geom_line(aes(y=exact), col="red")


# Run tests computing different types of power (complete and specific comparisons)
run_test <- function(data){
  
  anova_summary <- summary(aov(value ~ group, data))
  p_value_aov <- anova_summary [[1]]$`Pr(>F)`[1]
  
  p_value_t1 <- t.test(data$value[data$group=="1"], 
                       data$value[data$group=="2"],
                       alternative="less")$p.value
  
  p_value_t2 <- t.test(data$value[data$group=="2"],
                       data$value[data$group=="3"], 
                       alternative="less")$p.value
  
  res <- c(p_value_aov, p_value_t1, p_value_t2)
  names(res) <- c("aov", "1<2", "2<3")
  
  return(res)
  
}

# run simulations (warning: this may tbe slow to complete)
N <- seq(30, 300, 10)
N_sim <- 10^3
simres <- expand_grid(iteration=1:N_sim, N=N)
str(simres)
simres$p_aov <- NA
simres$p_1vs2 <- NA
simres$p_2vs3 <- NA
for(i in 1:nrow(simres)){
  res <- run_test(simulate_data(simres$N[i]))
  simres$p_aov[i] <- res['aov']
  simres$p_1vs2[i] <- res['1<2']
  simres$p_2vs3[i] <- res['2<3']
}

# format results for plotting
long_res <- pivot_longer(simres, 
                         cols = starts_with("p_"), 
                         names_prefix = "p_", 
                         names_to = "test", 
                         values_to = "p_value")

alpha <- 0.05

simres %>%
  group_by(N) %>%
  summarize(power = mean(p_1vs2<alpha & p_2vs3<alpha),
            se = 2*binomSE(p_1vs2<alpha & p_2vs3<alpha)) %>%
  mutate(test = 'complete') -> complete_res

long_res %>%
  group_by(N, test) %>%
  summarize(power = mean(p_value<alpha),
            se = 2*binomSE(p_value<alpha)) ->long_avg

long_avg <- rbind(long_avg, complete_res)

long_avg %>%
  ggplot(aes(x=N, y =power, 
             ymin=power-se, 
             ymax=power+se,
             color=test))+
  geom_errorbar(width=0)+
  geom_point()+
  geom_line()

# do plot
simres %>%
  ggplot(aes(x=p_1vs2, y=p_2vs3)) +
  geom_point(alpha=0.5) +
  facet_wrap(.~N, ncol=7) +
  labs(x='p-value group 1 vs. group 2',
       y='p-value group 1 vs. group 2')



