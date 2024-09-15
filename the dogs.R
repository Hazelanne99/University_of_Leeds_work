
################################################################################
setwd("C:/Users/hazel/OneDrive/Documents/UoL/Statistical Theory and methods/Coursework")
load("mysample")

#createsample(201202015)

#save(mysample, file = "mysample.RData")

table(mysample$Breed)

install.packages("viridis")
library(viridis)

#Data cleaning

mysample

colnames(mysample)

rownames(mysample)

#before cleaning
#d = 22
#s = 214
#h = 23

mysample1 = mysample[mysample$Rehomed != 99999,]
#from 265 obs to 256
#d = 21
#s = 207
#h = 23

table(mysample1$Breed)

mysample1$Rehomed

mysample2 = mysample1[!is.na(mysample1$Breed),]
# from 256 observations to 250, removed 6 NA dogs

cleaned_mysample = mysample2[mysample2$Visited >= 0,]
#from 250 obs to 248
#d = 21
#s = 205
#h = 22

table(cleaned_mysample$Breed)

################################################################################

#Data exploration 

#Separated each breed into their own datasets

dobermann = cleaned_mysample[cleaned_mysample$Breed == "Dobermann",]

staffy = cleaned_mysample[cleaned_mysample$Breed == 
                            "Staffordshire Bull Terrier",]

highland = cleaned_mysample[cleaned_mysample$Breed == 
                              "West Highland White Terrier",]



#summaries for each breeds information
#taking: mean of visiting time, rehoming time and health,
#and the most common value for age, reason for in the shelther and if returned

summary(dobermann)

summary(staffy)

summary((highland))

library(dplyr)

a = table(dobermann$Age)
b = table(dobermann$Returned)
c = table(dobermann$Reason)

summary_table_dobermann = dobermann %>%
  summarise(
    Mean_Visited = mean(Visited),
    Mean_Rehomed = mean(Rehomed),
    Mean_Health = mean(Health),
    Total_Breed = n(),
    Common_Age = names(which.max(a)),
    Common_Reason = names(which.max(c)),
    Common_Returned = names(which.max(b))
  
  )

d = table(staffy$Age)
e = table(staffy$Returned)
f = table(staffy$Reason)

summary_table_staffy = staffy %>%
  summarise(
    Mean_Visited = mean(Visited),
    Mean_Rehomed = mean(Rehomed),
    Mean_Health = mean(Health),
    Total_Breed = n(),
    Common_Age = names(which.max(d)),
    Common_Reason = names(which.max(f)),
    Common_Returned = names(which.max(e))
    
  )


g = table(highland$Age)
h = table(highland$Returned)
i = table(highland$Reason)

summary_table_highland = highland %>%
  summarise(
    Mean_Visited = mean(Visited),
    Mean_Rehomed = mean(Rehomed),
    Mean_Health = mean(Health),
    Total_Breed = n(),
    Common_Age = names(which.max(g)),
    Common_Reason = names(which.max(i)),
    Common_Returned = names(which.max(h))
    
  )

summary_table_highland
summary_table_staffy
summary_table_dobermann


#notable differences: total breeds, age, rehoming time
#possibly notable differences: health
#visiting time looks roughly the same for all breeds and the most common answer 
#for if returned is no for all three


summary_rehomed = cleaned_mysample$Rehomed

dobermann_rehomed = c(dobermann$Rehomed)
staffy_rehomed = c(staffy$Rehomed)
highland_rehomed = c(highland$Rehomed)

#boxplots for all the dogs comparing rehoming times
par(mfrow = c(1,3))

par(mar = c(5,5,4,2))

#?boxplot

#png(file = "myplot.png",  width = 1000, height = 1000, res = 200, units = "px")  

boxplot(dobermann_rehomed, col = magma(15)[6],
        main = "Rehomed time of Dobermann",
        xlab = "Dobermann",
        ylab = "Rehomed time in weeks")
axis(side =2, at = seq(0, 50, by = 5))
boxplot(staffy_rehomed, col = magma(15)[9],
        main = "Rehomed time of Staffordshire Bull Terrier",
        xlab = "Staffordshire Bull Terrier",
        ylab = "Rehomed time in weeks")
axis(side =2, at = seq(0, 50, by = 5))
boxplot(highland_rehomed, col = magma(15)[12],
        main = "Rehomed time of West Highland White Terrier",
        xlab = "West Highland White Terrier",
        ylab = "Rehomed time in weeks")
axis(side =2, at = seq(0, 50, by = 5))


par(mfrow = c(1,1))

all_dogs = list(dobermann_rehomed,staffy_rehomed,highland_rehomed)
all_fogs = rownames(c("Dobermann","Shaffordshire Bull Terrier","West Highland "))

boxplot(all_dogs, col = c(magma(9)[3],magma(9)[6],magma(9)[9]),
        main = "Rehomed time of all three dogs",
        xlab = "West Highland White Terrier",
        ylab = "Rehomed time in weeks")
axis(side =2, at = seq(0, 50, by = 5))



#dev.off()

#quantiles
quantile(dobermann_rehomed)
# 0%  25%  50%  75% 100% 
#4   11   16   22   42 
max(dobermann_rehomed)
#42
min(dobermann_rehomed)
#4
#median = 16

#IQR = 11


quantile(staffy_rehomed)
# 0%  25%  50%  75% 100% 
# 1   11   18   26   51 
#max = 51
#min = 1
#IQR = 15
#median = 18



quantile(highland_rehomed)
# 0%  25%  50%  75% 100% 
#  7.00 13.00 18.50 25.75 48.00 
#max = 48
#min = 7
#IQR = 12.75
#median = 18.5

quantile(transformed_h)
#0%      25%      50%      75%     100% 
#2.079442 2.639057 2.962128 3.286402 3.891820 
#max = 3.891820
#min = 2.079442
#IQR = 0.647
#median = 2.962128


#sample mean of each dog's rehoming time 

sample_mean_dobermann = mean(dobermann_rehomed)
sample_mean_dobermann
#17.85714 = 1st sample moment 


sample_variance_dobermann = sd(dobermann_rehomed)^2
sample_variance_dobermann
#5.291837 = 2nd sample moment ???????????? 111.1286


sd_dobermann = sd(dobermann_rehomed)
sd_dobermann
#2.357208 - sample sd ???????????????????? 10.54175


for (i in dobermann_rehomed) {
  xs = (i - sample_mean_dobermann)^2
  xsum = sum(xs)
  skewness_dobermann = ((1/21)*xsum)/sd_dobermann^3
}
skewness_dobermann
#0.0000741978 - 3rd sample moment ??????????????? 0.0000008295566


sample_mean_staffy = mean(staffy_rehomed)
sample_mean_staffy
#19.46341 = 1st sample moment 


sample_variance_staffy = sd(staffy_rehomed)^2
sample_variance_staffy
#0.5407694 - 2nd sample moment 110.8577

sd_staffy = sd(staffy_rehomed)
sd_staffy
#0.7360147 - sd 10.5289


for (i in staffy_rehomed) {
  xs = (i - sample_mean_staffy)^2
  xsum = sum(xs)
  skewness_staffy = ((1/21)*xsum)/sd_staffy^3
}
skewness_staffy
#3.828563 - 3rd sample moment 0.001250592



sample_mean_highland = mean(highland_rehomed)
sample_mean_highland
#20.22727 = 1st sample moment 


sample_variance_highland = sd(highland_rehomed)^2
sample_variance_highland
#4.649055 - 2nd sample moment 102.2792

sd_highland = sd(highland_rehomed)
sd_highland
#2.206907 - sd 10.11332


for (i in highland_rehomed) {
  xs = (i - sample_mean_highland)^2
  xsum = sum(xs)
  skewness_highland = ((1/21)*xsum)/sd_highland^3
}
skewness_highland
#1.246338 - 3rd sample moment  0.01295

################################################################################

#QQplots

par(mfrow = c(1,3))

par(mar = c(4,4,4,1))

mu_d = sample_mean_dobermann
sigma_d = sd_dobermann

qqnorm(dobermann_rehomed, col = magma(9)[3], lwd = 1,
       main = " Q-Q Plot of the Distribution for Dobermann"
       )
abline(a = mu_d, b = sigma_d, col = magma(9)[1], lwd =2)


#?colour

mu_s = sample_mean_staffy
sigma_s = sd(staffy_rehomed)

qqnorm(staffy_rehomed, col = magma(9)[6], lwd =1,
       main = " Q-Q Plot of the Distribution for Staffordshire Bull Terrier")
abline(a = mu_s, b = sigma_s, col = magma(9)[1], lwd = 2)


mu_h = sample_mean_highland
sigma_h = sd_highland

qqnorm(highland_rehomed, col = magma(9)[8], lwd = 1,
       main = " Q-Q Plot of the Distribution for West Highland White Terrier")
abline(a = mu_h, b = sigma_h, col = magma(9)[1], lwd =2)


transformed_h = log(highland_rehomed+1)
trans_mu = mean(transformed_h)
trans_sigma = sd(transformed_h)

qqnorm(transformed_h, lwd = 1, col = "green4",
       main = "Q-Q Plot for West Highland White Terrier after transformation")
abline(a = trans_mu, b = trans_sigma, col = "red")


#histogram of all the dogs rehome times


#mean and median lines added
par(mfrow = c(1,3))

par(mar = c(4,4,4,0))

#png(file = "histograms.png",  width = 1000, height = 1000, res = 200, 
 #   units = "px")  


hist(dobermann_rehomed, breaks = 10, right = FALSE, freq = FALSE,
     col = magma(15)[6],
     main = "Rehomed time for Dobermann",
     xlab = "Rehomed time in weeks",
     ylab = "Frequency Density")
abline(v=sample_mean_dobermann, lwd = 3, col =magma(15)[1])
abline(v=16, lty=2, lwd=3, col = magma(15)[4])
legend(x = 25, y = 0.04, 
       legend = c("Mean", "Median"), 
       lty = c(1, 2), 
       col = c(magma(15)[1], magma(15)[4]), 
       horiz = FALSE)


hist(staffy_rehomed, breaks = 15, right = FALSE, freq = FALSE,
     col = magma(15)[9],
     main = "Rehomed time for Staffordshire Bull Terrier",
     xlab = "Rehomed time in weeks",
     ylab = "Frequency Density")
abline(v=sample_mean_staffy, lwd = 3, col = magma(15)[1])
abline(v=18, lty=2, lwd=3, col=magma(15)[4])
legend(x = 25, y = 0.035, 
       legend = c("Mean", "Median"), 
       lty = c(1, 2), 
       col = c(magma(15)[1], magma(15)[4]), 
       horiz = FALSE)


hist(highland_rehomed, breaks = 5, right = FALSE, freq = FALSE,
     col = magma(15)[12],
     main = "Rehomed time for West Highland White Terrier",
     xlab = "Rehomed time in weeks",
     ylab = "Frequency Density")
abline(v=sample_mean_highland, lwd = 3, col = magma(15)[1])
abline(v=18.5, lty=2, lwd=3, col=magma(15)[4])
legend(x = 32, y = 0.06, 
       legend = c("Mean", "Median"), 
       lty = c(1, 2), 
       col = c(magma(15)[1], magma(15)[4]), 
       horiz = FALSE)

par(mar = c(4,4,4,1))

hist(transformed_h, breaks = 10, right = FALSE, freq = FALSE,
     col = "green4",
     main = "Histogram of Rehomed times for West Highland White Terrier
     after Transformation",
     xlab = "Rehomed time in weeks, log(x +1)",
     ylab = "Frequency Density")
abline(v=trans_mu, lwd = 2, col = "red")
abline(v=2.962128, lty=2, lwd=2, col="blue")
legend(x = 2.05, y = 1.4, 
       legend = c("Mean", "Median"), 
       lty = c(1, 2), 
       col = c("red", "blue"), 
       horiz = FALSE)


#dev.off()

#samples CDFs

set.seed(575)

par(mfrow = c(1,3), mar = c(4,4,4,1))


#png(file = "CDF.png",  width = 1000, height = 1000, res = 200,units = "px") 

Fn_d = ecdf(dobermann_rehomed)

Fn_d(100)

sum(dobermann_rehomed <= 100)/length(dobermann_rehomed) 

mu_d = sample_mean_dobermann
sigma_d = sd_dobermann

G <- function(x){
  
  return(pnorm(x, mean = mu_d, sd = sigma_d))
}

G(100)

plot(Fn_d, verticals = TRUE, pch = NA, col = magma(15)[6],
     main = "Sample CDF of Dobermann",
     xlab = "x (Rehomed time in weeks)",
     lwd = 2)
#abline(v=sample_mean_dobermann, lwd = 1, col = magma(15)[1])
#abline(v=16, lty=2, lwd=1, col=magma(15)[4])
#legend(x = 22, y = 0.40, 
       legend = c("Mean", "Median", "CDF of G(x)"), 
       lty = c(1, 2,1), 
       col = c(magma(15)[1], magma(15)[4], magma(15)[14]), 
       horiz = FALSE)


x <- 1:500   
lines(x, G(x), col = magma(15)[1], lwd = 2)


Fn_s = ecdf(staffy_rehomed)

Fn_s(100)

sum(staffy_rehomed <= 100)/length(staffy_rehomed) 

mu_s = sample_mean_staffy
sigma_s = sd_staffy

G <- function(x){
  
  return(pnorm(x, mean = mu_s, sd = sigma_s))
}

G(100)

plot(Fn_s, verticals = TRUE, pch = NA, col = magma(15)[9],
     lwd = 2,
     main = "Sample CDF of Staffordshire Bull Terrier",
     xlab = "x (Rehomed time in weeks)")
abline(v=sample_mean_staffy, lwd = 1, col = "red")
abline(v=18, lty=2, lwd=1, col="blue")
legend(x = 22, y = 0.40, 
       legend = c("Mean", "Median", "CDF of G(x)"), 
       lty = c(1, 2,1), 
       col = c("red", "blue", "purple"), 
       horiz = FALSE)


x <- 1:500   
lines(x, G(x), col = magma(15)[1], lwd =2)



Fn_h = ecdf(highland_rehomed)

Fn_h(100)

sum(highland_rehomed <= 100)/length(highland_rehomed) 

mu_h = sample_mean_highland
sigma_h = sd_highland

G <- function(x){
  
  return(pnorm(x, mean = mu_h, sd = sigma_h))
}

G(100)

plot(Fn_h, verticals = TRUE, pch = NA,, col = magma(15)[12],
     lwd = 2,
     main = "Sample CDF of West Highland White Terrier",
     xlab = "x (Rehomed time in weeks)")
abline(v=sample_mean_highland, lwd = 1, col = "red")
abline(v=18.5, lty=2, lwd=1, col="blue")
legend(x = 22, y = 0.40, 
       legend = c("Mean", "Median", "CDF of G(x)"), 
       lty = c(1, 2,1), 
       col = c("red", "blue", "purple"), 
       horiz = FALSE)

x <- 1:500   
lines(x, G(x), col = magma(15)[1], lwd = 2)



Fn_th= ecdf(transformed_h)

Fn_th(100)

sum(transformed_h <= 100)/length(transformed_h) 

mu_th = trans_mu
sigma_th = trans_sigma

G <- function(x){
  
  return(pnorm(x, mean = mu_th, sd =sigma_th))
}

G(100)

plot(Fn_th, verticals = TRUE, pch = NA,, col = "green4",
     lwd = 2,
     main = "Sample CDF of West Highland White Terrier after Transformation",
     xlab = "x (Rehomed time in weeks, log(x+1))")
legend(x = 22, y = 0.40, 
       legend = c("Mean", "Median", "CDF of G(x)"), 
       lty = c(1, 2,1), 
       col = c("red", "blue", "purple"), 
       horiz = FALSE)

x <- 1:500   
lines(x, G(x), col = "purple", lwd = 2)

#dev.off()


################################################################################

#Modelling and estimation

#MLE

set.seed(462)

x = rnorm(100)
x = x/sd(x) * sd_dobermann
x = x-mean(x) + sample_mean_dobermann
c('mean'=mean(x),'sd'=sd(x)) 


likelihood_d = function(x, par) {
  mu = par[1]  
  sigma = par[2]  
  n=length(dobermann_rehomed)
  
  log_likelihood_d = -(n/2)*(log(2*pi*sigma^2)) + (-1/(2*sigma^2)) * 
    sum((x-mu)^2)
  
  return(-log_likelihood_d)  
}

res_d = optim(par=c(.5,.5), likelihood_d, x=x)
res_d

cbind('direct'=c('mean'=sample_mean_dobermann,'sd'=sd_dobermann),
      'optim'=res_d$par)



set.seed(64)

x = rnorm(100)
x = x/sd(x) * sd_staffy
x = x-mean(x) + sample_mean_staffy
c('mean'=mean(x),'sd'=sd(x)) 


likelihood_s = function(x, par) {
  mu = par[1]  
  sigma = par[2]  
  n=length(staffy_rehomed)
  
  log_likelihood_s = -(n/2)*(log(2*pi*sigma^2)) + (-1/(2*sigma^2)) * 
    sum((x-mu)^2)
  
  return(-log_likelihood_s)  
}

res_s = optim(par=c(.5,.5), likelihood_s, x=x)
res_s

cbind('direct'=c('mean'=sample_mean_staffy,'sd'=sd_staffy),
      'optim'=res_s$par)



set.seed(4384)

x = rnorm(100)
x = x/sd(x) * sd_highland
x = x-mean(x) + sample_mean_highland
c('mean'=mean(x),'sd'=sd(x)) 


likelihood_h = function(x, par) {
  mu = par[1]  
  sigma = par[2]  
  n=length(highland_rehomed)
  
  log_likelihood_h = -(n/2)*(log(2*pi*sigma^2)) + (-1/(2*sigma^2)) * 
    sum((x-mu)^2)
  
  return(-log_likelihood_h)  
}

res_h = optim(par=c(.5,.5), likelihood_h, x=x)
res_h

cbind('direct'=c('mean'=sample_mean_highland,'sd'=sd_staffy),
      'optim'=res_h$par)



#estimators mean vs quartiles 

par(mar = c(3,3,3,3))

mu_d = sample_mean_dobermann
sigma_d = sd_dobermann

muhat1_d <- rep(NA, 1000)
muhat2_d<- rep(NA, 1000)
for(i in 1:1000){
  
  
  x <- rnorm(n = 21, mean = mu_d, sd = sigma_d)
  muhat1_d[i] <- mean(x)
  muhat2_d[i] <- quantile(x, type = 1)[3]
}

par(mfrow = c(1, 2))

hist(muhat1_d, xlim = range(c(muhat1_d, muhat2_d)))
abline(v = mu_d, col = "red3", lwd = 3)
abline(v = mean(muhat1_d), col = "blue", lty = 2, lwd = 3)
hist(muhat2_d, xlim = range(c(muhat1_d, muhat2_d)))
abline(v = mu_d, col = "red3", lwd = 3)
abline(v = mean(muhat2_d), col = "blue", lty = 2, lwd = 3)


sigma2hat1_d <- rep(NA, 1000)
sigma2hat2_d <- rep(NA, 1000)

for(i in 1:1000){
  x <- rnorm(n = 21, mean = mu_d, sd = sigma_d)
  
  sigma2hat1_d[i] <- sd(x)^2
  sigma2hat2_d[i] <- (9/10)*sd(x)^2
}

par(mfrow = c(1, 2))
hist(sigma2hat1_d, xlim = range(c(sigma2hat1_d, sigma2hat2_d)))
abline(v = sigma_d^2, col = "red3", lwd = 3)
abline(v = mean(sigma2hat1_d), col = "blue", lty = 2, lwd = 3)
hist(sigma2hat2_d, xlim = range(c(sigma2hat1_d, sigma2hat2_d)))
abline(v = sigma_d^2, col = "red3", lwd = 3)
abline(v = mean(sigma2hat2_d), col = "blue", lty = 2, lwd = 3)



mu_s = sample_mean_staffy
sigma_s = sd(staffy_rehomed)


muhat1_s <- rep(NA, 1000)
muhat2_s <- rep(NA, 1000)
for(i in 1:1000){
  
  
  x <- rnorm(n = 207, mean = mu_s, sd = sigma_s)
  muhat1_s[i] <- mean(x)
  muhat2_s[i] <- quantile(x, type = 1)[3]
}

par(mfrow = c(1, 2))

hist(muhat1_s, xlim = range(c(muhat1_s, muhat2_s)))
abline(v = mu_s, col = "red3", lwd = 3)
abline(v = mean(muhat1_s), col = "blue", lty = 2, lwd = 3)
hist(muhat2_s, xlim = range(c(muhat1_s, muhat2_s)))
abline(v = mu_s, col = "red3", lwd = 3)
abline(v = mean(muhat2_s), col = "blue", lty = 2, lwd = 3)


sigma2hat1_s <- rep(NA, 1000)
sigma2hat2_s <- rep(NA, 1000)

for(i in 1:1000){
  x <- rnorm(n = 207, mean = mu_s, sd = sigma_s)
  
  sigma2hat1_s[i] <- sd(x)^2
  sigma2hat2_s[i] <- (9/10)*sd(x)^2
}

par(mfrow = c(1, 2))
hist(sigma2hat1_s, xlim = range(c(sigma2hat1_s, sigma2hat2_s)))
abline(v = sigma_s^2, col = "red3", lwd = 3)
abline(v = mean(sigma2hat1_s), col = "blue", lty = 2, lwd = 3)
hist(sigma2hat2_s, xlim = range(c(sigma2hat1_s, sigma2hat2_s)))
abline(v = sigma_s^2, col = "red3", lwd = 3)
abline(v = mean(sigma2hat2_s), col = "blue", lty = 2, lwd = 3)


mu_h = sample_mean_highland
sigma_h = sd(highland_rehomed)

muhat1_h <- rep(NA, 1000)
muhat2_h <- rep(NA, 1000)
for(i in 1:1000){
  
  
  x <- rnorm(n = 22, mean = mu_h, sd = sigma_h)
  muhat1_h[i] <- mean(x)
  muhat2_h[i] <- quantile(x, type = 1)[3]
}

par(mfrow = c(1, 2))

hist(muhat1_h, xlim = range(c(muhat1_h, muhat2_h)))
abline(v = mu_h, col = "red3", lwd = 3)
abline(v = mean(muhat1_h), col = "blue", lty = 2, lwd = 3)
hist(muhat2_h, xlim = range(c(muhat1_h, muhat2_h)))
abline(v = mu_h, col = "red3", lwd = 3)
abline(v = mean(muhat2_h), col = "blue", lty = 2, lwd = 3)



sigma2hat1_h <- rep(NA, 1000)
sigma2hat2_h <- rep(NA, 1000)

for(i in 1:1000){
  x <- rnorm(n = 22, mean = mu_h, sd = sigma_h)
  
  sigma2hat1_h[i] <- sd(x)^2
  sigma2hat2_h[i] <- (9/10)*sd(x)^2
}

par(mfrow = c(1, 2))
hist(sigma2hat1_h, xlim = range(c(sigma2hat1_h, sigma2hat2_h)))
abline(v = sigma_h^2, col = "red3", lwd = 3)
abline(v = mean(sigma2hat1_h), col = "blue", lty = 2, lwd = 3)
hist(sigma2hat2_h, xlim = range(c(sigma2hat1_h, sigma2hat2_h)))
abline(v = sigma_h^2, col = "red3", lwd = 3)
abline(v = mean(sigma2hat2_h), col = "blue", lty = 2, lwd = 3)


#MOM

#for nd, muhat=mu and sigma^2hat=sigma^2



################################################################################


mu_d = sample_mean_dobermann
sigma_d = sd_dobermann

ks.test(x =  dobermann_rehomed, 
        y = "pnorm", 
        mean = mu_d, 
        sd = sigma_d)
# D = 0.16126, p-value = 0.6457

##

shapiro.test(dobermann_rehomed)
#W = 0.92157, p-value = 0.09327

##

library(nortest)

pearson.test(dobermann_rehomed)
#P = 2.6667, p-value = 0.6151




par(mar = c(2,2,2,2))

mu_s = sample_mean_staffy
sigma_s = sd(staffy_rehomed)

samplemeans <- rep(NA, 100)   # An empty vector to store our 100 sample means.

n <- 205

for(i in 1:100){
  
  simulatedsample <- rnorm(n, mu_s, sigma_s)  # Sample of size n = 10 from Exp(0.08).    
  
  samplemeans[i] <- mean(simulatedsample) # Store the i-th sample mean.
}

# Check if our sample means are normally distributed:

qqnorm(samplemeans)
abline(a = mu_s, b = sigma_s, col = "red")

##

ks.test(x =  staffy_rehomed, 
        y = "pnorm", 
        mean = mu_s, 
        sd = sigma_s)
##

shapiro.test(staffy_rehomed)

##

library(nortest)

pearson.test(staffy_rehomed)

##





samplemeans <- rep(NA, 100)   # An empty vector to store our 100 sample means.

n <- 22

for(i in 1:100){
  
  simulatedsample <- rnorm(n, mu_h, sigma_h)  # Sample of size n = 10 from Exp(0.08).    
  
  samplemeans[i] <- mean(simulatedsample) # Store the i-th sample mean.
}

# Check if our sample means are normally distributed:

qqnorm(samplemeans)

##

ks.test(x =  highland_rehomed, 
        y = "pnorm", 
        mean = mu_h, 
        sd = sigma_h)
##

shapiro.test(highland_rehomed)

##

library(nortest)

pearson.test(highland_rehomed)


ks.test(x =  transformed_h, 
        y = "pnorm", 
        mean = trans_mu, 
        sd = trans_sigma)
##

shapiro.test(transformed_h)

##

library(nortest)

pearson.test(transformed_h)


##

################################################################################


#Inference

#confidence interval
#Dobermann
#z test but using SAMPLE variance (not correct):
# mean = 27 and sigma = 74
mu_d
sigma_d
n = 21

var_d <- sqrt(sigma_d)
var_d

# Calculate the estimate of mu:
xbar_d <- mu_d                           

# Calculate the confidence interval:
CI_d =  xbar_d + c(-1, 1)*1.96*sqrt(var_d^2/n)   

CI_d
#16.46846 19.24583


#Staffy
#z test but using SAMPLE variance (not correct):
mu_s
sigma_s
n <- 205

var_s <- sqrt(sigma_s)


# Calculate the estimate of mu:
xbar_s <- mu_s                             

# Calculate the confidence interval:
CI_s =  xbar_s + c(-1, 1)*1.96*sqrt(var_s^2/n)   

CI_s
#18.83117 19.71955


#Highland

mu_h
sigma_h
n <- 22

var_h <- sqrt(sigma_h)


# Calculate the estimate of mu:
xbar_h <- mu_h            

# Calculate the confidence interval
#z test but using SAMPLE variance (not correct):
CI_h =  xbar_h + c(-1, 1)*1.96*sqrt(var_h^2/n)   

CI_h
#18.89837 21.55617

################################################################################


#dobermann
# Calculate the confidence interval
#t test, mean and var unknown, assume normality and IID:

n = 21


sample_mean_dobermann
sample_variance_dobermann

tscore_d = qt(p = 0.975, df = 20)
#2.056

lower_d = sample_mean_dobermann-(tscore_d*sqrt(sample_variance_dobermann/n))
upper_d = sample_mean_dobermann+(tscore_d*sqrt(sample_variance_dobermann/n))

lower_d
#13.05859

upper_d
#22.65569

t.test(dobermann_rehomed, mu = 27)
#13.05859 22.65569
#t = -3.9745, df = 20, p-value = 0.0007468




#staffy
# Calculate the confidence interval
#t test, mean and var unknown, assume normality and IID:

n = 205

sample_mean_staffy
sample_variance_staffy = var(staffy_rehomed)

tscore_s = qt(p = 0.975, df = 204)
#1.972

tscore_s

lower_s = sample_mean_staffy-(tscore_s*sqrt(sample_variance_staffy/n))
upper_s = sample_mean_staffy+(tscore_s*sqrt(sample_variance_staffy/n))

lower_s
#18.01351

upper_s
#20.91332

t.test(staffy_rehomed, mu = 27)
#18.01351 20.91332
#t = -10.249, df = 204, p-value < 2.2e-16



#staffy z test
w = qnorm(p= 0.975)

n_d = length(staffy_rehomed)

sigma_d = sd(staffy_rehomed)

xbar = mean(staffy_rehomed)


(z_test_d > w) | (z_test_d < -w)


CI=  xbar + c(-1, 1)*1.96*sqrt(sigma_d^2/n_d) 

CI
#15.06367 23.86316

z_test_d

p_value_d = 2*(1-pnorm(z_test_d)) 

p_value_d
#0.0007867942




#highland
# Calculate the confidence interval
#t test, mean and var unknown, assume normality and IID:

n = 22

sample_mean_highland
sample_variance_highland

tscore_h = qt(p = 0.975, df = 21)
#2.080

lower_h = sample_mean_highland-(tscore_h*sqrt(sample_variance_highland/n))
upper_h = sample_mean_highland+(tscore_h*sqrt(sample_variance_highland/n))

lower_h
#15.743

upper_h
#24.711

t.test(highland_rehomed, mu = 27)
#15.74328 24.71127
#t = -3.1411, df = 21, p-value = 0.004932


################################################################################

#comparison 

##seeing of varaince = variance => all TRUE
ds = sd_dobermann/sd_staffy

outcome1 = ifelse(ds >= (1/3) | ds <= 3, "same","not_same")
outcome1


outcome2 = ifelse(ds >= (1/4) | ds <= 4, "same","not_same")
outcome2

dh = sd_dobermann/sd_highland

outcome3 = ifelse(dh >= (1/4) | ds <= 4, "same","not_same")
outcome3

sh = sd_dobermann/sd_highland

outcome3 = ifelse(sh >= (1/3) | ds <= 3, "same","not_same")
outcome3

outcome4 = ifelse(sh >= (1/4) | ds <= 4, "same","not_same")
outcome4

j = length(dobermann_rehomed) 
k = length(staffy_rehomed) 
l = length(highland_rehomed) 

#dobermann vs staffy
meandiff_ds = sample_mean_dobermann - sample_mean_staffy

t_ds = qt(p = 0.975, df = j + k -2)

sp_ds = sqrt(((30-1)*sd_dobermann^2 + (k-1)*sd_staffy^2)/(j + k -2))

meandiff_ds+c(-1,1)*t_ds*sp_ds*sqrt(1/j + 1/k)
#-6.455487  3.242944


t.test(dobermann_rehomed,staffy_rehomed, mu = 0, paired = FALSE, 
       var.equal = TRUE, conf.level = 0.95)
# -6.587800  3.375256
#t = -0.6651, df = 24.272, p-value = 0.5123


#dobermann vs highland 
meandiff_dh = sample_mean_dobermann - sample_mean_highland

t_dh = qt(p = 0.975, df = j + l -2)

sp_dh = sqrt(((30-1)*sd_dobermann^2 + (l-1)*sd_highland^2)/(j + l -2))

meandiff_dh+c(-1,1)*t_dh*sp_dh*sqrt(1/j + 1/l)
# -9.421688  4.681428

t.test(dobermann_rehomed,highland_rehomed, mu = 0, paired = FALSE, 
       var.equal = TRUE, conf.level = 0.95)
#  -8.739112  3.998852
# t = -0.75173, df = 40.677, p-value = 0.4565


#staffy vs highland 
meandiff_sh = sample_mean_staffy - sample_mean_highland

t_sh = qt(p = 0.975, df = k + l -2)

sp_sh = sqrt(((30-1)*sd_staffy^2 + (l-1)*sd_highland^2)/(k + l -2))

sh = meandiff_sh+c(-1,1)*t_sh*sp_sh*sqrt(1/k+ 1/l)
#-2.922184  1.394468

t.test(staffy_rehomed,highland_rehomed, mu = 0, paired = FALSE, 
       var.equal = TRUE, conf.level = 0.95)
#   -5.445438  3.917722
# t = -0.3353, df = 26.133, p-value = 0.7401

p_value_sh = 2*(1-pnorm(abs(-0.32455)))
p_value_sh
sh

par(mfrow = c(2,2))

hist(staffy_rehomed); hist(highland_rehomed); qqnorm(staffy_rehomed);
qqnorm(highland_rehomed)


################################################################################

# Analysis labels for the left side:

analysis = c("Rehomed time for Dobermann (weeks)", 
             "Rehomed time for Staffordshire Bull Terrier (weeks)", 
             "Rehomed time for West Highland White Terrier (weeks)") 

# Results of each test (estimated mean, 
# upper CI limit, lower CI limit, p-value):

estimate  =  c(17.82, 19.46, 20.20)             
upper     =  c(22.66, 20.91, 24.71)
lower     =  c(13.06, 18.01, 15.74)
pval      =  c(2.2e-16,0.0007,0.005)

# Note that the order of the results in each vector
# must match the order of the labels in the 
# vector "analysis".

# Set the margin widths:

par(mar = c(5,5,1,5), mfrow = c(1,1))

# Create an empty plot of a suitable 
# size (considering the width of your
# confidence intervals):

plot(x = 0,                                  # One point at (0,0).
     xlim = c(0, 35), ylim=c(0, 4),        # Axis limits.
     type = "n", xaxt = "n", yaxt="n",       # No points, no axes drawn.
     xlab = NULL, ylab= NULL, ann = FALSE,   # No axis labels or numbers.
     bty="n")                                # No box.

# Add a horizontal (side = 1) axis:

axis(side = 1, cex.axis = 1) 

# Add an axis label 4 lines below the axis:

mtext("Estimated value of the mean across three breeds, with 95% 
      confidence interval", 
      side = 1, line = 3) 

# Add some grid lines, preferably lined up
# with the numbers on the horizontal axis:

for(i in c(5, 15, 25, 35)){
  
  lines(c(i, i), c(0, 4), lty = 2, col = "gray53")
  
}

# Add labels for each analysis on the left (side = 2) 
# at vertical heights of 1, 2, 3 and 4:

verticalpos = 1:3

mtext(text = analysis,  at = verticalpos, 
      side = 2, line = 10, outer = FALSE, las = 1, adj = 0)

# Try changing the "line" option to move these closer to 
# or away from the plotted intervals.

# Plot the four point estimates (centres 
# of the CIs for each analysis): 

points(estimate, verticalpos, pch = 16)


# Plot the four interval estimates:

for(i in 1:3 ){
  
  lines(c(lower[i], upper[i]), c(verticalpos[i], verticalpos[i]))
  
  lines(c(lower[i], lower[i]), c(verticalpos[i] + 0.2, verticalpos[i] - 0.2))
  
  lines(c(upper[i], upper[i]), c(verticalpos[i] + 0.2, verticalpos[i] - 0.2))
  
}

# Now we add numerical results on the right (side = 4), but we 
# need to put them into a nice form first. Note that
# paste() merges words and numbers, and formatC()
# allows us to control the number of decimal places.

est <- formatC(estimate, format='f', digits = 0)

P <- formatC(pval , format = 'f', digits = 3) 

pval <- paste("p =", P)    # Type pval to see what this does.

L <- formatC(lower, format = 'f', digits = 0)
U <- formatC(upper, format = 'f', digits = 0)

interval <- paste("(", L, ", ", U, "),", sep = "")   # Type interval to check.

# Putting it all together:

results <- paste(est, interval, pval)

# Add to the plot:

mtext(text = results, at = verticalpos, 
      side = 4, line = 4, outer = FALSE, las = 1, adj = 1)

# Like a Christmas present, an R 
# plot belongs in a box:

box("inner")
