### Section 4: Statistical Inference with Functional Data ###


# data import
dta <- readRDS("data/data_comb.RDS")
names(dta)



### 4.1 Functional ANOVA ###

# data for functional ANOVA
x_ankle <- dta[9:20]
x_hip <- dta[24:35]
x_knee <- dta[36:47]
x_conds <- c("slowbw","slowten","slowtwe")
tt <- (dta$cycle-1)/100

conds <- dta$cond[dta$sex=="m" & dta$cond%in%x_conds]
ids <- dta$id[dta$sex=="m" & dta$cond%in%x_conds]
whi <- c(1:5,21:25,41:45)
whc <- conds[whi]
cds <- unique(whc)

dta$cond[dta$study=="bernard_study1" & dta$sex=="m"]
dta$id[dta$study=="bernard_study1" & dta$sex=="m"]
dta$side[dta$study=="bernard_study1" & dta$sex=="m"]


# Figure 9

titles <- expand.grid(c("abd.-add.", "flex.-ext.","axial rot."),c("accel.","angle","moment","vel."))
titles <- paste("(", letters[1:12], ") ", titles$Var2, " ", titles$Var1, sep = "")
par(mfrow = c(4,3))
x <- x_knee
for (j in 1:length(x))
{
  xj <- x[[j]]
  xj <- xj[dta$sex=="m" & dta$cond%in%x_conds,]
  xj <- xj[whi,]
  matplot(tt, t(xj), type="l", col = as.numeric(as.factor(whc)),
          lty = 3,
          xlab = "relative time", ylab = "", bty = "n")
  title(titles[j])
  for (k in 1:length(cds)) {
    xjk <- xj[whc==cds[k],]
    lines(tt, apply(xjk,2,mean), lwd = 2, col = k)
  }
}
par(mfrow = c(1,1))



# functional ANOVA
library(fdANOVA)
help("fanova.tests")


# go through functional variables
x <- x_knee
grps <- factor(whc)
pvalsFP <- numeric(length(x))
names(pvalsFP) <- names(x)
pvalsGPF <- numeric(length(x))
names(pvalsGPF) <- names(x)
pvalsFmaxb <- numeric(length(x))
names(pvalsFmaxb) <- names(x)
pvalsCH <- numeric(length(x))
names(pvalsCH) <- names(x)
pvalsTRP <- numeric(length(x))
names(pvalsTRP) <- names(x)

set.seed(1234)
for (j in 1:length(x)) {
  print(j)
  xj <- x[[j]]
  xj <- xj[dta$sex=="m" & dta$cond%in%x_conds,]
  xj <- xj[whi,]
  faovj <- fanova.tests(x = t(xj), group.label = grps,
                        test = c("FP", "GPF", "Fmaxb", "CH", "TRP"))
  pvalsFP[j] <- faovj$FP[[2]]
  pvalsGPF[j] <- faovj$GPF[[2]]
  pvalsFmaxb[j] <- faovj$Fmaxb[[2]]
  pvalsCH[j] <- faovj$CH[[2]]
  pvalsTRP[j] <- faovj$TRP[[2]]
}


# p values
pvals <- cbind(pvalsFP,pvalsGPF,pvalsFmaxb,pvalsCH,pvalsTRP)


# Figure 10

par(mfrow = c(1,1))
barplot(t(pvals), beside = TRUE, ylim = c(0,1.1), col = c(1:4,7), ylab = "p-value", xlab = "variable",
        names.arg = paste("(", letters[1:12], ")", sep = ""),
        legend.text = c("Gorecki and Smaga (2015)", "Zhang and Liang (2014)", "Zhang et al. (2019)", 
                        "Cuevas et al. (2004)", "Cuesta-A. and Febrero-B. (2010)"))
abline(h=c(0.05,0.1), lty=2:3)




### DO NOT RUN: Supplementary Material ###

# check TRP
R <- 100
pvalsTRP_anova <- numeric(R)
pvalsTRP_ats <- numeric(R)
pvalsTRP_wtps <- numeric(R)
j <- 3

set.seed(2022)
for (r in 1:R) {
  print(r)
  xj <- x[[j]]
  xj <- xj[dta$sex=="m" & dta$cond%in%x_conds,]
  xj <- xj[whi,]
  faovr <- fanova.tests(x = t(xj), group.label = grps, test = c("TRP"))
  pvalsTRP_anova[r] <- faovr$TRP[[1]]
  pvalsTRP_ats[r] <- faovr$TRP[[2]]
  pvalsTRP_wtps[r] <- faovr$TRP[[3]]
}

hist(pvalsTRP_anova)
hist(pvalsTRP_ats)
hist(pvalsTRP_wtps)

alpha <- 0.05
mean(pvalsTRP_anova < alpha)
mean(pvalsTRP_ats < alpha)
mean(pvalsTRP_wtps < alpha)




# check type-I error rates

x <- x_knee$knee_angle_ap
dim(x)
grp <- sample(rep(1:5,nrow(x)/5), nrow(x), replace = FALSE)
table(grp)

grp <- sort(rep(1:5,30))
xr <- x[sample(1:nrow(x),length(grp)),]

system.time(faov <- fanova.tests(x = t(xr), group.label = grp,
                                 test = c("FP", "GPF", "Fmaxb", "CH", "TRP")))

R <- 1000
pvalsFP <- numeric(R)
pvalsGPF <- numeric(R)
pvalsFmaxb <- numeric(R)
pvalsCH <- numeric(R)
pvalsTRP_anova <- numeric(R)
pvalsTRP_ats <- numeric(R)
pvalsTRP_wtps <- numeric(R)

set.seed(2022)
for (r in 1:R) {
  print(r)
  xr <- x[sample(1:nrow(x),length(grp)),]
  faovr <- fanova.tests(x = t(xr), group.label = grp,
                        test = c("FP", "GPF", "Fmaxb", "CH", "TRP"))
  pvalsFP[r] <- faovr$FP[[2]]
  pvalsGPF[r] <- faovr$GPF[[2]]
  pvalsFmaxb[r] <- faovr$Fmaxb[[2]]
  pvalsCH[r] <- faovr$CH[[2]]
  pvalsTRP_anova[r] <- faovr$TRP[[1]]
  pvalsTRP_ats[r] <- faovr$TRP[[2]]
  pvalsTRP_wtps[r] <- faovr$TRP[[3]]
}


# type-I error rates
alpha <- 0.1
mean(pvalsFP < alpha)
mean(pvalsGPF < alpha)
mean(pvalsFmaxb < alpha)
mean(pvalsCH < alpha)
mean(pvalsTRP_anova < alpha)
mean(pvalsTRP_ats < alpha)
mean(pvalsTRP_wtps < alpha)


# declare significance if any p-value is below alpha
pvals <- cbind(pvalsFP,pvalsGPF,pvalsFmaxb,pvalsCH,pvalsTRP_wtps)
mean(apply(pvals,1,min) < alpha)



# TRP "with cheating"

R <- 100
pvalsTRP_anova3 <- numeric(R)
pvalsTRP_ats3 <- numeric(R)
pvalsTRP_wtps3 <- numeric(R)
pvalsTRP_anova5 <- numeric(R)
pvalsTRP_ats5 <- numeric(R)
pvalsTRP_wtps5 <- numeric(R)


set.seed(2022)
for (r in 1:R) {
  print(r)
  xr <- x[sample(1:nrow(x),length(grp)),]
  pvalsTRP_aovr <- numeric(5)
  pvalsTRP_atsr <- numeric(5)
  pvalsTRP_wtpsr <- numeric(5)
  for (j in 1:5) {
    faovrj <- fanova.tests(x = t(xr), group.label = grp, test = c("TRP"))
    pvalsTRP_aovr[j] <- faovrj$TRP[[1]]
    pvalsTRP_atsr[j] <- faovrj$TRP[[2]]
    pvalsTRP_wtpsr[j] <- faovrj$TRP[[3]]
  }
  pvalsTRP_anova3[r] <- min(pvalsTRP_aovr[1:3])
  pvalsTRP_ats3[r] <- min(pvalsTRP_atsr[1:3])
  pvalsTRP_wtps3[r] <- min(pvalsTRP_wtpsr[1:3])
  pvalsTRP_anova5[r] <- min(pvalsTRP_aovr)
  pvalsTRP_ats5[r] <- min(pvalsTRP_atsr)
  pvalsTRP_wtps5[r] <- min(pvalsTRP_wtpsr)
}

# type-I error rates
alpha <- 0.1
mean(pvalsTRP_anova3 < alpha)
mean(pvalsTRP_ats3 < alpha)
mean(pvalsTRP_wtps3 < alpha)
mean(pvalsTRP_anova5 < alpha)
mean(pvalsTRP_ats5 < alpha)
mean(pvalsTRP_wtps5 < alpha)
