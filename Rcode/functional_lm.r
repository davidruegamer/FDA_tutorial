# data import
dta <- readRDS("data_comb.RDS")
names(dta)


# data for the functional linear model

# functional predictor
# make sure that there is only one observation per subject
x <- dta$knee_accl_vt[dta$cond=="slowbw",]
x <- x/10000
tt <- (dta$cycle-1)/100

# scalar response
y <- dta$knee_moment_vt[dta$cond=="slowbw",]
y <- apply(y, 1, max) 

# scalar covariate
sx <- dta$sex[dta$cond=="slowbw"]

# plot the functional data with colors corresponding to y
cx <- rainbow(12)[7:12]
cy <- floor(10*y)-1
pdf("kav.pdf", width = 5, height = 6)
matplot(tt,t(x), type="l", lty = 1, col = cx[cy], xlab = "percentage cycle points",
        ylab = "knee acceleration vertical/10000", bty = "n", lwd = 2)
points((41:100)/101, rep(-2,60), pch = 15, col = rep(cx,rep(10,6)))
text((4:10)/10, rep(-2.2,7), seq(0.2,0.8,by=0.1), cex = 1)
dev.off()



# model fitting
library(refund)
help("pfr")
flm1 <- pfr(y ~ lf(x, k = 20) + sx, method = "REML")
summary(flm1)

# estimated beta function
pdf("kavb.pdf", width = 5, height = 6)
plot(flm1, rug = FALSE, xlab = "percentage cycle points", shade=TRUE,
     ylab = "estimated coefficient function", bty = "n", lwd = 2)
abline(h = 0, lty=3)
dev.off()



# data for the functional logit model

# change functional predictor
x <- dta$hip_accl_ap[dta$cond=="slowbw",]/10000

# use sex as response
y <- dta$sex[dta$cond=="slowbw"]
table(y)
y <- as.numeric(as.factor(y))-1
table(y)


# plot the functional data with colors corresponding to y
pdf("haa.pdf", width = 5, height = 6)
matplot(tt,t(x), type="l", lty = y+1, col = -3*y + 4, xlab = "percentage cycle points",
        ylab = "hip acceleration anterior-posterior/10000", bty = "n", lwd = 2)
legend("topright", legend = c("female","male"), col = c(4,1), lty = 1:2, lwd = 2, bty = "n")
dev.off()


# model fitting
flm2 <- pfr(y ~ lf(x, k = 20), family = binomial(link = logit))
summary(flm2)


# estimated beta
pdf("haab.pdf", width = 5, height = 6)
plot(flm2, rug = FALSE, xlab = "percentage cycle points", shade=TRUE,
     ylab = "estimated coefficient function", bty = "n", lwd = 2)
abline(h = 0, lty=3)
dev.off()

