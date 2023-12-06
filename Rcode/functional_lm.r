# data import
dta <- readRDS("data/data_comb.RDS")
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
matplot(tt,t(x), type="l", lty = 1, col = cx[cy], xlab = "relative time",
        ylab = "knee axial rotation acceleration / 10000", bty = "n", lwd = 2)
points((41:100)/101, rep(-2,60), pch = 15, col = rep(cx,rep(10,6)))
text((4:10)/10, rep(-2.2,7), seq(0.2,0.8,by=0.1), cex = 1)


# model fitting
library(refund)
help("pfr")

# functional linear model
flm1 <- pfr(y ~ lf(x, k = 15, bs = "ps", m=2) + sx, method = "REML")
summary(flm1)

# estimated beta function
plot(flm1, rug = FALSE, xlab = "relative time", main = "Response: maximum moment", shade=TRUE,
     ylab = "estimated coefficient function", bty = "n", lwd = 2, ylim = c(-2.2,2))
abline(h = 0, lty=3)




# data for the functional logit model

# change functional predictor
x <- dta$hip_accl_ap[dta$cond=="slowbw",]/10000

# use sex as response
y <- dta$sex[dta$cond=="slowbw"]
table(y)
y <- as.numeric(as.factor(y))-1
table(y)


# plot the functional data with colors corresponding to y
matplot(tt,t(x), type="l", lty = y+1, col = -3*y + 4, xlab = "relative time",
        ylab = "hip abduction-adduction acceleration / 10000", bty = "n", lwd = 2)
legend("topright", legend = c("female","male"), col = c(4,1), lty = 1:2, lwd = 2,
       bty = "n")



# model fitting
flm2 <- pfr(y ~ lf(x, k = 15, bs = "ps", m = 2), family = binomial(link = logit),
            method = "REML")
summary(flm2)


# estimated beta
plot(flm2, rug = FALSE, xlab = "relative time", main = "Response: sex", shade=TRUE, 
     ylab = "estimated coefficient function", bty = "n", lwd = 2, ylim = c(-100,250))
abline(h = 0, lty=3)




# function-on-scalar regression

# switch x and y
y <- dta$hip_accl_ap[dta$cond=="slowbw",]/10000
x <- dta$sex[dta$cond=="slowbw"]
table(x)
x <- as.numeric(as.factor(x)) -1
table(x)


# compare pffr manual/examples
help("pffr")

# fit model under independence assumption:
m0 <- pffr(y ~ x, yind=tt)

# get first eigenfunctions of residual covariance 
# (i.e. first functional PCs of empirical residual process
# with 95% variance explained)
help("fpca.sc")
Ehat <- resid(m0)
fpcE <- fpca.sc(Ehat, nbasis = 20, pve = 0.95)
efunctions <- fpcE$efunctions
evalues <- fpcE$evalues
id <- factor(1:nrow(y))

# refit model with fpc-based residuals
m1 <- pffr(y ~ x + pcre(id=id, efunctions=efunctions, evalues=evalues, yind=tt),
           method = "ML",
           bs.yindex = list(bs="ps", k=20, m=c(2, 1)),
           bs.int = list(bs="ps", k=25, m=c(2, 1)),
           yind=tt)
t1 <- predict(m1, type="terms", se = TRUE)



# intercept
plot(tt, t1$fit[[1]][1,], type = "l", xlab = "relative time t",
     ylab = expression(alpha(t)), lwd = 2, bty = "n", ylim = c(-0.5, 0.7),
     col = 4)
lines(tt, t1$fit[[1]][1,] + 2*t1$se.fit[[1]][1,], type = "l", lty = 3, col = 4)
lines(tt, t1$fit[[1]][1,] - 2*t1$se.fit[[1]][1,], type = "l", lty = 3, col = 4)


# beta
plot(tt, t1$fit[[2]][1,], type = "l", xlab = "relative time t",
     ylab = expression(beta(t)), ylim = c(-0.1,0.2), lty = 2, lwd = 2, bty = "n")
lines(tt, t1$fit[[2]][1,] + 2*t1$se.fit[[2]][1,], type = "l", lty = 3)
lines(tt, t1$fit[[2]][1,] - 2*t1$se.fit[[2]][1,], type = "l", lty = 3)

lines(tt, t1$fit[[2]][2,], col = 4, lwd = 2)

legend("topright", legend = c("female","male"), col = c(4,1), lty = 1:2, lwd = 2,
       bty = "n")



