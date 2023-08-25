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
matplot(tt,t(x), type="l", lty = 1, col = cx[cy], xlab = "relative time",
        ylab = "knee axial rotation acceleration / 10000", bty = "n", lwd = 2)
points((41:100)/101, rep(-2,60), pch = 15, col = rep(cx,rep(10,6)))
text((4:10)/10, rep(-2.2,7), seq(0.2,0.8,by=0.1), cex = 1)
dev.off()


pdf("kav_d.pdf", width = 4, height = 6)
matplot(tt,t(x), type="l", lty = 1, col = cx[cy], xlab = "Relative Zeit",
        ylab = "Vertikale Kniebeschleunigung", bty = "n", lwd = 2)
points((41:100)/101, rep(-2,60), pch = 15, col = rep(cx,rep(10,6)))
text((4:10)/10, rep(-2.2,7), seq(0.2,0.8,by=0.1), cex = 0.8)
dev.off()




# model fitting
library(refund)
help("pfr")
flm1 <- pfr(y ~ lf(x, k = 15, bs = "ps", m=2) + sx, method = "REML")
summary(flm1)

# estimated beta function
pdf("kavb151.pdf", width = 5, height = 6)
plot(flm1, rug = FALSE, xlab = "percentage cycle points", main = "Response: maximum moment", shade=TRUE,
     ylab = "estimated coefficient function", bty = "n", lwd = 2, ylim = c(-2.2,2))
abline(h = 0, lty=3)
dev.off()


pdf("kavb_d.pdf", width = 4, height = 6)
plot(flm1, rug = FALSE, xlab = "Relative Zeit", main = "Response: Maximales Moment", shade=FALSE,
     ylab = "Geschätzte Koeffizientenfunktion", bty = "n", lwd = 2, col = 4)
abline(h = 0, lty=3)
dev.off()


# check alpha level (approximate significance of smooth terms)
R <- 1000
pvals <- numeric(R)
for (r in 1:R)
{
  y0 <- sample(y, 31, replace = TRUE)
  flmr <- pfr(y0 ~ lf(x, k = 20) + sx, method = "REML")
  pvals[r] <- summary(flmr)$s.table[1,4]
}

mean(pvals < 0.05)
hist(pvals)





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
matplot(tt,t(x), type="l", lty = y+1, col = -3*y + 4, xlab = "relative time",
        ylab = "hip abduction-adduction acceleration / 10000", bty = "n", lwd = 2)
legend("topright", legend = c("female","male"), col = c(4,1), lty = 1:2, lwd = 2,
       bty = "n")
dev.off()


pdf("haa_d.pdf", width = 4, height = 6)
matplot(tt,t(x), type="l", lty = y+1, col = -2*y + 3, xlab = "Relative Zeit",
        ylab = "Hüftbeschleunigung (anterior/posterior)", bty = "n", lwd = 2)
legend("topright", legend = c("Frauen","Männer"), col = c(3,1), lty = 1:2, lwd = 2,
       bty = "n")
dev.off()



# model fitting
flm2 <- pfr(y ~ lf(x, k = 15, bs = "ps", m = 2), family = binomial(link = logit),
            method = "REML")
summary(flm2)


# estimated beta
pdf("haab152.pdf", width = 5, height = 6)
plot(flm2, rug = FALSE, xlab = "percentage cycle points", main = "Response: sex", shade=TRUE, 
     ylab = "estimated coefficient function", bty = "n", lwd = 2, ylim = c(-100,250))
abline(h = 0, lty=3)
dev.off()


pdf("haab_d.pdf", width = 4, height = 6)
plot(flm2, rug = FALSE, xlab = "Relative Zeit", main = "Response: Geschlecht", shade=FALSE,
     ylab = "Geschätzte Koeffizientenfunktion", bty = "n", lwd = 2, col = 2)
abline(h = 0, lty=3)
dev.off()


# check alpha level (approximate significance of smooth terms)
R <- 1000
pvals <- numeric(R)
for (r in 1:R)
{
  y0 <- sample(c(0,1), 31, replace = TRUE, prob = c(1-mean(y),mean(y)))
  flmr <- pfr(y0 ~ lf(x, k = 20), family = binomial(link = logit))
  pvals[r] <- summary(flmr)$s.table[1,4]
}

mean(pvals < 0.05)
hist(pvals)



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
summary(m0)
plot(m0)

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
summary(m1)
plot(m1)
t1 <- predict(m1, type="terms", se = TRUE)
names(t1)
names(t1$fit)
names(t1$se.fit)

head(t1$fit[[1]])
head(t1$fit[[2]])
head(t1$se.fit[[1]])
head(t1$se.fit[[2]])


# intercept
pdf("fosra.pdf", width = 5, height = 6)
plot(tt, t1$fit[[1]][1,], type = "l", xlab = "t",
     ylab = expression(alpha(t)), lwd = 2, bty = "n", ylim = c(-0.5, 0.7),
     col = 4)
lines(tt, t1$fit[[1]][1,] + 2*t1$se.fit[[1]][1,], type = "l", lty = 3, col = 4)
lines(tt, t1$fit[[1]][1,] - 2*t1$se.fit[[1]][1,], type = "l", lty = 3, col = 4)
dev.off()


pdf("fosra_d.pdf", width = 4, height = 6)
plot(tt, t1$fit[[1]][1,], type = "l", xlab = "t",
     ylab = expression(alpha(t)), lwd = 2, bty = "n", ylim = c(-0.5, 0.7),
     col = 3)
lines(tt, t1$fit[[1]][1,] + 2*t1$se.fit[[1]][1,], type = "l", lty = 3, col = 4)
lines(tt, t1$fit[[1]][1,] - 2*t1$se.fit[[1]][1,], type = "l", lty = 3, col = 4)
dev.off()


# beta
pdf("fosrb.pdf", width = 5, height = 6)
plot(tt, t1$fit[[2]][1,], type = "l", xlab = "t",
     ylab = expression(beta(t)), ylim = c(-0.1,0.2), lty = 2, lwd = 2, bty = "n")
lines(tt, t1$fit[[2]][1,] + 2*t1$se.fit[[2]][1,], type = "l", lty = 3)
lines(tt, t1$fit[[2]][1,] - 2*t1$se.fit[[2]][1,], type = "l", lty = 3)

lines(tt, t1$fit[[2]][2,], col = 4, lwd = 2)
lines(tt, t1$fit[[2]][2,] + 2*t1$se.fit[[2]][2,], type = "l", col = 4, lty = 3)
lines(tt, t1$fit[[2]][2,] - 2*t1$se.fit[[2]][2,], type = "l", col = 4, lty = 3)

legend("topright", legend = c("female","male"), col = c(4,1), lty = 1:2, lwd = 2,
       bty = "n")
dev.off()



pdf("fosrb_d.pdf", width = 4, height = 6)
plot(tt, t1$fit[[2]][1,], type = "l", xlab = "t",
     ylab = expression(beta(t)), ylim = c(-0.1,0.2), lty = 2, lwd = 2, bty = "n")
lines(tt, t1$fit[[2]][1,] + 2*t1$se.fit[[2]][1,], type = "l", lty = 3)
lines(tt, t1$fit[[2]][1,] - 2*t1$se.fit[[2]][1,], type = "l", lty = 3)

lines(tt, t1$fit[[2]][2,], col = 3, lwd = 2)
lines(tt, t1$fit[[2]][2,] + 2*t1$se.fit[[2]][2,], type = "l", col = 3, lty = 3)
lines(tt, t1$fit[[2]][2,] - 2*t1$se.fit[[2]][2,], type = "l", col = 3, lty = 3)

legend("topright", legend = c("Frauen","Männer"), col = c(3,1), lty = 1:2, lwd = 2,
       bty = "n")
dev.off()


pdf("fosrre.pdf", width = 5, height = 6)
matplot(tt, t(t1$fit[[3]]), type = "l", col = -3*x + 4, xlab = "t", lty = 1, bty = "n",
        ylab = expression(E[i](t)))
dev.off()


# white noise error? not really...
matplot(t(resid(m1)), type = "l",  col = -3*x + 4)

