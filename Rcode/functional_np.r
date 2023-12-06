# data import
dta <- readRDS("data/data_comb.RDS")
names(dta)

# R functions needed
# download from
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/
source("npfda.r")


# the data to use
x <- dta$knee_accl_vt[dta$cond=="slowbw",]
x <- x/10000
tt <- (dta$cycle-1)/100

y <- dta$knee_moment_vt[dta$cond=="slowbw",]
y <- apply(y, 1, max) 


# loo nonparametric regression
n <- length(y)
preds_np <- numeric(n)
for (i in 1:n)
{
  preds_np[i] <- funopare.kernel.cv(y[-i],x[-i,],matrix(x[i,],1,101), q=5,
                                    semimetric = "pca")$Predicted.values
}
mse_np <- mean((y-preds_np)^2)


# compare to pfr/refund
library(refund)
help("pfr")
help("predict.pfr")
preds_pfr <- numeric(n)
for (i in 1:n)
{
  y_i <- y[-i]
  x_i <- x[-i,]
  flmi <- pfr(y_i ~ lf(x_i, k = 20, bs = "ps", m = 2), method = "ML")
  preds_pfr[i] <- predict(flmi, newdata = list("x_i" = matrix(x[i,],1,101)))
}
points(y,preds_pfr, col = 2, pch = 3, lwd =5)
mse_pfr <- mean((y-preds_pfr)^2)



plot(y, preds_pfr, lwd = 5, bty = "n", xlab = "true values", ylab="predicted",
     xlim = c(0.2,0.75), ylim = c(0.2,0.7), col = 2, pch = 3)
abline(c(0,1), lty=2, col = 1, lwd = 2)
points(y,preds_np, col = 4, pch = 1, lwd =5)
legend("bottomright", legend = paste(c("npfda (MSE = ","pfr (MSE = "), round(c(mse_np,mse_pfr), digits = 4), 
                                     c(")", ")"), sep = ""),
       col = c(4,2), pch = c(1,3),
       pt.lwd = 5, bty = "n")




# binary response
x <- dta$hip_accl_ap[dta$cond=="slowbw",]/10000
y <- dta$sex[dta$cond=="slowbw"]
table(y)
y <- as.numeric(as.factor(y))-1
table(y)


# loo nonparametric regression
n <- length(y)
preds_np <- numeric(n)
for (i in 1:n)
{
  preds_np[i] <- funopare.kernel.cv(y[-i],x[-i,],matrix(x[i,],1,101), q=5,
                                    semimetric = "pca")$Predicted.values
}


# compare to pfr/refund
library(refund)
help("pfr")
help("predict.pfr")
preds_pfr <- numeric(n)
for (i in 1:n)
{
  y_i <- y[-i]
  x_i <- x[-i,]
  flmi <- pfr(y_i ~ lf(x_i, k = 20, bs = "ps", m = 2), family = binomial(link = logit),
              method = "ML")
  preds_pfr[i] <- predict(flmi, newdata = list("x_i" = matrix(x[i,],1,101)),
                          type = "response")
}

preds <- c(preds_np,preds_pfr)
yy <- c(y,y)
mthd <- rep(c("npfda","pfr"),c(n,n))

boxplot(preds ~ mthd*yy, xlab = "method : true class",
        ylab = "estimated probability of class 1", col = c(4,2), ylim = c(0,1.05))
abline(v=2.5, lty = 3)
text(c(1.5,3.5),c(1.05,1.05), c("class 0","class 1"))

