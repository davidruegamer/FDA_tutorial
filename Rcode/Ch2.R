### data import ###
setwd('/Users/sgreven/Dropbox/Publications/Paper/Tutorial_FDA/Bewegungsdaten_Bernard')
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dta <- readRDS(file="data_comb.RDS")
#names(dta)
#x <- dta$knee_accl_vt[dta$cond=="slowbw",] / 10000
x <- dta$knee_accl_vt[dta$cond=="slowbw",] / 10000
y <- x[1,]
tt <- (dta$cycle-1)/100
n <- nrow(x)

### smoothing of first observation ###
library(mgcv)
oneobs <- gam(y ~ s(tt, k = 30))

plot(oneobs, resid = TRUE, xlab= 'relative time', 
     ylab = 'knee acceleration vertical')
lines(tt, y - oneobs$coefficients[1], col= 2)


library(fda)
nbasis = 10
basisobj = create.bspline.basis(c(0,1),nbasis)
pdf("Bsplinebasis.pdf")
plot(basisobj)
dev.off()
nbasis = 30
basisobj = create.bspline.basis(c(0,1),nbasis)
ys = smooth.basis(argvals=tt, y=y, fdParobj=basisobj)
pdf("fittedcurve.pdf")
plotfit.fd(y, tt, ys$fd, main = "")
dev.off()

### FPCA ###
library(refund)
myfpcs <- fpca.face(Y = x, argvals = tt, pve = 0.99, p = 3, m = 1)
# Eigenvalues
print('explained variance')
round(myfpcs$evalues / sum(myfpcs$evalues),2)
plot(myfpcs$evalues)
# Eigenfunctions

par(mfrow = c(2,2))
mu <- apply(x, 2, mean)
plot(tt, mu, ylim = range(mu + myfpcs$efunctions), type = "l", xlab = 'relative time', 
     ylab = 'Mean +/- 1st Eigenfunction', main = "FPC1 (41% Variance)")
points(tt, mu + myfpcs$efunctions[,1], pch = "+", col = 2)
points(tt, mu - myfpcs$efunctions[,1], pch = "-", col = 2)
plot(tt, mu, ylim = range(mu + myfpcs$efunctions), type = "l", xlab = 'relative time', 
     ylab = 'Mean +/- 2nd Eigenfunction', main = "FPC2 (21% Variance)")
points(tt, mu + myfpcs$efunctions[,2], pch = "+", col = 2)
points(tt, mu - myfpcs$efunctions[,2], pch = "-", col = 2)
plot(tt, tt, ylim = range(myfpcs$efunctions), col=0, xlab = 'relative time', 
     ylab = 'Eigenfunctions', main = "First four Eigenfunctions")
for (k in 1:4){
  lines(tt, myfpcs$efunctions[,k], col = k)
}
#legend('bottomright', col = 1:4, lty = 1, legend = 1:4) #, title = 'Eigenfunctions')

# Scores
xi <- myfpcs$scores
sex <- ifelse(dta$sex == 'f', 1, 2) # 1 female, 2 male
plot(xi[,1], xi[,2], col = sex, pch = sex, xlab = 'FPC score 1', ylab = 'FPC score 2',
     main = "First two FPC scores")
legend('bottomright', col = 1:2, pch = 1:2, legend = c('f','m'))



### Descriptives ###
library(rainbow)
z <- t(x)
colnames(z) <- 1:ncol(z)
fdsdta <- fds(x = tt, y = z, xname = "tt", yname = "y")
# rainbow plot, red high depth, yellow low depth
pdf(file = "descriptives1.pdf")
plot.fds(fdsdta, plot.type = "depth", colorchoice = "heat.colors", 
         xlab= 'relative time', 
         ylab = 'knee acceleration vertical')
dev.off()
pdf(file = "descriptives2.pdf")
# functional boxplot
fboxplot(fdsdta,  plot.type = "functional", legendpos = "bottomright",
         xlab= 'relative time', 
         ylab = 'knee acceleration vertical')
dev.off()

par(mfrow = c(1,1))

### Phase and Amplitude Variation ###
library(fda)
nbasis = 30
basisobj = create.bspline.basis(c(0,1),nbasis)
ys = smooth.basis(argvals=tt, y=t(x), fdParobj=basisobj)
plot(ys)

yfd <- as.fd(ys)


Wnbasis   <- 6
Wbasis    <- create.bspline.basis(c(0,1), Wnbasis)
Wfd0      <- fd(matrix(0,Wnbasis,1),Wbasis)
WfdParobj <- fdPar(Wfd0, Lfdobj=2, lambda=0.01)

# smooth alignment
yalign <- register.fd(yfd, crit = 2, WfdParobj=WfdParobj)
#plotreg.fd(yalign)
plot(yfd)
plot(yalign$regfd)
plot(yalign$warpfd)

#landmark registration
diff <- function(x){
  mylen <- length(x)
  x[2:mylen] - x[1:(mylen-1)]
}

plot(x[1,])
plot(diff(x[1,]))

landm <- function(x){
  mylen <- length(x)
  minmax <- c()
  for (k in 3:(mylen-2)){
    if (x[k] > x[k-1] && x[k-1] > x[k-2] && x[k] > x[k+1] && x[k+1] > x[k+2]){
      minmax <- c(minmax, tt[k])
    }
    if (x[k] < x[k-1] && x[k-1] < x[k-2] && x[k] < x[k+1] && x[k+1] < x[k+2]){
      minmax <- c(minmax, tt[k])
    }
  }
  if (minmax[1] > 0.12){
    minmax <- c(0.01,minmax)
  }
  if (max(minmax) < 0.82){
    minmax <- c(minmax,0.99)
  }
  return(minmax)
}

myind <- c(1:2, 5, 7:11, 13:14, 17:18, 20, 22:30)
myx <- x[myind,]
ximarks <- t(apply(myx,1,landm))
x0marks <- apply(ximarks,2,mean)
nbasis = 30
basisobj = create.bspline.basis(c(0,1),nbasis)
ys = smooth.basis(argvals=tt, y=t(myx), fdParobj=basisobj)
plot(ys)
yfd <- as.fd(ys)


lmkreg <- landmarkreg(yfd, ximarks, x0marks)





