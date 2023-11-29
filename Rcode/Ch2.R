### data import ###
#setwd('/Users/sgreven/Dropbox/Publications/Paper/Tutorial_FDA/Bewegungsdaten_Bernard')
dta <- readRDS(file="data_comb.RDS")
x <- dta$knee_accl_vt[dta$cond=="slowbw",] / 10000
y <- x[1,]
tt <- (dta$cycle-1)/100
n <- nrow(x)
L <- length(tt) 

### 2.1 Descriptives ###
library(rainbow)
z <- t(x)
colnames(z) <- 1:ncol(z)
fdsdta <- fds(x = tt, y = z, xname = "tt", yname = "y")
# rainbow plot, red high depth, yellow low depth
pdf(file = "descriptives1.pdf", width = 5, height = 6)
plot.fds(fdsdta, plot.type = "depth", colorchoice = "heat_hcl", lwd = 2, bty = "n",
         xlab= 'relative time', ylab = 'knee axial rotation acceleration / 10000') 
dev.off()
pdf(file = "descriptives2.pdf", width = 5, height = 6)
# functional boxplot
par(bty = "n")
fboxplot(fdsdta,  plot.type = "functional", legendpos = "bottomright", lwd = 2, 
         xlab= 'relative time', ylab = 'knee axial rotation acceleration / 10000')
dev.off()


### 2.2 smoothing of first observation ###
library(fda)
nbasis = 10
basisobj = create.bspline.basis(c(0,1),nbasis)
pdf("Bsplinebasis.pdf", width = 5, height = 6)
  plot(basisobj, lwd = 2, bty = "n", xlab = "relative time", ylab = "B-spline basis")
dev.off()
nbasis = 30
basisobj = create.bspline.basis(c(0,1),nbasis)
ys = smooth.basis(argvals=tt, y=y, fdParobj=basisobj)
pdf("fittedcurve.pdf", width = 5, height = 6)
  plotfit.fd(y, tt, ys$fd, main = "", lwd = 2, bty = "n", xlab = "relative time", ylab = "reconstruction of knee axial rotation acceleration / 10000") 
dev.off()


### 2.3 Phase and Amplitude Variation ###
library(fdasrvf)
library(rainbow)
alignedknee <- time_warping(f = t(x), time = tt)

pdf(file = "unaligned.pdf", width = 5, height = 6)
  fdsdta <- fds(x = tt, y = alignedknee$f0, xname = "tt", yname = "y")
  plot.fds(fdsdta, plot.type = "functions", cex=1.5, cex.lab = 1.5, lwd = 2, bty = "n", #colorchoice = "heat_hcl",
         xlab= 'relative time', ylab = 'knee axial rotation acceleration / 10000')
  unalignedmean <- apply(alignedknee$f0, 1, mean)
  lines(tt, unalignedmean, lwd =6, lty = 1)
dev.off()
pdf(file = "aligned.pdf", width = 5, height = 6)
  fdsdta <- fds(x = tt, y = alignedknee$fn, xname = "tt", yname = "y")
  plot.fds(fdsdta, plot.type = "functions", cex=1.5, cex.lab = 1.5, lwd = 2, bty = "n",
         xlab= 'aligned time', ylab = 'knee axial rotation acceleration / 10000')
  lines(tt, alignedknee$fmean, lwd =6, lty = 1)
dev.off()
pdf(file = "warping.pdf", width = 5, height = 6)
  fdsdta <- fds(x = tt, y = alignedknee$warping_functions, xname = "tt", yname = "y")
  plot.fds(fdsdta, plot.type = "functions", cex=1.5, cex.lab = 1.5, lwd = 2, bty = "n",
         xlab = 'aligned time', ylab= 'relative time')
  lines(tt, tt, lwd =6, lty = 1)
dev.off()


### 2.4 FPCA ###
library(refund)
myfpcs <- fpca.face(Y = x, pve = 0.99, p = 3, m = 1) #,  argvals = tt

# Eigenfunctions
phi <- myfpcs$efunctions
phi <- phi * sqrt(L) # make eigenfunctions orthonormal w.r.t. the L2 and not vector inner product
# Eigenfunctions orthonormal? yes, matrix of inner products yields identity matrix
round(t(phi) %*% phi / L,6)

# Eigenvalues
print('explained variance')
nu <- myfpcs$evalues / L # rescale eigenvalues/variances according to rescaling of eigenfunctions
myfracs <- round(nu / sum(nu),2) * 100
plot(nu)

pdf("FPCA.pdf")
par(mfrow = c(2,2))
mu <- apply(x, 2, mean)
mylims <- range(cbind(mu + sqrt(nu[1]) * phi[,1], mu + sqrt(nu[2]) * phi[,2]))
plot(tt, mu, ylim = mylims, type = "l", xlab = 'relative time', 
     ylab = 'Mean +/- multiple of 1st Eigenfunction', main = "FPC1 (41% Variance)", lwd = 2, bty = "n") 
points(tt, mu + sqrt(nu[1]) * phi[,1], pch = "+", col = 2)
points(tt, mu - sqrt(nu[1]) * phi[,1], pch = "-", col = 2)
plot(tt, mu, ylim = mylims, type = "l", xlab = 'relative time', 
     ylab = 'Mean +/- multiple of 2nd Eigenfunction', main = "FPC2 (21% Variance)", lwd = 2, bty = "n") 
points(tt, mu + sqrt(nu[2]) * phi[,2], pch = "+", col = 2)
points(tt, mu - sqrt(nu[2]) * phi[,2], pch = "-", col = 2)
plot(tt, tt, ylim = c(-3, max(phi[,1:2])), col=0, xlab = 'relative time', 
     ylab = 'Eigenfunctions', main = "First four Eigenfunctions", bty = "n") 
for (k in 1:4){
  lines(tt, phi[,k], col = k)
}
abline(h=0, lty = 2)
legend('bottomright', col = 1:4, lty = 1, 
       legend = paste(myfracs[1:4], "%"), ncol = 2, bty = "n", 
       title = 'Eigenfunctions (%Var.)') 

# Scores
sex <- ifelse(dta$sex == 'f', 1, 2) # 1 female, 2 male
xi <- myfpcs$scores / sqrt(L) # rescale scores according to rescaling of eigenfunctions
plot(xi[,1], xi[,2], col = sex, pch = sex, xlab = 'FPC score 1', 
     ylab = 'FPC score 2', #ylim = c(-6,7), 
     main = "First two FPC scores", bty = "n")
legend('topright', col = 1:2, pch = 1:2, legend = c('female','male')) 
dev.off()

# check that fitted functions with 10 FPCs closely match the observed ones
par(mfrow = c(1,2))
myfit <- mu + phi %*% t(xi)
fdsdta2 <- fds(x = tt, y = myfit, xname = "tt", yname = "y")
plot.fds(fdsdta2, plot.type = "functions", colorchoice = "heat_hcl", lwd = 2, bty = "n",
         xlab= 'relative time', ylab = 'knee axial rotation acceleration / 10000', ylim = c(-2.5,2.5)) 
plot.fds(fdsdta, plot.type = "functions", colorchoice = "heat_hcl", lwd = 2, bty = "n",
         xlab= 'relative time', ylab = 'knee axial rotation acceleration / 10000', ylim = c(-2.5,2.5)) 

