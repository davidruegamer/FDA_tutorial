### Section 1: Introduction ###


# data import
dta <- readRDS("data/data_comb.RDS")
names(dta)


### 1.3 Running Data Example(s) ###


# knee axial rotation acceleration
# (make sure that there is only one observation per subject)
x <- dta$knee_accl_vt[dta$cond=="slowbw",]
x <- x/10000
tt <- (dta$cycle-1)/100

# moment
y <- dta$knee_moment_vt[dta$cond=="slowbw",]

# maximum
y <- apply(y, 1, max) 


# Figure 1 (left)
# plot x with colors corresponding to y
cx <- rainbow(12)[7:12]
cy <- floor(10*y)-1
matplot(tt,t(x), type="l", lty = 1, col = cx[cy], xlab = "relative time",
        ylab = "knee axial rotation acceleration / 10000", bty = "n", lwd = 2)
points((41:100)/101, rep(-2,60), pch = 15, col = rep(cx,rep(10,6)))
text((4:10)/10, rep(-2.2,7), seq(0.2,0.8,by=0.1), cex = 1)



# hip abduction–adduction acceleration
x <- dta$hip_accl_ap[dta$cond=="slowbw",]/10000

# sex
y <- dta$sex[dta$cond=="slowbw"]
y <- as.numeric(as.factor(y))-1


# Figure 1 (right)
# plot x with colors corresponding to y
matplot(tt,t(x), type="l", lty = y+1, col = -3*y + 4, xlab = "relative time",
        ylab = "hip abduction-adduction acceleration / 10000", bty = "n", lwd = 2)
legend("topright", legend = c("female","male"), col = c(4,1), lty = 1:2, lwd = 2,
       bty = "n")
