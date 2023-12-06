# data import
dta <- readRDS("data/data_comb.RDS")
names(dta)

# choose functional variable (compare functional linear model)
x <- dta$knee_accl_vt[dta$cond=="slowbw",]
x <- x/10000
tt <- (dta$cycle-1)/100


# calculate euclidian distances
help("dist")
dst <- dist(x)

# hierarchical clustering
help(hclust)

hs <- hclust(dst, method = "single")
ha <- hclust(dst, method = "average")
hc <- hclust(dst, method = "complete")
hw <- hclust(dst, method = "ward.D2")

# dendrograms
plot(hs)
plot(ha)
plot(hc)
plot(hw)

# use complete linkage

# two cluster solutions
cc <- cutree(hc, k=2)

# compare clusters with sex
ys <- dta$sex[dta$cond=="slowbw"]
tyc <- table(ys,cc)
barplot(tyc, col = c(7,4), xlab = "cluster", ylab = "frequency")
legend("topright", legend = c("female","male"), fill = c(7,4))

# age per cluster
ya <- dta$age[dta$cond=="slowbw"]
boxplot(ya ~ ys*cc, varwidth = TRUE, col = c(7,4))

# show results as dot plot
cx <- rep(4,length(ys))
cx[ys=="f"] <- 7

set.seed(1234)
plot(cc + rnorm(31, sd=0.15), ya, col = cx, pch = cx+1, lwd = 5, xaxp = c(1,2,1),
     xlab = "cluster", ylab = "age")
legend("topright", legend = c("female","male"), col = c(7,4), pch = c(7,4)+1, pt.lwd = 5)


# FPCA
library(refund)
help("fpca.face")
fpcax <- fpca.face(Y = x, pve = 0.95, p=3, m=1)

# make eigenfunctions orthonormal w.r.t. the L2 and not vector inner product
phi <- fpcax$efunctions
L <- length(tt) 
phi <- phi * sqrt(L) 

# rescale scores according to rescaling of eigenfunctions
xi <- fpcax$scores / sqrt(L) 

# k-means using scores

# within sum of squares
K <- 20
wss <- numeric(K)
for (k in 1:K) {
  wss[k] <- kmeans(xi, centers = k)$tot.withinss
}
plot(wss, type = "b")

# use two cluster solution
kmx <- kmeans(xi, centers = 2)

# and compare to complete linkage
table(kmx$cluster, cc)



# model-based clustering
# using funHDDC
library(funHDDC)
help(funHDDC)


# basis representation
basis<- create.bspline.basis(c(0,1), nbasis=25)
var1<-smooth.basis(argvals=seq(0,1,length.out = 101),y=t(x),fdParobj=basis)$fd

# clustering with 2 clusters 
res.uni<-funHDDC(var1,K=2,
                 model=c('AkjBkQkDk','AkjBQkDk','AkBkQkDk','ABkQkDk','AkBQkDk','ABQkDk'),
                 init="kmeans",threshold=0.2)

# compare with complete linkage
names(res.uni)
table(res.uni$class,cc)



# plot of pc scores 
cx <- rep(3,length(cc))
cx[cc==2] <- 1 
plot(xi[,1:2], type = "n",
     xlab = "component 1", ylab = "component 2",
     xlim = c(-0.7,0.7), ylim = c(-0.6,0.6),  bty = "n")
text(xi[,1], xi[,2], 1:length(cx), col = cx)
text(xi[xi[,2]>0.4,1],xi[xi[,2]>0.4,2],
     "12", col = 4,)
legend("bottomright", legend = c("cluster 1","cluster 2"),
       fill = c(3,1), bty = "n")


# clustering/plot of functions
matplot(tt,t(x), type="l", lty = cc, col = cx, xlab = "relative time",
        ylab = "knee axial rotation acceleration / 10000", bty = "n", lwd = 2,
        ylim = c(-3,2))
lines(tt,x[xi[,2]>0.4,], col = 4, lwd = 2)
lines(tt, phi[,1], lty = 1, lwd = 2, col = 6)
legend("bottomright", legend = c("1st. eigenfunction","cluster 1","cluster 2"),
       col = c(6,3,1), lty = c(1,1,2), lwd = 2,
       bty = "n")
