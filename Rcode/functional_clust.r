# data import
dta <- readRDS("data/data_comb.RDS")
names(dta)

# choose functional variable (compare functional linear model)
x <- dta$knee_accl_vt[dta$cond=="slowbw",]
x <- x/10000
tt <- (dta$cycle-1)/100
sx <- dta$sex[dta$cond=="slowbw"]

matplot(tt,t(x), type="l", lty = 1, col = (sx == "f") + 1,
        xlab = "percentage cycle points",
        ylab = "knee acceleration vertical/10000", bty = "n", lwd = 2)




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

# two cluster solutions
cc <- cutree(hc, k=2)
cw <- cutree(hw, k=2)
table(cc,cw)

# use complete linkage
pdf("hclust1.pdf", width = 5, height = 6)
plot(hc)
dev.off()

pdf("hclust1_d.pdf", width = 4, height = 5)
plot(hc, main = "", ylab = "Höhe", xlab = "Beobachtungsnummer", sub = "", cex = 0.7)
dev.off()



# compare clusters with sex
ys <- dta$sex[dta$cond=="slowbw"]
tyc <- table(ys,cc)
tyw <- table(ys,cw)
barplot(tyc, col = c(7,4), xlab = "cluster", ylab = "frequency")
legend("topright", legend = c("female","male"), fill = c(7,4))
chisq.test(tyc)

# age per cluster
ya <- dta$age[dta$cond=="slowbw"]
t.test(ya ~ cc)
boxplot(ya ~ ys*cc, varwidth = TRUE, col = c(7,4))

# show results as dot plot
cx <- rep(4,length(ys))
cx[ys=="f"] <- 7

pdf("hclust2.pdf", width = 5, height = 6)
set.seed(1234)
plot(cc + rnorm(31, sd=0.15), ya, col = cx, pch = cx+1, lwd = 5, xaxp = c(1,2,1),
     xlab = "cluster", ylab = "age")
legend("topright", legend = c("female","male"), col = c(7,4), pch = c(7,4)+1, pt.lwd = 5)
dev.off()


cx <- rep(1,length(ys))
cx[ys=="f"] <- 3

pdf("hclust2_d.pdf", width = 4, height = 5)
set.seed(54321)
plot(cc + rnorm(31, sd=0.1), ya, col = cx, pch = cx, lwd = 5, xaxp = c(1,2,1),
     xlab = "Cluster", ylab = "Alter", bty = "n")
legend("top", legend = c("Frauen","Männer"), col = c(3,1), pch = c(3,1), pt.lwd = 5)
dev.off()




# FPCA
library(refund)
help("fpca.face")
fpcax <- fpca.face(Y = x, pve = 0.95)
fpcax$scores
plot(fpcax$efunctions[,1], type = "l")

matplot(tt,t(x), type="l", lty = 1, col = (sx == "m") + 1,
        xlab = "percentage cycle points",
        ylab = "knee acceleration vertical/10000", bty = "n", lwd = 2)
lines(tt,fpcax$efunctions[,1], lwd=5, col = 4)

plot(fpcax$scores[,1:2], col = (sx == "m") + 1, pch = (sx == "m") + 1,
     lwd = 3)


# k-means using scores

# within sum of squares
K <- 20
wss <- numeric(K)
for (k in 1:K) {
  wss[k] <- kmeans(fpcax$scores, centers = k)$tot.withinss
}
plot(wss, type = "b")

# use two cluster solution
kmx <- kmeans(fpcax$scores, centers = 2)

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

# scores with col corresponding to funHDDC and style to complete linkage 
cx <- rep(3,length(cc))
cx[cc==2] <- 1 
plot(fpcax$scores[,1:2], col = res.uni$class, pch = cx, lwd = 5,
     xlab = "component 1", ylab = "component 2",
     xlim = c(-7,7), ylim = c(-7,7),  bty = "n")

# the point switching between clusters
ya[fpcax$scores[,2]>4]
ys[fpcax$scores[,2]>4]
which(fpcax$scores[,2]>4)



# plot for the paper 
pdf("kmeans1a.pdf", width = 5, height = 6)
#plot(fpcax$scores[,1:2], col = cx, pch = cx, lwd = 5,
plot(fpcax$scores[,1:2], type = "n",
     xlab = "component 1", ylab = "component 2",
     xlim = c(-7,7), ylim = c(-7,7),  bty = "n")
text(fpcax$scores[,1], fpcax$scores[,2], 1:length(cx), col = cx)
#points(fpcax$scores[fpcax$scores[,2]>4,1],fpcax$scores[fpcax$scores[,2]>4,2],
#       col = 4, lwd = 5, pch = 3)
text(fpcax$scores[fpcax$scores[,2]>4,1],fpcax$scores[fpcax$scores[,2]>4,2],
     "12", col = 4,)
legend("bottomright", legend = c("cluster 1","cluster 2"),
       fill = c(3,1),
       #pch = c(3,1), pt.lwd = 5,
       bty = "n")
dev.off()


# clustering/plot of functions
pdf("kmeans2.pdf", width = 5, height = 6)
matplot(tt,t(x), type="l", lty = cc, col = cx, xlab = "percentage cycle points",
        ylab = "knee acceleration vertical/10000", bty = "n", lwd = 2)
lines(tt,x[fpcax$scores[,2]>4,], col = 4, lwd = 2)
lines(tt, fpcax$efunctions[,1], lty = 1, lwd = 4, col = 6)
legend("bottomright", legend = c("1st. eigenfunction","cluster 1","cluster 2"),
       col = c(6,3,1), lty = c(1,1,2), lwd = 2,
       bty = "n")
dev.off()
