source("./codes/new-codes.R")

## Density regression
n.obs <- 2000
b0x <- cbind(b0, threshold(rt(order, df=6),1))
x.obs <- cbind(1, rnorm(n.obs))
y.obs <- simulate.fx(n.obs, b0x, x.obs, a=rep(0,ncol(x.obs)))

fit.x <- bdreg(y=y.obs, x=x.obs, burn=100, thin=9)
cat("Runtime", round(fit.x$runtime[3]), "seconds\n")
b.est <- fit.x$b

par(mfrow = c(3,2))

for(j in 1:2){
  boxplot(t(b.est[,j,]), col = tcol(1+1:order,.3), border=tcol(1+1:order,.7))
  grid()
  points(b0x[,j], pch="+", cex=2)
}

xnew <- cbind(1, seq(-2,2,1.25))
xb0 <- tcrossprod(xnew, b0x)
f0x <- apply(xb0, 1, get.f)

for(cl in 1:4){
  ix.cl <- abs(x.obs[,2] - xnew[cl,2]) < 0.5
  hist(y.obs[ix.cl], freq=FALSE)
  xb <- apply(b.est, 3, function(b) tcrossprod(xnew[cl,,drop=FALSE], b))
  fx <- apply(xb, 2, get.f)
  for(i in 1:100) lines(y.grid, fx[,i], col=tcol(2,.5))
  lines(y.grid, f0x[,cl], lwd=2)
}



