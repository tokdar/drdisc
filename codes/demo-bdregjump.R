adapt <- FALSE
bunch <- FALSE
source("./codes/new-codes.R")
p <- 2
  
simulate <- TRUE
## Density regression with jump
if(simulate){
  n.obs <- 2e3
  b0x <- cbind(b0, threshold(rt(order, df=6),1))
  while(ncol(b0x) < p) b0x <- cbind(b0x, 0)
  a0 <- c(1,0.5,rep(0,p-2))
  shapes0 <- c(4,1)
  pers0 <- 1/2
  jbw0 <- 0.16*2*pers0
  pos.synth <- TRUE
  pos.est <- TRUE
  
  x.obs <- cbind(1, matrix(rnorm(n.obs*(p-1)), n.obs, p-1))
  y.obs <- simulate.fx(n.obs, b0x, x.obs, a0, shapes0, jbw0, positive=pos.synth)
}

cat("a0:", round(a0,2), "\nb0:\n")
print(signif(b0x,3))

fit.x <- bdregjump(y=y.obs, x=x.obs, 
                   burn=100, nsamp=100, thin=9, 
                   jump=list(positive=pos.est, update.jbw=TRUE))
cat("Runtime", round(fit.x$runtime[3]), "seconds\n")
b.est <- fit.x$b
a.est <- fit.x$a
shapes.est <- fit.x$shapes
jbw.est <- fit.x$jbw

par(mfrow = c(4,2), mar=c(4,3,2,2)+.1, mgp=c(2,.5,0))

for(j in 1:2){
  boxplot(data.frame(b=t(b.est[,j,]),a=a.est[j,]), col=tcol(1+1:(1+order),.3), border=tcol(1+1:(1+order),.7))
  grid()
  points(c(b0x[,j], a0[j]), pch="+", cex=2)
}

xnew <- cbind(1, seq(-2,2,1.25))
nnew <- nrow(xnew)
xb0 <- tcrossprod(xnew, b0x)
xa0 <- c(xnew %*% a0)
if(pos.synth) xa0 <- pmax(0, xa0)
f0x <- sapply(1:nnew, function(i) get.f(b=xb0[i,],a=xa0[i], sh=shapes0, jbw=jbw0))

nsamp <- ncol(a.est)
for(cl in 1:4){
  ix.cl <- abs(x.obs[,2] - xnew[cl,2]) < 0.25
  hist(y.obs[ix.cl], freq=FALSE, main="", xlab="Y", xlim=c(-1,1))
  ylim <- par("usr")[3:4]
  xb <- apply(b.est, 3, function(b) tcrossprod(xnew[cl,,drop=FALSE], b))
  xa <- c(xnew[cl,,drop=FALSE] %*% a.est)
  if(pos.est) xa <- pmax(xa, 0)
  fx <- sapply(1:nsamp, function(s) get.f(b=xb[,s],a=xa[s],sh=shapes.est[,s],jbw=jbw.est[s]))
  for(s in 1:nsamp) lines(y.grid, fx[,s], col=tcol(2,.5))
  lines(y.grid, rowMeans(fx), col=4, lwd=2)
  lines(y.grid, f0x[,cl], lwd=2)
  dropx <- 1-fx[100,]/fx[101,]
  drop.est <- round(100*median(dropx))
  drop.ci <- round(100*quantile(dropx, pr=c(.025,.975)))
  text(-1, ylim[2]*0.9, bquote('Drop%'==.(drop.est)[paste("[",.(drop.ci[1]),",",.(drop.ci[2]),"]")]), pos=4)
}

cat("accperatnce rate:", round(100*fit.x$acpt), "%\n")
plot(a.est[1,],ty="l", ylim=range(a.est), ylab="a", xlab="MCMC draw")
for(i in 1:length(a0)) lines(a.est[i,], col=tcol(i,.5), lwd=2)
abline(h = a0, lty=2, col=1:length(a0))

plot(jbw.est, ty="l", xlab="MCMC draw", ylab="jump bw", col=tcol(1,.5), lwd=2)
abline(h = jbw0, lty=2)

cat("Estimated jump persistence (true=", pers0, ")\n", round(quantile(jbw.est, c(0.025, .25,.5,.75, 0.975))/0.32,2), "\n")

