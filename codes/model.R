vote <- read.csv("./codes/vote_data.csv")
str(vote)
vote <- subset(vote, abs(from_requirement_threshold) < 0.5)
set.seed(1234)
train.ss <- sample(nrow(vote), 15e3)

preds <- model.matrix(~ #major.requirement 
                        + ISS.against.recommendation 
                      #+ shares.oustanding.base 
                      #+ special.meeting 
                      + analyst.coverage 
                      #+ institutional.ownership 
                      + past.stock.return 
                      + q 
                      + firm.size, 
                      data = vote)[,-1]
preds.train <- scale(preds[train.ss,])
x <- cbind(1, preds.train)
x.test <- cbind(1, scale(preds[-train.ss,],
                         center=attr(preds.train,"scaled:center"),
                         scale=attr(preds.train,"scaled:scale")))
dimnames(x)[[2]][1] <- "Intercept"
dimnames(x.test)[[2]][1] <- "Intercept"
x.names <- dimnames(x)[[2]]

y <- 2*vote[train.ss,"from_requirement_threshold"]
y.test <- 2*vote[-train.ss,"from_requirement_threshold"]

adapt <- FALSE
bunch <- FALSE
pos.est <- TRUE
source("./codes/new-codes.R")
fit.x <- bdregjump(y=y, x=x, burn=100, nsamp=100, thin=9, 
                   jump=list(positive=pos.est, update.jbw=TRUE))
save(train.ss, fit.x, file=paste0("run-", Sys.time(), ".Rd"))
cat("Runtime", round(fit.x$runtime[3]), "seconds\n")
b.est <- fit.x$b; dimnames(b.est)[[2]] <- x.names
a.est <- fit.x$a; dimnames(a.est)[[1]] <- x.names
shapes.est <- fit.x$shapes
jbw.est <- fit.x$jbw

a.dt <- data.frame(t(round(apply(a.est, 1, quantile, pr = c('Est'=0.5, 'CI95%Lo'=0.025, 'CI95%Hi'=0.975),names=FALSE), 2)))
names(a.dt) <- c("Estimate", "CI95%Lo", "CI95%Hi")
par(mfrow = c(4,2), mar=c(4,3,2,2)+.1, mgp=c(2,.5,0))
cat("Estimates of alpha\n")
print(a.dt)

n <- nrow(x)
p <- ncol(x)
for(j in 1:p){
  boxplot(data.frame(b=t(b.est[,j,]),a=a.est[j,]), col=tcol(1+1:(1+order),.3), border=tcol(1+1:(1+order),.7))
  grid()
  abline(h=0,lty=2)
  title(main=x.names[j], xlab="Parameter", ylab="Posterior sample")
  #points(c(b0x[,j], a0[j]), pch="+", cex=2)
}

cat("accperatnce rate:", round(100*fit.x$acpt), "%\n")
plot(a.est[1,],ty="n", ylim=range(a.est), main="Traceplot: a", xlab="MCMC iteration", ylab="Posterior draw")
for(i in 1:p) lines(a.est[i,], col=tcol(i,.5), lwd=2)
#abline(h = a0, lty=2, col=1:length(a0))

plot(jbw.est, ty="l", xlab="MCMC draw", ylab="Posterior draw", main="Traceplot: bw", col=tcol(1,.5), lwd=2)
#abline(h = jbw0, lty=2)

cat("Estimated jump persistence\n", round(quantile(jbw.est, c(0.025, .25,.5,.75, 0.975))/0.32,2), "\n")

x.uniq.rows <- which(!duplicated(x))
n.unique <- length(x.uniq.rows)
rows.to.assess <- sort(sample(x.uniq.rows, min(8, n.unique)))
xnew <- x[rows.to.assess,,drop=FALSE]
nnew <- nrow(xnew)
prednew <- apply(round(preds[train.ss[rows.to.assess],],1), 1, paste, collapse=", ")
nsamp <- ncol(a.est)
xanew <- xnew %*% a.est
if(pos.est) xanew <- matrix(pmax(0, xanew), nrow=nnew)

for(cl in 1:nnew){
  ix.cl <- apply(x.test, 1, function(z) all(abs(z - xnew[cl,]) < 1))
  hist(y.test[ix.cl], 20, freq=FALSE, main=paste0("X=(",prednew[cl],")"), xlab="Holdout Y", xlim=c(-1,1))
  ylim <- par("usr")[3:4]
  xb <- apply(b.est, 3, function(b) tcrossprod(xnew[cl,,drop=FALSE], b))
  xa <- c(xnew[cl,,drop=FALSE] %*% a.est)
  if(pos.est) xa <- pmax(xa, 0)
  fx <- sapply(1:nsamp, function(s) get.f(b=xb[,s],a=xa[s], sh=shapes.est[,s],jbw=jbw.est[s]))
  for(s in 1:nsamp) lines(y.grid, fx[,s], col=tcol(2,.5))
  lines(y.grid, rowMeans(fx), col=4, lwd=2)
  dropx <- 1-fx[100,]/fx[101,]
  drop.est <- round(100*median(dropx))
  drop.ci <- round(100*quantile(dropx, pr=c(.025,.975)))
  text(-1, ylim[2]*0.9, bquote('Drop%'==.(drop.est)[paste("[",.(drop.ci[1]),",",.(drop.ci[2]),"]")]), pos=4)
}


stop("stop here")
