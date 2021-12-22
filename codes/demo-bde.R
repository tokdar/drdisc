source("./codes/new-codes.R")

## Density estimation (no jump)
## Simulate Data
n.obs <- 2e3
ix <- sample(3, n.obs, TRUE, p0)
y.obs <- 2*rbeta(n.obs, shape1[ix], shape2[ix]) - 1
#hist(y.obs, 20, freq=FALSE, xlim=c(-1,1))

## Fit Model
fit <- bde(y.obs, "reject", burn=100, thin=9)
cat("Run time", round(fit$runtime[3]), "seconds\n")
fsamp <- with(fit, apply(b, 2, get.f))

## Visualize Fit
par(mfrow=c(1,2))
hist(y.obs, 20, freq=FALSE, xlim=c(-1,1))
grid()
for(i in 1:100) lines(y.grid, fsamp[,i], col=tcol(2, .3))
lines(y.grid, f0, col=1, lwd=2)
lines(y.grid, f0b, col=4, lwd=2, lty=2)

with(fit, boxplot(t(b), col=tcol(1+1:order,.5), border=tcol(1+1:order,.7), ylim=c(-4,4)))
grid()
points(1:order, b0, pch="+", cex=2)

