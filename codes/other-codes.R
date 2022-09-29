latex_matrix_display <- function(result, N_para, nsamp){
  all_latex <- NULL
  p <- ncol(result[[1]]$Sigma_adap[,,nsamp])
  for(i in 1:N_para){
    mat_string <- NULL
    result_curr <- result[[i]]$Sigma_adap[,,nsamp]
    
    for(j in 1:p){
      current_row <- paste0(format(result_curr[j,], digits = 2), collapse = " & ")
      mat_string <- c(mat_string, current_row)
    }
    mat_latex <- paste0(mat_string, collapse = " \\\\ ")
    mat_latex <- paste0("\\begin{bmatrix} \n", mat_latex, "\n \\end{bmatrix}")
    
    all_latex <- c(all_latex, mat_latex)
  }
  
  ret <- paste0(all_latex, collapse = ", \n ")
  cat(ret)
}


multi_chain_coda <- function(result_multi, burn = 1000, N_para, p, order, moreburn = 0){
  nlast <- length(result_multi[[1]]$b[1,1,])
  coda_result <- list()
  matrix_result <- NULL
  for(i in 1:N_para){
    matrix_result <- NULL
    for(j in 1:order){
      for(k in 1:p){
        bijk <- matrix(result_multi[[i]]$b[j,k,(moreburn+1):nlast], ncol = 1)
        colnames(bijk) <- paste0("b",j,k)
        matrix_result <- cbind(matrix_result, bijk)
      }
    }
    a_result <- t(result_multi[[i]]$a)
    a_result <- a_result[(moreburn+1):nlast,]
    colnames(a_result) <- paste0("a",1:p)
    lambda <- result_multi[[i]]$jbw
    lambda <- lambda[(moreburn+1):nlast]
    matrix_result <- cbind(matrix_result, a_result, lambda)
    coda_result[[i]] <- mcmc(matrix_result, start = burn+moreburn+1)
  }
  coda_result_list <- mcmc.list(coda_result)
  return(coda_result_list)
}


multi_chain_coda_nol <- function(result_multi, burn = 1000, N_para, p, order,
                                 moreburn = 0, nlast = NULL){
  if(is.null(nlast)) nlast <- length(result_multi[[1]]$b[1,1,])
  coda_result <- list()
  matrix_result <- NULL
  for(i in 1:N_para){
    matrix_result <- NULL
    for(j in 1:order){
      for(k in 1:p){
        bijk <- matrix(result_multi[[i]]$b[j,k,(moreburn+1):nlast], ncol = 1)
        colnames(bijk) <- paste0("b",j,k)
        matrix_result <- cbind(matrix_result, bijk)
      }
    }
    a_result <- t(result_multi[[i]]$a)
    a_result <- a_result[(moreburn+1):nlast,]
    colnames(a_result) <- paste0("a",1:p)
    matrix_result <- cbind(matrix_result, a_result)
    coda_result[[i]] <- mcmc(matrix_result, start = burn+moreburn+1)
  }
  coda_result_list <- mcmc.list(coda_result)
  return(coda_result_list)
}

multi_chain_coda_nol_onlya <- function(result_multi, burn = 1000, N_para, p, order,
                                 moreburn = 0, nlast = NULL){
  if(is.null(nlast)) nlast <- length(result_multi[[1]]$b[1,1,])
  coda_result <- list()
  matrix_result <- NULL
  for(i in 1:N_para){
    a_result <- t(result_multi[[i]]$a)
    a_result <- a_result[(moreburn+1):nlast,]
    colnames(a_result) <- paste0("a",1:p)
    matrix_result <- a_result
    coda_result[[i]] <- mcmc(matrix_result, start = burn+moreburn+1)
  }
  coda_result_list <- mcmc.list(coda_result)
  return(coda_result_list)
}

multi_chain_parallel_coda <- function(result_multi, burn = 1000, N_para, p, order, moreburn = 0){
  nlast <- length(result_multi[[1]]$b[1,1,])
  coda_result <- list()
  matrix_result <- NULL
  for(i in 1:N_para){
    matrix_result <- NULL
    for(j in 1:order){
      for(k in 1:p){
        bijk <- matrix(result_multi[[i]]$b[j,k,(moreburn+1):nlast], ncol = 1)
        colnames(bijk) <- paste0("b",j,k)
        matrix_result <- cbind(matrix_result, bijk)
      }
    }
    a_result <- t(result_multi[[i]]$a_PT[,1,])
    a_result <- a_result[(moreburn+1):nlast,]
    colnames(a_result) <- paste0("a",1:p)
    jbw_result <- result_multi[[i]]$jbw
    jbw_result <- jbw_result[(moreburn+1):nlast]
    matrix_result <- cbind(matrix_result, a_result)
    coda_result[[i]] <- mcmc(matrix_result, start = burn+moreburn+1)
  }
  coda_result_list <- mcmc.list(coda_result)
  return(coda_result_list)
}

a_parallel_plot <- function(PT_order15_jbw_L3_new, i, p, L, first = TRUE){
  par(mfrow = c(p, 1))
  if(first){
    for(j in 1:p){
      plot(PT_order15_jbw_L3_new[[i]]$a_PT[j,1,], type="l",
           ylim=range(PT_order15_jbw_L3_new[[i]]$a_PT[j,1:2,]))
      lines(PT_order15_jbw_L3_new[[i]]$a_PT[j,2,], col=tcol(2,.5), lwd=2)
    }
  }
  else{
    for(j in 1:p){
      plot(PT_order15_jbw_L3_new[[i]]$a_PT[j,1,], type="l",
           ylim=range(PT_order15_jbw_L3_new[[i]]$a_PT[j,,]))
      for(k in 2:L){
        lines(PT_order15_jbw_L3_new[[i]]$a_PT[j,k,], col=tcol(k,.5), lwd=2)
      }
    }
  }
}

get.result.bdregjump <- function(fit.x, x.obs, y.obs, b0x, a0,
                                 trim = 1, thin_index=NULL){
  y.grid <- seq(-trim, trim, .01)
  cat("Runtime", round(fit.x$runtime[3]), "seconds\n")
  if(is.null(thin_index)){
    b.est <- fit.x$b
    a.est <- fit.x$a
    shapes.est <- fit.x$shapes
    jbw.est <- fit.x$jbw
  }else{
    b.est <- fit.x$b[,,thin_index]
    a.est <- fit.x$a[,thin_index]
    shapes.est <- fit.x$shapes[,thin_index]
    jbw.est <- fit.x$jbw[thin_index]
  }
  
  p <- dim(fit.x$b)[2]
  
  # plot1
  par(mfrow = c(p,2))
  
  for(j in 1:p){
    boxplot(data.frame(b=t(b.est[,j,]),a=a.est[j,]), col=tcol(1+1:(1+order),.3), border=tcol(1+1:(1+order),.7),
            main = paste0("predictor", j))
    grid()
    points(c(b0x[,j], a0[j]), pch="+", cex=2)
  }
  
  # plot2
  
  plot(a.est[1,],ty="l", ylim=range(a.est), ylab="a", xlab="MCMC draw", main = TeX("traceplot of $\\alpha$"))
  for(i in 1:length(a0)) lines(a.est[i,], col=tcol(i,.5), lwd=2)
  abline(h = a0, lty=2, col=1:length(a0))
  
  plot(jbw.est, ty="l", xlab="MCMC draw", ylab="jump bw", col=tcol(1,.5), lwd=2, main = TeX("traceplot of $\\lambda$"))
  abline(h = jbw0, lty=2)
  
  cat("Estimated jump persistence (true=", pers0, ")\n", round(quantile(jbw.est, c(0.025, .25,.5,.75, 0.975))/0.32,2), "\n")
  
  # plot3
  par(mfrow = c(2,2))
  xnew <- rep(1,4)
  for(i in 1:(p-1)) xnew <- cbind(xnew, seq(-2,2,1.25))
  nnew <- nrow(xnew)
  xb0 <- tcrossprod(xnew, b0x)
  xa0 <- c(xnew %*% a0)
  if(pos.synth) xa0 <- pmax(0, xa0)
  f0x <- sapply(1:nnew, function(i) get.f(b=xb0[i,],a=xa0[i], sh=shapes0, jbw=jbw0,
                                          trim=trim))
  
  nsamp <- ncol(a.est)
  for(cl in 1:4){
    ix.cl <- abs(x.obs[,2] - xnew[cl,2]) < 0.25
    hist(y.obs[ix.cl], freq=FALSE, main="", xlab="Y", xlim=c(-trim,trim))
    ylim <- par("usr")[3:4]
    xb <- apply(b.est, 3, function(b) tcrossprod(xnew[cl,,drop=FALSE], b))
    xa <- c(xnew[cl,,drop=FALSE] %*% a.est)
    if(pos.est) xa <- pmax(xa, 0)
    fx <- sapply(1:nsamp, function(s) get.f(b=xb[,s],a=xa[s],sh=shapes.est[,s],jbw=jbw.est[s],
                                            trim=trim))
    for(s in 1:nsamp) lines(y.grid, fx[,s], col=tcol(2,.5))
    lines(y.grid, rowMeans(fx), col=4, lwd=2)
    lines(y.grid, f0x[,cl], lwd=2)
    drop_pos <- trim/0.01
    dropx <- 1-fx[drop_pos,]/fx[drop_pos+1,]
    drop.est <- round(100*median(dropx))
    drop.ci <- round(100*quantile(dropx, pr=c(.025,.975)))
    text(-1, ylim[2]*0.9, bquote('Drop%'==.(drop.est)[paste("[",.(drop.ci[1]),",",.(drop.ci[2]),"]")]), pos=4)
  }
  
  cat("accperatnce rate:", round(100*fit.x$acpt), "%\n")
}

get.result.bdregjump.new <- function(fit.x, x.obs, y.obs, b0x, a0,
                                 trim = 1, thin_index=NULL){
  ############################################################  
  ###### accommodate true order not equal to fitted order#####
  ###### add true parameters = 0 for redundant dimensions#####
  ############################################################
  
  y.grid <- seq(-trim, trim, .01)
  cat("Runtime", round(fit.x$runtime[3]), "seconds\n")
  if(is.null(thin_index)){
    b.est <- fit.x$b
    a.est <- fit.x$a
    shapes.est <- fit.x$shapes
    jbw.est <- fit.x$jbw
  }else{
    b.est <- fit.x$b[,,thin_index]
    a.est <- fit.x$a[,thin_index]
    shapes.est <- fit.x$shapes[,thin_index]
    jbw.est <- fit.x$jbw[thin_index]
  }
  
  p <- dim(fit.x$b)[2]
  order <- dim(b.est)[1]
  # plot1
  par(mfrow = c(p,2))
  
  for(j in 1:p){
    boxplot(data.frame(b=t(b.est[,j,]),a=a.est[j,]), col=tcol(1+1:(1+order),.3), border=tcol(1+1:(1+order),.7),
            main = paste0("predictor", j))
    grid()
    points(c(b0x[,j], rep(0, order - nrow(b0x)) ,a0[j]), pch="+", cex=2)
  }
  
  # plot2
  
  plot(a.est[1,],ty="l", ylim=range(a.est), ylab="a", xlab="MCMC draw", main = TeX("traceplot of $\\alpha$"))
  for(i in 1:length(a0)) lines(a.est[i,], col=tcol(i,.5), lwd=2)
  abline(h = a0, lty=2, col=1:length(a0))
  
  plot(jbw.est, ty="l", xlab="MCMC draw", ylab="jump bw", col=tcol(1,.5), lwd=2, main = TeX("traceplot of $\\lambda$"))
  abline(h = jbw0, lty=2)
  
  cat("Estimated jump persistence (true=", pers0, ")\n", round(quantile(jbw.est, c(0.025, .25,.5,.75, 0.975))/0.32,2), "\n")
  
  # plot3
  par(mfrow = c(2,2))
  xnew <- rep(1,4)
  for(i in 1:(p-1)) xnew <- cbind(xnew, seq(-2,2,1.25))
  nnew <- nrow(xnew)
  xb0 <- tcrossprod(xnew, b0x)
  xa0 <- c(xnew %*% a0)
  if(pos.synth) xa0 <- pmax(0, xa0)
  lpoly <- list(
    #  Vectorize(function(x) return(1/sqrt(2)), "x"), ## k = 0
    function(x) return(sqrt(3/2)*x), ## k = 1
    function(x) return(sqrt(5/8)*(3*x^2 - 1)), ## k = 2
    function(x) return(sqrt(7/8)*(5*x^3 - 3*x)) ## k = 3
  )
  lpoly <- lpoly[1:ncol(b0x)]
  get.poly.mat <- function(y) return(sapply(lpoly, function(f) f(y)))
  yFn <- get.poly.mat
  wFn <- function(y, beta, alpha=0, bw=0.16) return(c(matrix(yFn(y), nrow=length(y)) %*% beta) + alpha*half.kern(y,bw))
  fFn <- function(y, shapes=c(1,1), trim=1, ...) return((1/2)*dtrunc((y+1)/2, "beta", a = (1-trim)/2, b = (1+trim)/2,
                                                                     shape1 = shapes[1], shape2 = shapes[2])*Phi(wFn(y, ...)))
  get.f <- function(b, a=0, sh=c(1,1), jbw=0.16, trim=1) {
    y.grid <- seq(-trim, trim, .01)
    f.grid <- fFn(y.grid, trim=trim, beta=b, alpha=a, shapes=sh, bw=jbw)
    normalization.const <- (integrate(fFn, -trim, 0, beta=b, alpha=a, shapes=sh, bw=jbw, trim = trim)$value
                            + integrate(fFn, 0, trim, beta=b, alpha=a, shapes=sh, bw=jbw, trim = trim)$value)

    return(f.grid/normalization.const)
  }
  f0x <- sapply(1:nnew, function(i) get.f(b=xb0[i,],a=xa0[i], sh=shapes0, jbw=jbw0,
                                          trim=trim))
  
  nsamp <- ncol(a.est)
  for(cl in 1:4){
    ix.cl <- abs(x.obs[,2] - xnew[cl,2]) < 0.25
    hist(y.obs[ix.cl], freq=FALSE, main="", xlab="Y", xlim=c(-trim,trim))
    ylim <- par("usr")[3:4]
    xb <- apply(b.est, 3, function(b) tcrossprod(xnew[cl,,drop=FALSE], b))
    xa <- c(xnew[cl,,drop=FALSE] %*% a.est)
    if(pos.est) xa <- pmax(xa, 0)
    lpoly <- list(
      #  Vectorize(function(x) return(1/sqrt(2)), "x"), ## k = 0
      function(x) return(sqrt(3/2)*x), ## k = 1
      function(x) return(sqrt(5/8)*(3*x^2 - 1)), ## k = 2
      function(x) return(sqrt(7/8)*(5*x^3 - 3*x)) ## k = 3
    )
    lpoly <- lpoly[1:order]
    fx <- sapply(1:nsamp, function(s) get.f(b=xb[,s],a=xa[s],sh=shapes.est[,s],jbw=jbw.est[s],
                                            trim=trim))
    for(s in 1:nsamp) lines(y.grid, fx[,s], col=tcol(2,.5))
    lines(y.grid, rowMeans(fx), col=4, lwd=2)
    lines(y.grid, f0x[,cl], lwd=2)
    drop_pos <- trim/0.01
    dropx <- 1-fx[drop_pos,]/fx[drop_pos+1,]
    drop.est <- round(100*median(dropx))
    drop.ci <- round(100*quantile(dropx, pr=c(.025,.975)))
    text(-1, ylim[2]*0.9, bquote('Drop%'==.(drop.est)[paste("[",.(drop.ci[1]),",",.(drop.ci[2]),"]")]), pos=4)
  }
  
  cat("accperatnce rate:", round(100*fit.x$acpt), "%\n")
}
                 
                 
get.result.data <- function(fit.x, x.obs, y.obs, b0x, a0, thin_index=NULL){
  x <- x.obs
  cat("Runtime", round(fit.x$runtime[3]), "seconds\n")
  if(is.null(thin_index)){
    b.est <- fit.x$b
    a.est <- fit.x$a
    shapes.est <- fit.x$shapes
    jbw.est <- fit.x$jbw
  }else{
    b.est <- fit.x$b[,,thin_index]
    a.est <- fit.x$a[,thin_index]
    shapes.est <- fit.x$shapes[,thin_index]
    jbw.est <- fit.x$jbw[thin_index]
  }
  
  p <- dim(fit.x$b)[2]
  
  # plot1
  par(mfrow = c(p,1))
  
  for(j in 1:p){
    boxplot(data.frame(b=t(b.est[,j,]),a=a.est[j,]), col=tcol(1+1:(1+order),.3), border=tcol(1+1:(1+order),.7))
    grid()
    points(c(b0x[,j], a0[j]), pch="+", cex=2)
  }
  
  cat("accperatnce rate:", round(100*fit.x$acpt), "%\n")
  
  
  # plot3
  par(mfrow = c(2,1))
  plot(a.est[1,],ty="l", ylim=range(a.est), ylab="a", xlab="MCMC draw")
  for(i in 1:length(a0)) lines(a.est[i,], col=tcol(i,.5), lwd=2)
  abline(h = a0, lty=2, col=1:length(a0))
  
  plot(jbw.est, ty="l", xlab="MCMC draw", ylab="jump bw", col=tcol(1,.5), lwd=2)
  abline(h = jbw0, lty=2)
  
  cat("Estimated jump persistence (true=", pers0, ")\n", round(quantile(jbw.est, c(0.025, .25,.5,.75, 0.975))/0.32,2), "\n")
}

plot_multi_a <- function(result_M1_chol, a0){
  for(i in 1:length(a0)){
    plot(result_M1_chol[[1]]$a[i,],ty="l", ylim = c(a0[i]-4, a0[i]+10),
         ylab=paste0("a",(i-1)), xlab="MCMC draw")
    for(j in 2:N_para){
      lines(result_M1_chol[[j]]$a[i,], col=rainbow(j)[j], lwd=2)
    }
    abline(h = a0[i], lty=2)
  }
  acpt <- NULL
  for(i in 1:N_para){
    acpt <- rbind(acpt, result_M1_chol[[i]]$acpt)
  }
  return(acpt)
}


real_data_plot <- function(fit.x, x, x.test, y.test, yFn){
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
}

get.swap.ratio <- function(adapt_PT_new_rho, N_parallel = 6, L = 5, moreburn = 0){
  nlast <- ncol(adapt_PT_new_rho[[1]]$swap_info)
  result <- matrix(NA, N_para, L-1)
  for(i in 1:N_para){
    current_info <- adapt_PT_new_rho[[i]]$swap_info[,(1+moreburn):nlast]
    for(j in 1:(L-1)){
      swap_ind <- (current_info[1,] == j)
      result[i,j] <- mean(current_info[3,swap_ind])
    }
  }
  return(result)
}
