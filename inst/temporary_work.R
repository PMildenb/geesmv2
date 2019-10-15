mat.sqrt_pm <- function(A) {
  ei  <- eigen(A)
  d   <- ei$values
  d   <- (d + abs(d))/2
  if(isSymmetric(A)) ei_vec_inv <- t(ei$vectors)
  else               ei_vec_inv <- solve(ei$vectors)
  ans <- ei$vectors %*% diag(sqrt(d)) %*% ei_vec_inv
  return(ans)  
}

C <- matrix(c(2,0,0,2),nrow=2)
B <- matrix(c(20,1,-3,2),nrow=2)

mat.sqrt_pm(C)
geesmv::mat.sqrt(C)

BB <- mat.sqrt_pm(B)
geesmv::mat.sqrt(B)

A <- B

BB <- geesmv2::mat.sqrt(B)
BB <- geesmv::mat.sqrt(B)
BB %*% BB


###################################

file.copy("F://R-Code/gmds2019/modified_geesmv/GEE.var.wl_pm.R",
          "F://R-Code/geesmv2/R/GEE.var.wl_bin.aggr.R")
