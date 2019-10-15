mat.sqrt <-
function(A)
{
   ei<-eigen(A)
   d<-ei$values
   d<-(d+abs(d))/2
   d2<-sqrt(d)
   if(isSymmetric(A)) ei_vec_inv <- t(ei$vectors)
   else               ei_vec_inv <- solve(ei$vectors)
   ans <- ei$vectors %*% diag(sqrt(d)) %*% ei_vec_inv
   return(ans)
}
