mat.sqrt.inv <-
function(A)
{
   ei<-eigen(A)
   d<-ei$values
   d<-(d+abs(d))/2
   d2<-1 / sqrt(d)
   d2[d == 0]<-0
   if(isSymmetric(A)) ei_vec_inv <- t(ei$vectors)
   else               ei_vec_inv <- solve(ei$vectors)
   ans <- ei$vectors %*% diag(d2) %*% ei_vec_inv
   return(ans)
}
