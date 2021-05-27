## modified version 
## - output includes gee and se for kc and md method
## - specifically for binomial setting with complete data
## - needs a specific auxiliary function "mat.sqrt_pm" defined below

require(matrixcalc)
# require(expm)
GEE.var.kc.md.fg.pan.wl_pm <-
  function(formula,id,family="binomial",data,corstr="independence",b=3/4,verbose=FALSE){
    #########################################################################
    # Arguments:
    # formula  specify the model of interest
    # family   "gaussian", "binomial" or "poisson"
    # data     data frame
    # corstr   Working correlation structure: "independence", "AR-M", "exchangeable", "unstructured".
    # value:   GEE returns a list containing the following elements
    #          cov.beta     estimate of robust variance for \hat{\beta}
    #          cov.var      estimate of the variance-covariance matrix for robust variance.
    #########################################################################

    # Delete the records with missing data in predictors or outcomes;
    if (is.null(data$id)){
      index <- which(names(data)==id)
      data$id <- data[,index]}
    
    ### Get the design matrix;
    m   <- model.frame(formula, data)
    mt  <- attr(m, "terms") 
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
    # }
    
    ### Fit the GEE model to get the estimate of parameters \hat{\beta};
    gee.fit  <- suppressMessages(gee(formula,data=data,id=id,family=family,corstr=corstr))
    beta_est <- gee.fit$coefficient
    alpha    <- gee.fit$working.correlation[1,2]
    len      <- length(beta_est)
    len_vec  <- len^2
    
    ### Estimate the robust variance for \hat{\beta}
    data$id  <- gee.fit$id
    cluster  <- cluster.size(data$id)
    ncluster <- max(cluster$n)
    size     <- cluster$m
    mat$subj <- rep(unique(data$id), cluster$n)
    
    if(is.character(corstr)){
      var <- switch(corstr,
                    "independence"=cormax.ind(ncluster),
                    "exchangeable"=cormax.exch(ncluster, alpha),
                    "AR-M"=cormax.ar1(ncluster, alpha),
                    "unstructured"=summary(gee.fit)$working.correlation)
    }else{
      print(corstr)
      stop("'working correlation structure' not recognized")
    }  

    ## gaussian/poisson REMOVED
    
    family <- "binomial"  ## added
    
    ### for-loop added for step01 in WL
    step01.wl <- matrix(0, nrow=len, ncol=len)
    for (i in 1:size){
      y<-as.matrix(data$response[data$id==unique(data$id)[i],1])  # ,1 added
      N<-y+as.matrix(data$response[data$id==unique(data$id)[i],2])  # added
      covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
      var_i=var[1:cluster$n[i],1:cluster$n[i]]
      
      ## gaussian/poisson REMOVED
      
      D<-mat.prod(covariate, exp(covariate%*%beta_est)/((1+exp(covariate%*%beta_est))^2))
      D<-diag(c(N))%*%D #added
      Vi <- diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),cluster$n[i])%*%var_i%*%diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),cluster$n[i])
      Vi<-diag(sqrt(c(N)))%*%Vi%*%diag(sqrt(c(N))) #added
      xx<-t(D)%*%solve(Vi)%*%D
      step01.wl<-step01.wl+xx 
    }
    
    cov.beta <- matrix(0, nrow=len, ncol=len)
    step.pan <- step.wl <- matrix(0, nrow = cluster$n[1], ncol = cluster$n[1])             ## added for PAN, binomial
    step11   <- matrix(0, nrow=len, ncol=len)

    for (i in 1:size){
      y<-as.matrix(data$response[data$id==unique(data$id)[i],1])  # ,1 added
      N<-y+as.matrix(data$response[data$id==unique(data$id)[i],2])  # added
      covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
      var_i=var[1:cluster$n[i],1:cluster$n[i]]

      ## added for middle part in PAN:
      residual_tmp <- (y - N*exp(covariate %*% beta_est)/(1 + exp(covariate %*% beta_est))) ## N* added
      resid <- residual_tmp %*% t(residual_tmp)
      B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
      diag(B) <- 1/sqrt(exp(covariate %*% beta_est)/(1 + exp(covariate %*% beta_est))^2)
      step.pan <- step.pan + B %*% resid %*% B

      ## gaussian/poisson case REMOVED 
      
      D  <- mat.prod(covariate, exp(covariate %*% beta_est)/((1+exp(covariate %*% beta_est))^2))
      D  <- diag(c(N)) %*% D #added
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1+exp(covariate %*% beta_est))^2)),cluster$n[i]) %*%
            var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1+exp(covariate %*% beta_est))^2)),cluster$n[i])
      Vi <- diag(sqrt(c(N)))%*%Vi%*%diag(sqrt(c(N))) #added

      ## added for middle part in WL:
      resid<-B%*%solve(cormax.ind(cluster$n[i])-D%*%solve(step01.wl)%*%t(D)%*%solve(Vi))%*%
        (y-N*exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est)))  ## N* added
      step.wl  <- step.wl + resid %*% t(resid)         
      
      xx <- t(D) %*% solve(Vi) %*% D
      step11 <- step11+xx 
    }
    unstr.pan  <- step.pan/size  ## added
    unstr.wl   <- step.wl /size  ## added
    step11_inv <- solve(step11)  ## added

        
    step12.md <- step12.lz <- step12.kc <- step12.fg <- step12.pan <- step12.wl <- matrix(0,nrow=len,ncol=len)
    step13.md <- step13.lz <- step13.kc <- step13.fg <- step13.pan <- step13.wl <- matrix(0,nrow=len_vec,ncol=1)
    step14.md <- step14.lz <- step14.kc <- step14.fg <- step14.pan <- step14.wl <- matrix(0,nrow=len_vec,ncol=len_vec)
    p.md      <- p.lz      <- p.kc      <- p.fg      <- p.pan      <- p.wl      <- matrix(0,nrow=len_vec,ncol=size)
 
    for (i in 1:size){
      y<-as.matrix(data$response[,1][data$id==unique(data$id)[i]])  # [,1] added
      N<-y+as.matrix(data$response[data$id==unique(data$id)[i],2])  # added
      covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
      var_i=var[1:cluster$n[i],1:cluster$n[i]]
      
    ## gaussian/poisson case REMOVED 
      
        D  <- mat.prod(covariate, exp(covariate %*% beta_est)/((1+exp(covariate %*% beta_est))^2))
        D  <- diag(c(N)) %*% D #added
        Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1+exp(covariate %*% beta_est))^2)),cluster$n[i]) %*%
              var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1+exp(covariate %*% beta_est))^2)),cluster$n[i])
        Vi <- diag(sqrt(c(N))) %*% Vi %*% diag(sqrt(c(N))) #added
        Vi_inv <- solve(Vi) #added
        residual_tmp <- (y-N*exp(covariate %*% beta_est)/(1+exp(covariate %*% beta_est))) # N* added
        
        xy.lz     <- t(D) %*% Vi_inv %*% residual_tmp
        xy.lz_sq  <- xy.lz %*% t(xy.lz)
        step12.lz <- step12.lz + xy.lz_sq
        step13.lz <- step13.lz + vec(xy.lz_sq)
        p.lz[,i]  <- vec(xy.lz_sq) 
        
        I_H_inv   <- solve(cormax.ind(cluster$n[i])-D %*% step11_inv %*% t(D) %*% Vi_inv)
        xy.md     <- t(D) %*% Vi_inv %*% I_H_inv %*% residual_tmp #N* added
        xy.md_sq  <- xy.md %*% t(xy.md)
        step12.md <- step12.md + xy.md_sq
        step13.md <- step13.md + vec(xy.md_sq)
        p.md[,i]  <- vec(xy.md_sq) 
        
        xy.kc     <- t(D) %*% Vi_inv %*% mat.sqrt_pm(I_H_inv)%*% residual_tmp #N* added, function and function argument substituted
        xy.kc_sq  <- xy.kc %*% t(xy.kc)
        step12.kc <- step12.kc + xy.kc_sq
        step13.kc <- step13.kc + vec(xy.kc_sq)
        p.kc[,i]  <- vec(xy.kc_sq) 
        
        xx <- t(D) %*% Vi_inv %*% D
        Qi <- xx %*% step11_inv
        Ai <- diag((1-pmin(b,diag(Qi)))^(-0.5))
        xy.fg     <- Ai %*% t(D) %*% Vi_inv %*% residual_tmp #N* added
        xy.fg_sq  <- xy.fg%*%t(xy.fg)
        step12.fg <- step12.fg + xy.fg_sq
        step13.fg <- step13.fg + vec(xy.fg_sq)
        p.fg[,i]  <- vec(xy.fg_sq) 
        
        A.pan       <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
        diag(A.pan) <- exp(covariate %*% beta_est)/(1 + exp(covariate %*%  beta_est))^2
        xy.pan      <- t(D) %*% Vi_inv %*% sqrt(A.pan) %*% unstr.pan %*% sqrt(A.pan) %*% Vi_inv %*% D
        step12.pan  <- step12.pan + xy.pan
        step13.pan  <- step13.pan + vec(xy.pan)
        p.pan[, i]  <- vec(xy.pan)
        
        B         <- matrix(0,nrow=cluster$n[i],ncol=cluster$n[i])
        diag(B)   <- exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2
        xy.wl     <- t(D)%*%solve(Vi)%*%sqrt(B)%*%unstr.wl%*%sqrt(B)%*%solve(Vi)%*%D
        xx.wl     <- t(D)%*%solve(Vi)%*%D
        step12.wl <- step12.wl + xy.wl
        step13.wl <- step13.wl + vec(xy.wl)
        p.wl[,i]  <- vec(xy.wl) 
        
    }

    for (i in 1:size){
      
      dif.lz <- (p.lz[,i]-step13.lz/size) %*% t(p.lz[,i]-step13.lz/size)
      step14.lz<-step14.lz + dif.lz
      
      dif.md <- (p.md[,i]-step13.md/size) %*% t(p.md[,i]-step13.md/size)
      step14.md<-step14.md + dif.md

      dif.kc <- (p.kc[,i]-step13.kc/size) %*% t(p.kc[,i]-step13.kc/size)
      step14.kc<-step14.kc + dif.kc

      dif.fg <- (p.fg[,i]-step13.fg/size) %*% t(p.fg[,i]-step13.fg/size)
      step14.fg<-step14.fg + dif.fg

      dif.pan <- (p.pan[,i]-step13.pan/size) %*% t(p.pan[, i]-step13.pan/size)
      step14.pan <- step14.pan + dif.pan
      
      dif.wl <- (p.wl[,i]-step13.wl/size) %*% t(p.wl[,i]-step13.wl/size)
      step14.wl <- step14.wl + dif.wl
    }
    
    
    diag_cols          <- 1:len + (0:(len-1))*len  ## only variances of main diagonal of var(beta) needed
    kroneck_step11_inv <- kronecker(step11_inv, step11_inv) ## added to prevent multiple computation
    kroneck_cols       <- kroneck_step11_inv[,diag_cols]
    kroneck_rows       <- kroneck_step11_inv[diag_cols,]
    
    cov.beta.lz <- step11_inv %*% (step12.lz) %*% step11_inv
    cov.var.lz  <- size/(size-1) * colSums(t(kroneck_rows) * step14.lz %*% kroneck_cols) ## minor change, only diag calculated
    
    cov.beta.md <- step11_inv %*% (step12.md) %*% step11_inv
    cov.var.md  <- size/(size-1) * colSums(t(kroneck_rows) * step14.md %*% kroneck_cols) ## minor change, only diag calculated

    cov.beta.kc <- step11_inv %*% (step12.kc) %*% step11_inv
    cov.var.kc  <- size/(size-1) * colSums(t(kroneck_rows) * step14.kc %*% kroneck_cols)

    cov.beta.fg <- step11_inv %*% (step12.fg) %*% step11_inv
    cov.var.fg  <- size/(size-1) * colSums(t(kroneck_rows) * step14.fg %*% kroneck_cols)
    
    cov.beta.pan <- step11_inv %*% (step12.pan) %*% step11_inv
    cov.var.pan  <- size/(size-1) * colSums(t(kroneck_rows) * step14.pan %*% kroneck_cols)
    
    cov.beta.wl <- step11_inv %*% (step12.wl) %*% step11_inv
    cov.var.wl  <- size/(size-1) * colSums(t(kroneck_rows) * step14.wl %*% kroneck_cols)
    
    if(!verbose){
      out <- c(GEE_beta        = beta_est[2],
                   se.beta.naive   = sqrt(diag(gee.fit$naive.variance)[2]), # ["tx.var","tx.var"]),
                   se.beta.robust  = sqrt(diag(gee.fit$robust.variance)[2]),# ["tx.var","tx.var"]),
                   se.beta.lz      = sqrt(diag(cov.beta.lz)[2]),   ## should equal se.beta.robust 
                   se.beta.md      = sqrt(diag(cov.beta.md)[2]),  #  [2] returns only se of tx.var 
                   se.beta.kc      = sqrt(diag(cov.beta.kc)[2]),  #  [2] returns only se of tx.var 
                   se.beta.fg      = sqrt(diag(cov.beta.fg)[2]),  #  [2] returns only se of tx.var 
                   se.beta.pan     = sqrt(diag(cov.beta.pan)[2]), #  [2] returns only se of tx.var 
                   se.beta.wl      = sqrt(diag(cov.beta.wl)[2]),  #  [2] returns only se of tx.var 
                   var.se.beta.lz      = (cov.var.lz)[2],  #  [2] returns only se of tx.var 
                   var.se.beta.md      = (cov.var.md)[2],  #  [2] returns only se of tx.var 
                   var.se.beta.kc      = (cov.var.kc)[2],  #  [2] returns only se of tx.var 
                   var.se.beta.fg      = (cov.var.fg)[2],  #  [2] returns only se of tx.var 
                   var.se.beta.pan     = (cov.var.pan)[2],  #  [2] returns only se of tx.var 
                   var.se.beta.wl      = (cov.var.wl)[2]  #  [2] returns only se of tx.var 
      )
    } else {
      out <- list(GEE_beta        = beta_est,
                   se.beta.naive   = sqrt(diag(gee.fit$naive.variance)), 
                   se.beta.robust  = sqrt(diag(gee.fit$robust.variance)),
                   se.beta.lz      = sqrt(diag(cov.beta.lz)),   ## should equal se.beta.robust 
                   se.beta.md      = sqrt(diag(cov.beta.md)),   
                   se.beta.kc      = sqrt(diag(cov.beta.kc)),   
                   se.beta.fg      = sqrt(diag(cov.beta.fg)),   
                   se.beta.pan     = sqrt(diag(cov.beta.pan)),  
                   se.beta.wl      = sqrt(diag(cov.beta.wl)),   
                   var.se.beta.lz      = (cov.var.lz),   
                   var.se.beta.md      = (cov.var.md),   
                   var.se.beta.kc      = (cov.var.kc),   
                   var.se.beta.fg      = (cov.var.fg),   
                   var.se.beta.pan     = (cov.var.pan),  
                   var.se.beta.wl      = (cov.var.wl)  
      )                   
    }
    
    return(out) 
  }


### auxiliary function for non-symmetric matrices
mat.sqrt_pm <- function(A) {
    ei  <- eigen(A)
    d   <- ei$values
    d   <- (d + abs(d))/2
    if(isSymmetric(A)) ei_vec_inv <- t(ei$vectors)
    else               ei_vec_inv <- solve(ei$vectors)
    ans <- ei$vectors %*% diag(sqrt(d)) %*% ei_vec_inv
    return(ans)  
  }


