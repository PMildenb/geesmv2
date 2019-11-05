require(matrixcalc)
#' GEE.var.md_bin.aggr
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): 
#' a formula expression as for other regression models to be fitted, 
#' of the form c(successes,failures) ~ predictors. The details of formula specification can be seen in glm() and gee().
#' @param id a vector which identifies the clusters. The length of id should be the same as the total number of observations. 
#' Data is assumed to be sorted so that observations on a cluster are contiguous rows for all entities in the formula.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) 
#' containing the variables in the model. If not found in data, the variables are taken from environment(formula), 
#' typically the environment from which glm is called.
#' @param corstr a character string specifying working correlation structure: 
#' "independence", "AR-M","exchangeable", "unstructured"  are possible.
#'
#' @return
#' @export
#'
#' @examples
GEE.var.md_bin.aggr  <-
  function(formula,id,data,corstr="independence"){
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
    
    ### na.action: only na.omit is used for gee;
    init <- model.frame(formula, data)
    init$num <- 1:length(init[,2]) #### 1 changed to 2
    if(any(is.na(init))){
      index <- na.omit(init)$num
      data <- data[index,]
      ### Get the design matrix;
      m <- model.frame(formula, data)
      mt <- attr(m, "terms") 
      data$response <- model.response(m, "numeric")
      mat <- as.data.frame(model.matrix(formula, m))
    }else{
      ### Get the design matrix;
      m <- model.frame(formula, data)
      mt <- attr(m, "terms") 
      data$response <- model.response(m, "numeric")
      mat <- as.data.frame(model.matrix(formula, m))
    }
    
    ### Fit the GEE model to get the estimate of parameters \hat{\beta};
    gee.fit <- gee(formula,data=data,id=id,family=family,corstr=corstr)
    beta_est <- gee.fit$coefficient
    alpha <- gee.fit$working.correlation[1,2]
    len <- length(beta_est)
    len_vec <- len^2
    
    ### Estimate the robust variance for \hat{\beta}
    data$id <- gee.fit$id
    cluster<-cluster.size(data$id)
    ncluster<-max(cluster$n)
    size<-cluster$m
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
    
    family <- "binomial"
    
    
    cov.beta<-unstr<-matrix(0,nrow=len,ncol=len)
    step11<-matrix(0, nrow=len, ncol=len)
    for (i in 1:size){
      y<-as.matrix(data$response[data$id==unique(data$id)[i],1])  # ,1 added
      N<-y+as.matrix(data$response[data$id==unique(data$id)[i],2])  # added
      covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
      var_i=var[1:cluster$n[i],1:cluster$n[i]]
      if (family=="gaussian"){ 
        xx<-t(covariate)%*%solve(var_i)%*%covariate
        step11<-step11+xx  
      }else if (family=="poisson") {
        D<-mat.prod(covariate, exp(covariate%*%beta_est))
        Vi <- diag(sqrt(c(exp(covariate%*%beta_est))),cluster$n[i])%*%
              var_i%*%diag(sqrt(c(exp(covariate%*%beta_est))),cluster$n[i])
        xx<-t(D)%*%solve(Vi)%*%D
        step11<-step11+xx
      }else if (family=="binomial"){
        D <- mat.prod(covariate, exp(covariate%*%beta_est)/((1+exp(covariate%*%beta_est))^2))
        D <- diag(c(N))%*%D #added
        Vi <- diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),cluster$n[i])%*%
              var_i%*%diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),cluster$n[i])
        Vi <- diag(sqrt(c(N)))%*%Vi%*%diag(sqrt(c(N))) #added
        xx<-t(D)%*%solve(Vi)%*%D
        step11<-step11+xx 
      }
    }
    step12<-matrix(0,nrow=len,ncol=len)
    step13<-matrix(0,nrow=len_vec,ncol=1)
    step14<-matrix(0,nrow=len_vec,ncol=len_vec)
    p<-matrix(0,nrow=len_vec,ncol=size)
    for (i in 1:size){
      y<-as.matrix(data$response[,1][data$id==unique(data$id)[i]])  # [,1] added
      N<-y+as.matrix(data$response[data$id==unique(data$id)[i],2])  # added
      covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
      var_i=var[1:cluster$n[i],1:cluster$n[i]]
      if (family=="gaussian"){ 
        xy<-t(covariate)%*%solve(var_i)%*%solve(cormax.ind(cluster$n[i])-covariate%*%solve(step11)%*%t(covariate)%*%solve(var_i))%*%(y-covariate%*%beta_est)
        step12<-step12+xy%*%t(xy)
        step13<-step13+vec(xy%*%t(xy))
        p[,i]<-vec(xy%*%t(xy))
      }else if (family=="poisson") {
        D<-mat.prod(covariate, exp(covariate%*%beta_est))
        Vi <- diag(sqrt(c(exp(covariate%*%beta_est))),cluster$n[i])%*%var_i%*%diag(sqrt(c(exp(covariate%*%beta_est))),cluster$n[i])
        xy<-t(D)%*%solve(Vi)%*%solve(cormax.ind(cluster$n[i])-D%*%solve(step11)%*%t(D)%*%solve(Vi))%*%(y-exp(covariate%*%beta_est))
        step12<-step12+xy%*%t(xy)
        step13<-step13+vec(xy%*%t(xy))
        p[,i]<-vec(xy%*%t(xy))
      }else if (family=="binomial"){
        D<-mat.prod(covariate, exp(covariate%*%beta_est)/((1+exp(covariate%*%beta_est))^2))
        D<-diag(c(N))%*%D #added
        Vi <- diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),cluster$n[i])%*%var_i%*%diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),cluster$n[i])
        Vi <- diag(sqrt(c(N)))%*%Vi%*%diag(sqrt(c(N))) #added
        xy<-t(D)%*%solve(Vi)%*%solve(cormax.ind(cluster$n[i])-D%*%solve(step11)%*%t(D)%*%solve(Vi))%*%(y-N*exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))) #N* added
        step12<-step12+xy%*%t(xy)
        step13<-step13+vec(xy%*%t(xy))
        p[,i]<-vec(xy%*%t(xy)) 
      }    
    }
    for (i in 1:size){
      dif<-(p[,i]-step13/size)%*%t(p[,i]-step13/size)
      step14<-step14+dif
    }
    cov.beta<-solve(step11)%*%(step12)%*%solve(step11)
    cov.var<-size/(size-1)*kronecker(solve(step11), solve(step11))%*%step14%*%kronecker(solve(step11), solve(step11))
    return(list(cov.beta=diag(cov.beta), cov.var=cov.var,D=D,Vi=Vi,beta_est=beta_est,xy=xy))
  }