# simulation data for multivariate joint model linear
sim_mjm_linear = function(I,  obstime = 1:10, miss = FALSE, miss.rate = 0.1, opt = c("none","interaction")){
    
    # I : number of subjects
    # J : number of visits
    # obstime: observation times
    # miss: whether introduce missing (missing complete at random) in longitudinal data. Different from drop-out
    # miss.rate: missing rate.
    
    J = length(obstime)
    N = I*J
    
    ## longitudinal submodel
    beta0 = c(1.5,2,0.5)
    beta1 = c(2,-1,1)
    betat = c(1.5, -1, 0.6) 
    b.var = c(1,1.5,2)
    e.var = c(1,1,1)
    rho = c(-0.2,0.1,-0.3)
    b.Sigma = diag(b.var)
    b.Sigma[1,2] = b.Sigma[2,1] = sqrt(b.var[1]*b.var[2])*rho[1]
    b.Sigma[1,3] = b.Sigma[3,1] = sqrt(b.var[1]*b.var[3])*rho[2]
    b.Sigma[2,3] = b.Sigma[3,2] = sqrt(b.var[2]*b.var[3])*rho[3]
    
    # sample covariate
    X = rep(rnorm(I, 3, 1), each=J)
    # sample random effect
    ranef = mvrnorm(I, c(0,0,0), b.Sigma)
    id = rep(1:I,each=J)
    ranef = ranef[id,]
    # construct longitudinal submodel
    eta.long = matrix(0, nrow=N, ncol=3)
    for(i in 1:3)
        eta.long[,i] = beta0[i] + beta1[i]*X + ranef[,i]
    
    
    ## survival submodel
    #interaction
    if(opt=="interaction"){
        gamma = c(-4,-2,4)
        alpha = c(0.2, -0.2, 0.4)
        x1 = rbinom(I, size = 1, prob=.5)
        x2 = rnorm(I)
        x3 = x1*x2
        W = cbind(x1,x2,x3)
        eta.surv = W%*%gamma + c(alpha%*%t(eta.long[!duplicated(id),]))
    }
        
    else{
        gamma = c(-4,-2)
        alpha = c(0.2, -0.2, 0.4)
        x1 = rbinom(I, size = 1, prob=.5)
        x2 = rnorm(I)
        W = cbind(x1,x2)
        eta.surv = W%*%gamma + c(alpha%*%t(eta.long[!duplicated(id),]))
    }
    
    ## simulate survival time
    scale = exp(-7)
    S = runif(I)
    Ti = rep(NA, I)
    alpha.beta = alpha%*%betat
    f = function(tau){
        h = function(t) {
            scale *exp(eta.surv[i] + c(alpha.beta)*t)
        }
        S[i] - exp(-stats::integrate(h, 0, tau)$value)
    }
    f = Vectorize(f)
    
    for(i in 1:I){
        Ti[i] = uniroot(f, c(0, 100))$root
    }
    
    ## simulate true survival probability
    pre.surtime = function(tau){
        h = function(t) {
            scale *exp(eta.surv[i] + c(alpha.beta)*t)
        }
        
        exp(-stats::integrate(h, 0, tau)$value)
    }
    
    true.prob = matrix(NA, nrow=I, ncol=length(obstime[-1]))
    for(i in 1:I){
        ith = 0
        for(tau in obstime[-1]){
            ith = ith + 1
            true.prob[i, ith] = pre.surtime(tau)
        }
    }
    
    colnames(true.prob) = as.character(obstime[-1])
    
    #--------------------------------
    # simulate censor time
    C = runif(I,min=obstime[3], max=obstime[5]+20)
    time <- pmin(Ti, C) #observed time is min of censored and true
    event <- ifelse(time==Ti, 1, 0) #0: censored ; 1: event; 
    
    # prepare data
    visit = rep(1:J, I)
    obstime = rep(obstime, I) 
    erro = mvrnorm(N, c(0,0,0), diag(e.var))
    Y = matrix(0, nrow=N, ncol=3)
    for(i in 1:3)
        Y[,i] = eta.long[,i] + betat[i]*obstime  + erro[,i]
    
    
    long.all = data.frame(id=id, visit=visit, time = rep(time, each=J), event = rep(event, each=J),
                          Y1=Y[,1], Y2=Y[,2], Y3=Y[,3],obstime=obstime, X=X, ranef=ranef, W  = W[rep(1:I, each=J)], erro=I(erro))
    
    long = long.all
    
    # introduce missing complete at random
    if(miss){
        miss.index = sample(which(long$obstime>obstime[2]), miss.rate*N)
        long = long[!c(1:N)%in%miss.index, ]
    }
    
    surv = data.frame(id = c(1:I),time=time, event=event, W = W, true.prob=I(true.prob))
    
    # remove observations after event or censoring
    long = long[long$obstime<long$time, ]
    
    return(list(long=long, surv=surv, long.all=long.all))
}





# univariate FPCA via principal analysis by conditional estimation(PACE)
uPACE = function(testData, domain, predData=NULL, nbasis = 10, pve = 0.9, npc = NULL){
    
    tmp = funData(domain, testData)
    if(is.null(predData)){
        tmp2 = NULL
    }else{
        tmp2 = funData(domain, predData)
    }
    res = PACE(tmp, tmp2, pve=pve, npc= npc, nbasis=nbasis)
    return(res)
} 

# multivariate FPCA based on results from uPACE
mFPCA = function(Xi, phi, p , L, I=I){
    
    # eigenanalysis on matrix M
    M = t(Xi)%*%Xi/(I-1)
    eigen.M = eigen(M)
    values = eigen.M$values
    pve = cumsum(values)/sum(values)
    Cms = eigen.M$vectors
    index = unlist(lapply(1:length(L), function(x) rep(x, L[x])))
    
    # MFPCA score
    rho = mfpca.score(Xi, Cms)
    
    # MFPCA eigenfunction
    psis = NULL
    for(j in 1:p){
        psi = NULL
        for(m in 1:dim(Cms)[2]){
            psi = cbind(psi, phi[[j]]%*%Cms[which(index==j),m])
        }
        psis[[j]] = psi
    }
    
    out = list(eigenvalue = values, Cms = Cms, pve = pve, index=index, rho = rho, psis=psis)
    
    return(out)
}

# mfpc score calculation
mfpca.score = function(predXi, Cms){
    rho = matrix(NA, nrow = nrow(predXi), ncol=dim(Cms)[2])
    for(i in 1:nrow(predXi)){
        for(m in 1:dim(Cms)[2]){
            rho[i,m] = predXi[i,]%*%Cms[,m]
        }
    }
    return(rho)
}


# mfpc trajectories prediction
mfpca.pred = function(score, meanf, psi, n.rho=NULL){
    p = length(psi)
    n = nrow(score)
    
    if(is.null(n.rho)){
        n.rho = ncol(score)
    }
    
    pred = array(NA, c(n, length(meanf[[1]]), p))
    for(m in 1:p){
        pred[,,m] = matrix(meanf[[m]], nrow=n, ncol =length(meanf[[m]]), byrow = T ) + score[,1:n.rho]%*%t(psi[[m]][, 1:n.rho])
    }
    
    out = pred
    return(out)
}


#risk predict using predictSurvProb instead of survfit to handle RSF. (REQUIRES PEC)
cond.prob.pec = function(model, newdata, Tstart, Tpred){
    risk.Tstart = as.numeric(predictSurvProb(model,newdata=newdata,times=Tstart))
    risk.Tpred = as.numeric(predictSurvProb(model,newdata=newdata,times=Tpred))
    return(risk.Tpred/risk.Tstart)
}

