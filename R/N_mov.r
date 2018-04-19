N_mov <-
function(Z_mat=Z_mat, first.order=TRUE,home.field=home.field, control=control){
#save(list = ls(all = TRUE), file = "test21318.RData",envir=environment())    
home.field=TRUE
control$Hessian<-FALSE
#Hessian portion of the code not updated after adding REML option
#It is not really necessary since standard error for home effect is produced without it.
if(home.field){
y_fixed_effects <- formula(~1)
}else{
y_fixed_effects <- formula(~1)
}

home_field<-home.field

#------------------
#First derivative of fn.eta (score)
gr.eta <- function(eta, G.inv, n_eta,X, Y, Z, ybetas, R_inv) {
    gr.y <- crossprod(Z, R_inv) %*% (Y - X %*% ybetas - Z %*% eta)
    gr.p.eta <- -G.inv %*% eta
    -as.vector(gr.y + gr.p.eta)
}
#--------------------
#Second derivative of fn.eta (Hessian)
H.eta <- function(eta,G.inv, n_eta,Z,R_inv) {
    h.eta <- G.inv
    h.y <- crossprod(Z, R_inv) %*% Z



    symmpart(Matrix(h.eta +h.y))
}
#---------------------------------------------------
#ltriangle: extract lower triangle
#and rebuild from elements
#note ltriangle(ltriangle(x))=x
ltriangle <- function(x) {
    if (!is.null(dim(x)[2])) {
        resA <- as.vector(x[lower.tri(x, diag = TRUE)])
        resA
    } else {
        nx <- length(x)
        d <- 0.5 * (-1 + sqrt(1 + 8 * nx))
        resB <- Diagonal(d)
        resB[lower.tri(resB, diag = TRUE)] <- x
        if (nx > 1) {
            resB <- resB + t(resB) - diag(diag(resB))
        }
        resB
    }
}
#--------------
reduce.G<-function(G) ltriangle(as.matrix(G[1,1]))
#-------------------
#finds the values of the random effects that maximize the
#EM-complete data likelihood. These are the EBLUPS (must be
#corrected if non-normal data used)
update.eta <- function(eta, G, n_eta, cons.logLik,X,Y,Z,cross_Z,R_inv,ybetas) {
    G.chol <- chol(G)
    #For G, chol2inv is faster than solve
    G.inv <- symmpart(chol2inv(G.chol))
    eta <- as.vector(eta)
    H <- H.eta(eta = eta, G.inv = G.inv, n_eta = n_eta,Z=Z,R_inv=R_inv)
    chol.H <- chol(H)
    H.inv <- symmpart(chol2inv(chol.H))
    
      c.temp <- crossprod(X, R_inv ) %*% Z
      c.1 <- rbind(crossprod(X, R_inv ) %*% X, t(c.temp))
      c.2 <- rbind(c.temp, H)
      C_inv <- cbind(c.1, c.2)
      chol.C_inv <- chol(forceSymmetric(symmpart(C_inv)))
      cs <- forceSymmetric(symmpart(chol2inv(chol.C_inv)))
      C.mat <- cs[-c(1:length(ybetas)),-c(1:length(ybetas))]
      betacov <- as.matrix(cs[c(1:length(ybetas)),c(1:length(ybetas))])
      C12<-as.matrix(cs[1:length(ybetas),(length(ybetas)+1):ncol(cs)])
      if (control$REML.N) {
      var.eta <- C.mat
    } else {
      var.eta <- H.inv
    }
    rm(H) 
    eta<-H.inv%*% as.vector(crossprod(Z, R_inv)  %*%(Y-X%*%ybetas))
    

    log.p.eta <- -(length(eta)/2) * log(2 * pi) - sum(log(diag(G.chol))) - 0.5 * crossprod(eta, as(G.inv,"generalMatrix")) %*% eta
    log.p.y <- -((Ny)/2) * log(2 * pi) + sum(log(diag(chol(R_inv)))) - 0.5 * crossprod(Y - X %*% ybetas - Z %*% eta, R_inv) %*% (Y - X %*% ybetas - Z %*% eta)
    res <- var.eta
    if (control$REML.N) {
      attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta + log.p.y - 0.5 * (2 * sum(log(diag(chol.C_inv)))))
                                          #   log.p.y - 0.5 * (2 * sum(log(diag(chol.C_inv)))))
      #V<-Z%*%G%*%t(Z)+solve(R_inv)
      #Vi<-solve(V)
      #r2<-Y-X%*%ybetas                                    
      #attr(res, "likelihood") <--.5*(sum(log(eigen(V)$values))+log(det(t(X)%*%Vi%*%X))+t(r2)%*%Vi%*%r2+(Ny-ncol(X))*log(2*pi))
      
    } else {
      attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta + 
                                             log.p.y - 0.5 * (2 * sum(log(diag(chol.H)))))
    }  
    #reml sum(log(eigen(V)$values))+log(det(t(X)%*%Vi%*%X))+t(r2)%*%Vi%*%r2+(Ny-ncol(X))*log(2*pi)
    attr(res, "eta") <- eta
    attr(res, "betacov") <- betacov
    attr(res, "C12") <- C12
    res
}


#-----------------



    update.ybeta <- function(X, Y, Z, R_inv, eta.hat) {
        A.ybeta <- crossprod(X, R_inv) %*% X
        B.ybeta <- crossprod(X, R_inv) %*% (Y - Z %*% eta.hat)
        as.vector(solve(A.ybeta, B.ybeta))
    }
    
    
Score <- function(thetas) {
    ybetas <- thetas[1:n_ybeta]
    if(home.field){
    sigup<-R.tri<-R.i.parm <-  thetas[(n_ybeta+1):(n_ybeta+1)]
    }else{
    ybetas<-numeric(n_ybeta)
    sigup<-R.tri<-R.i.parm <-  thetas[(n_ybeta+1):(n_ybeta+1)] 
    }

    #R_i.parm.constrained <- R.i.parm[-LRI]
    R_i<-sigup
    R<-suppressMessages(sigup*(Diagonal(Ny)))
    R_inv<-suppressMessages((1/sigup)*(Diagonal(Ny)))
    if(home.field){
    G <- thetas[(n_ybeta+2):length(thetas)]
    }else{
     G <- thetas[(n_ybeta+2):length(thetas)]
    }
    G<-as.numeric(G)*Diagonal(length(teams))
    #update.eta returns new var.eta with eta and likelihood as attr()

    new.eta <- update.eta(eta = eta.hat,  G = G, n_eta = n_eta, cons.logLik = cons.logLik,X=X,Y=Y,Z=Z,cross_Z=cross_Z,R_inv=R_inv,ybetas=ybetas)
    eta <- attr(new.eta, "eta")
    betacov <- matrix(attr(new.eta, "betacov"),nrow=length(ybetas))
    C12<-matrix(attr(new.eta, "C12"),nrow=length(ybetas))
    eta.hat<-eta
    var.eta <- var.eta.hat <- new.eta
    eta.hat <- as.vector(eta)
    temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
   # temp_mat_R <- attr(new.eta, "h.inv") + tcrossprod(eta.hat,
    #        eta.hat)
    rm(new.eta)
    if(control$REML.N){   
        score.R <- as.numeric(-.5*(Ny/sigup-(1/sigup^2)*crossprod(Y-X%*%ybetas-Z%*%eta.hat)-(1/sigup^2)*sum(diag(crossprod(Z)%*%var.eta.hat))-(1/sigup^2)*sum(diag(crossprod(X)%*%betacov))-(1/sigup^2)*sum(diag(2*t(Z)%*%X%*%C12))))
      }else{
      score.R <- as.numeric(-.5*(Ny/sigup-(1/sigup^2)*crossprod(Y-X%*%ybetas-Z%*%eta.hat)-(1/sigup^2)*sum(diag(crossprod(Z)%*%var.eta.hat))))
    }
    
  
     gam_t_sc <- list()
        index1 <- 0
        score.G <- Matrix(0, 0, 0)
         gam_t_sc <- matrix(0, 1,1)
         index2 <- c(1)
        
                gam_t_sc <- sum(diag(temp_mat))
       
            gam_t <- G[1:1, 1:1]
            sv_gam_t <- chol2inv(chol(gam_t))
        der <- -0.5 * (nteams * sv_gam_t - sv_gam_t %*% 
                gam_t_sc %*% sv_gam_t)
            if (is.numeric(drop(sv_gam_t))) {
                score.eta.t <- der
            }else {
                score.eta.t <- 2 * der - diag(diag(der))
            }
            
           # for (k in 1:nteams) {
           #     score.G <- bdiag(score.G, score.eta.t)
           # }
        
      score.G<-ltriangle(score.eta.t)   
                                             


     -c( score.R, score.G)
  
}

    
#save(list = ls(all = TRUE), file = "test21218.RData",envir=environment())    
#-------------------------------
#Data format

Z_mat$home<-as.character(Z_mat$home)
Z_mat$away<-as.character(Z_mat$away)
teams <- sort(unique(c(Z_mat$home,Z_mat$away)))
nteams<-length(teams)
#R_mat contains data for missing data mechanism
R_mat <- Z_mat
#number of measured scores
Ny<- length(R_mat$home_win)

Y<-Z_mat$Score.For-Z_mat$Score.Against

RE_mat <- Matrix(0,Ny,length(teams))
colnames(RE_mat)<-teams

for(i in 1:length(teams)){
RE_mat[R_mat$home==teams[i],i]<-rep(1,length(RE_mat[R_mat$home==teams[i],i]))
RE_mat[R_mat$away==teams[i],i]<-rep(-1,length(RE_mat[R_mat$away==teams[i],i]))
}

n_eta <- length(teams)


J_mat<-as.data.frame(R_mat)
J_mat$J_Y<-as.numeric(Y)


X_mat <- sparse.model.matrix(y_fixed_effects, R_mat, drop.unused.levels = TRUE)
X_mat<-X_mat[,colSums(abs(X_mat))!=0,drop=FALSE]

    if (rankMatrix(X_mat,method = 'qrLINPACK')[1] != dim(X_mat)[2]) {
        stop("WARNING: Fixed-effects design matrix for scores not full-rank. May need
        to remove home field effect.")
    }
if(any(Z_mat$neutral.site==1))X_mat[Z_mat$neutral.site==1,]<-0

    
n_ybeta<-ncol(X_mat)



FE.count<-0

X <- Matrix(X_mat)
Y <- as.numeric(as.vector(J_mat$J_Y))
Z <- Matrix(RE_mat)
cross_Z<-crossprod(Z)

#initialize parameters
eta.hat <- trc.y1 <- numeric(n_eta)
var.eta.hat <- Matrix(0, n_eta, n_eta)
G <- Diagonal(n_eta)
if (control$REML.N) {
  cons.logLik <- 0.5 * (n_eta + ncol(X)) * log(2 * pi)
} else {
  cons.logLik <- 0.5 * (n_eta) * log(2 * pi)
}

#these next few lines are used to populate
#comp.list, a list of components needed for
#the E-step update (not all n_eta^2 are needed)
dummy_mat <- Diagonal(length(teams))
dummy_mat<-dummy_mat + t(Z)%*%Z
dummy_mat <- drop0(triu(dummy_mat))
comp.list <- which(abs(dummy_mat) > 1e-10, arr.ind = TRUE)
comp.list <- unique(comp.list)
iter <- control$iter.EM
y.mat <- Matrix(0, iter, n_ybeta)
time.mat <- Matrix(0, iter, 1)
G.mat <- Matrix(0, iter, 6)
R<-R_inv<-Diagonal(Ny)
ybetas <- update.ybeta(X=X, Y=Y, Z=Z, R_inv=R_inv, eta.hat=eta.hat)
if(!home.field) ybetas<-0
names(ybetas)<-colnames(X_mat)
R.mat<-Matrix(0,iter,1)
lgLik <- numeric(iter)
L1.conv <- FALSE
L2.conv <- FALSE
L1.conv.it <- 0
#save(list = ls(all = TRUE), file = "test21318.RData",envir=environment())   
#Begin EM algorithm
for (it in 1:iter) {
    ptm <- proc.time()
    #rm(var.eta.hat)

    new.eta <- update.eta(eta = eta.hat, G = G, n_eta = n_eta, cons.logLik = cons.logLik,X=X,Y=Y,Z=Z,cross_Z=cross_Z,R_inv=R_inv,ybetas=ybetas)

    #save parameter values in matrix
    y.mat[it,] <- c(ybetas)
    lgLik[it] <- attr(new.eta, "likelihood")
    trc.y1 <- numeric(n_eta)
    trc.y2 <- Matrix(0, n_eta, n_eta)
    G.mat[it, ] <- reduce.G(G)
    R.mat[it,]<-R[1,1]
    eta <- as.vector(attr(new.eta, "eta"))
    betacov <- matrix(attr(new.eta, "betacov"),nrow=length(ybetas))

    C12<-matrix(attr(new.eta, "C12"),nrow=length(ybetas))
    var.eta <- new.eta
    rm(new.eta)
    if(home.field){
    thets1 <- c(G.mat[it - 1, ],R.mat[it - 1, ],y.mat[it - 1, ])
    thets2 <- c(G.mat[it, ],R.mat[it, ],y.mat[it, ])
    }else{
        thets1 <- c(G.mat[it - 1, ],R.mat[it - 1, ],y.mat[it - 1, ])
    thets2 <- c(G.mat[it, ],R.mat[it, ],y.mat[it, ])
    }
    # print results if verbose
    if ((control$verbose) & (it > 1)) {
        cat("\n\niter:", it, "\n")
        cat("log-likelihood:", lgLik[it], "\n")
        cat("change in loglik:", sprintf("%.7f", lgLik[it] - 
                                           lgLik[it - 1]), "\n")
        cat("n.mean:", round(ybetas, 4), "\n")
        cat("G:", reduce.G(G),"\n")
        cat("R:", R.mat[it,],"\n")
    }
    #check for convergence of first order Laplace approximation
    if ((it > 5) & (L1.conv == FALSE)) {
      check.lik <- abs(lgLik[it] - lgLik[it - 1])/abs(lgLik[it] + 
                                                        control$tol1) < control$tol1
        if (check.lik)  {
            L1.conv <- TRUE
            L1.conv.it <- it
            eblup.L1 <- round(cbind(eta, sqrt(diag(var.eta))), 4)
            if (control$verbose) {
                cat("\n Algorithm has converged \n")
               
                cat("log-likelihood:", lgLik[it], "\n")
                cat("change in loglik:", sprintf("%.7f", lgLik[it] - 
                                                   lgLik[it - 1]), "\n")
        cat("G:", reduce.G(G),"\n")
        cat("R:", R.mat[it,],"\n")
            }
      
first.order.eblup<-cbind(1:length(teams),eblup[order(eblup[,1],decreasing=TRUE),])
            if(first.order) break
        }
    }
    #E-step
    #Fully exponential corrections, calculated after 1st order
    #Laplace approximation converges
            
      
      
    eta.hat <- as.vector(eta)
    flush.console()
    ptm.rbet <- proc.time()[3]


    rm(eta)
    var.eta.hat <- var.eta
    eblup <- as.matrix(cbind(eta.hat, sqrt(diag(var.eta.hat))))
    colnames(eblup) <- c("eblup", "std. error")
    rownames(eblup) <- teams
    #rm(var.eta)
      if(first.order==TRUE&L1.conv==TRUE){
        
                if (control$verbose) {
                cat("\nModel has converged!\n")
                cat("\n\niter:", it, "\n")
                cat("log-likelihood:", lgLik[it], "\n")
                cat("change in loglik:", sprintf("%.7f",lgLik[it] - 
                                                   lgLik[it - 1]), "\n")
        cat("G:", reduce.G(G),"\n")
        cat("R:", R.mat[it,],"\n")
            }
            break
                
      }
    

    # M-step
    #The following steps update the G
    # matrix

    temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
    gt1<-mean(diag(temp_mat))
    Gn<-gt1*Diagonal(length(teams))

    if(control$REML.N){   
      sigup<-(1/Ny)*(crossprod(Y-X%*%ybetas-Z%*%eta.hat)+sum(diag(crossprod(Z)%*%var.eta.hat))+sum(diag(crossprod(X)%*%betacov))+sum(diag(2*t(Z)%*%X%*%C12)))
    }else{
sigup<-(1/Ny)*(crossprod(Y-X%*%ybetas-Z%*%eta.hat)+sum(diag(crossprod(Z)%*%var.eta.hat)))
}
#sigup<-mean(diag(tcrossprod(Y-X%*%ybetas)-(Y-X%*%ybetas)%*%t(Z%*%eta.hat)-(Z%*%eta.hat)%*%t(Y-X%*%ybetas)+ Z%*%temp_mat%*%t(Z))) 

    R<-suppressMessages(sigup*(Diagonal(Ny)))
    R_inv<-suppressMessages((1/sigup)*(Diagonal(Ny)))
    if(!control$REML.N){
    ybetas<-update.ybeta(X=X, Y=Y, Z=Z, R_inv=R_inv, eta.hat=eta.hat)
    }
    if(!home.field) ybetas<-0
    G <- Gn
    
    if(control$REML.N){   
      G.chol <- chol(G)
      G.inv <- chol2inv(G.chol)
      R.inv.Z <- R_inv %*% Z
      V.1 <- symmpart(chol2inv(chol(G.inv + t(Z) %*% R.inv.Z)))
      tX.Rinv.Z <- t(X) %*% R.inv.Z
      tX.Rinv.X <- t(X) %*% R_inv %*% X
      
      ybetas <-
        as.vector(chol2inv(chol(forceSymmetric(
          symmpart(tX.Rinv.X -
                     tX.Rinv.Z %*% V.1 %*% t(tX.Rinv.Z))
        ))) %*% (t(X) %*% R_inv -
                   tX.Rinv.Z %*% V.1 %*% t(R.inv.Z)) %*% Y)
    } 
    
    rm(Gn)
    it.time <- (proc.time() - ptm)[3]
    time.mat[it, ] <- c(it.time)
    cat("Iteration", it, "took", it.time, "\n")
    eblup <- cbind(eta.hat, sqrt(diag(var.eta.hat)))

    colnames(eblup) <- c("eblup", "std. error")
    rownames(eblup) <- teams
    #write.csv(eblup,file="NB_mvglmmRank_ratings.csv")
}  #end EM

if(home.field){
thetas <- c(ybetas, as.numeric(sigup), reduce.G(G))
names(thetas)<-c("Mean","R[1,1]","G[1,1]")
}else{
thetas <- c( as.numeric(sigup), reduce.G(G))
names(thetas)<-c("R[1,1]","G[1,1]")
}

Hessian<-NULL
gradient<-Score(thetas)
if(control$Hessian){
cat("\nCalculating Hessian with a central difference approximation...\n")
flush.console()
Hessian <- symmpart(jacobian(Score, thetas, method="simple"))
rownames(Hessian)<-colnames(Hessian)<-names(thetas)
#std_errors <- c(sqrt(diag(solve(Hessian))))
if(!all(eigen(Hessian)$values>0)) cat("\nWarning: Hessian not positive-definite\n")
}

#if(control$REML.N){
  c.temp <- crossprod(X, R_inv) %*% Z
  c.1 <- rbind(crossprod(X, R_inv) %*% X, t(c.temp))
  G.inv <- chol2inv(chol(G))
  c.2 <- rbind(c.temp,H.eta(eta = eta, G.inv = G.inv, n_eta = n_eta,Z=Z,R_inv=R_inv))
              
  C_inv <- cbind(c.1, c.2)
  C <- solve(C_inv)
  eblup_stderror <- sqrt(diag(C)[-c(1:ncol(X))])
  ybetas_stderror <- sqrt(diag(C)[1:ncol(X)])
  ybetas_asycov<-C[1:ncol(X),1:ncol(X)]
  ybetas_eblup_asycov<-C
  rm(C, C_inv, c.2, c.1, c.temp)
  eblup <- as.matrix(cbind(eta.hat, eblup_stderror))
  eblup <- as.data.frame(eblup)
  eblup <- as.data.frame(cbind(colnames(Z), eblup)) 
#}else{
  
#  ybetas_asycov<-ybetas_eblup_asycov<-ybetas_stderror<-c("This feature only enabled for REML.N=TRUE")
#}

FE.X<-cbind(X,Z)

colnames(FE.X)<-c("Intercept",teams)
#FE.fit<-lm.fit(as.matrix(FE.X),Y)
FE.beta<-as.vector(ginv(as.matrix(t(FE.X)%*%FE.X))%*%t(FE.X)%*%Y)
names(FE.beta)<-c("Intercept",teams)
FE.G<-ginv(as.matrix(t(FE.X)%*%FE.X))
FE.H<-FE.G%*%t(FE.X)%*%FE.X
FE.c<-numeric(nrow(FE.H))
FE.c[1]<-1
FE.mean.estimable=all(zapsmall(t(FE.c)%*%FE.H)==t(FE.c))
FE.pred<-as.vector(FE.X%*%FE.beta)
FE.resid<-as.vector(Y-FE.pred)
FE.var<-as.numeric((t(FE.resid)%*%FE.resid)/(Ny-length(FE.beta)))
FE.beta.cov<-FE.var*ginv(as.matrix(t(FE.X)%*%FE.X))
colnames(FE.beta.cov)<-c("Intercept",teams)
FE.hfe.se<-sqrt(diag(FE.beta.cov))
fixed.effect.model.output<-list(X=FE.X,beta=FE.beta,is.mean.estimable=FE.mean.estimable,pred=FE.pred,resid=FE.resid,sigma.sq=as.numeric(FE.var),beta.covariance=FE.beta.cov)

G.res<-as.matrix(G[1,1])
colnames(G.res)<-c("Margin of Victory")
G.res.cor<-1
R.res<-as.matrix(R[1,1])
colnames(R.res)<-c("Variance")
R.res.cor<-cov2cor(R.res)
names(ybetas)<-"LocationHome"
N.output<-list(Z=Z,Y=Y,X=X,G=G,R=R,eta=eta.hat,var.eta=var.eta,ybetas_eblup_asycov=ybetas_eblup_asycov, ybetas_asycov=ybetas_asycov,ybetas_stderror=ybetas_stderror,gradient=gradient )
   res<-list(n.ratings.mov=eblup[,1:2],n.ratings.offense=NULL,n.ratings.defense=NULL,p.ratings.offense=NULL,p.ratings.defense=NULL,b.ratings=NULL,n.mean=ybetas,p.mean=NULL,b.mean=NULL,G=G.res,G.cor=G.res.cor,R=R.res,R.cor=R.res.cor,home.field=home.field,actual=Y,pred=X%*%ybetas+Z%*%eta.hat,Hessian=Hessian,parameters=thetas,N.output=N.output,fixed.effect.model.output=fixed.effect.model.output)
}
