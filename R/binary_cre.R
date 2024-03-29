binary_cre <-
function(Z_mat=Z_mat, first.order=first.order,home.field, control=control){
y_fixed_effects <- formula(~Location)
r_fixed_effects <- formula(~neutral.site)
home_field<-home.field
#------------------
#First derivative of fn.eta (score)
gr.eta <- function(eta, R_X, R_Y, R_Z, rbetas, G.inv, n_eta) {
    gr.r <- colSums(R_Z * as.vector((-1)^(1 - R_Y) * rder1((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta))))
    gr.p.eta <- -G.inv %*% eta
    -as.vector(gr.r + gr.p.eta)
}
#--------------------
#Second derivative of fn.eta (Hessian)
H.eta <- function(eta, rbetas, R_X, R_Y, R_Z,G.inv, n_eta) {
    h.eta <- G.inv
    temp <- -rder2((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta))
    h.r <- Matrix(crossprod(R_Z, temp * R_Z))

    Matrix(h.eta + h.r)
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
#These functions are derivatives of the standard normal CDF that are used in the missing data mechanism
rder0 <- function(q) {
    q <- as.vector(q)
    pnorm(q, log.p = TRUE)
}
rder1 <- function(q) {
    q <- as.vector(q)
    dnorm(q)/pmax(pnorm(q), 1e-15)
}
rder2 <- function(q) {
    q <- as.vector(q)
    (-(q/sqrt(2 * pi) * exp(-q^2/2)) * pnorm(q) - dnorm(q)^2)/pmax(pnorm(q)^2, 1e-15)
}
rder3 <- function(q) {
    q <- as.vector(q)
    ((q^2 - 1)/sqrt(2 * pi) * exp(-q^2/2) * pnorm(q)^2 + 3 * (q/sqrt(2 * pi) * exp(-q^2/2)) * pnorm(q) * dnorm(q) + 2 * dnorm(q)^3)/pmax(pnorm(q)^3, 1e-15)
}
rder4 <- function(q) {
    q <- as.vector(q)
    ((-q * (q^2 - 3)/sqrt(2 * pi) * exp(-q^2/2)) * pnorm(q)^3 - 4 * ((q^2 - 1)/sqrt(2 * pi) * exp(-q^2/2)) * pnorm(q)^2 * dnorm(q) - 12 * (q/sqrt(2 * pi) * exp(-q^2/2)) * pnorm(q) * dnorm(q)^2 - 3 * (q/sqrt(2 * pi) * exp(-q^2/2))^2 * pnorm(q)^2 - 6 * 
        dnorm(q)^4)/pmax(pnorm(q)^4, 1e-15)
}
reduce.G<-function(G) G[1,1]
#-------------------
#finds the values of the random effects that maximize the
#EM-complete data likelihood. These are the EBLUPS (must be
#corrected if non-normal data used)
update.eta <- function(eta, R_X, R_Y, R_Z, rbetas, G, nyear, n_eta, cons.logLik) {
    G.chol <- chol(G)
    #For G, chol2inv is faster than solve
    G.inv <- chol2inv(G.chol)
    eta <- as.vector(eta)
    H <- H.eta(eta = eta, rbetas = rbetas,  R_X = R_X, R_Y = R_Y, R_Z = R_Z, G.inv = G.inv, n_eta = n_eta)
    #For H, solve is much faster than chol2inv
    var.eta <- as.matrix(solve(H))
    product <- matrix(rep(1, n_eta), n_eta, 1)
    gr <- rep(1, n_eta)
    gr.loop.it<-1
    while (as.numeric(gr %*% product) > 1e-08&gr.loop.it<10) {
        gr <- gr.eta(eta = eta, R_X = R_X, R_Y = R_Y, R_Z = R_Z, rbetas = rbetas,  G.inv = G.inv, n_eta = n_eta)
        product <- var.eta %*% gr
        eta <- eta - product
        H <- H.eta(eta = eta, rbetas = rbetas, R_X = R_X, R_Y = R_Y, R_Z = R_Z, G.inv = G.inv,  n_eta = n_eta)
        var.eta <- as.matrix(solve(H))
        gr.loop.it<-gr.loop.it+1
    }
    rm(product, gr)
    chol.H <- chol(H)
    log.p.eta <- -(length(eta)/2) * log(2 * pi) - sum(log(diag(G.chol))) - 0.5 * crossprod(eta, as(G.inv,"generalMatrix")) %*% eta
    log.p.r <- sum(pnorm(as.vector((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta)), log.p = TRUE))
    res <- var.eta
    attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta  + log.p.r - 0.5 * (2 * sum(log(diag(chol.H)))))
    attr(res, "eta") <- eta
    res
}


#-----------------
Sc.rbetas.f <- function(eta, rbetas, R_X, R_Y, R_Z, n_eta, Svar, Svar2) {
    E.rbetas <- ((-1)^(1 - R_Y) * rder1((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta)))
    R2 <- rder2((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta))
    R3 <- (-1)^(1 - R_Y) * rder3((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta))
    temp.C <- (R2 * R_Z)
    block <- (Svar %*% t(temp.C)) * R3
    trc.rbetas <- rowSums(Svar2 * R3) + colSums(rowSums(Svar2) * block)
    score.rbetas <- as.vector(t(R_X) %*% (E.rbetas + 0.5 * trc.rbetas))
    rm(temp.C, block)
    score.rbetas
}
Sc.rbetas.f.first.order <- function(eta, rbetas,  R_X, R_Y, R_Z) {
      E.rbetas <- ((-1)^(1 - R_Y) * rder1((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta)))
    score.rbetas <- as.vector(t(R_X) %*% E.rbetas )
    score.rbetas
}

#update rbetas
#update for fixed effects in missing data mechansim
#requires var.eta, not var.eta.hat
update.rbetas <- function(eta, var.eta, rbetas, R_X, R_Y, R_Z, n_eta) {
    Svar <- as.matrix(R_Z %*% var.eta)
    Svar2 <- Svar * R_Z
    product <- matrix(rep(1, length(rbetas)), length(rbetas), 1)
    scrbetas <- rep(1, length(rbetas))
    rbeta.it.count<-1
    while (as.numeric(scrbetas %*% product) > 1e-08& rbeta.it.count<10) {
        Hrbetas <- as.matrix(jacobian(Sc.rbetas.f, rbetas, method = "simple", eta = eta, R_X = R_X, R_Y = R_Y, R_Z = R_Z, n_eta = n_eta, Svar = Svar, Svar2 = Svar2))
        scrbetas <- Sc.rbetas.f(eta, rbetas, R_X, R_Y, R_Z, n_eta, Svar = Svar, Svar2 = Svar2)
        product <- c(solve(Hrbetas, scrbetas))
        rbetas <- rbetas - product
        rbeta.it.count<- rbeta.it.count + 1
    }
    rm(Svar, Svar2, scrbetas, product)
    rbetas
}
update.rbetas.first.order <- function(eta, rbetas, R_X, R_Y, R_Z, n_eta) {
    
    product <- matrix(rep(1, length(rbetas)), length(rbetas), 1)
    scrbetas <- rep(1, length(rbetas))
    while (as.numeric(scrbetas %*% product) > 1e-08) {
        Hrbetas <- as.matrix(jacobian(Sc.rbetas.f.first.order, rbetas, method = "simple", eta = eta, R_X = R_X, R_Y = R_Y, R_Z = R_Z))
        scrbetas <- Sc.rbetas.f.first.order(eta, rbetas, R_X, R_Y, R_Z)
        product <- c(solve(Hrbetas, scrbetas))
        rbetas <- rbetas - product
    }
    rm( scrbetas, product)
    rbetas
    }
    
    
    Score <- function(thetas) {
   
    if(home.field){
      rbetas <- thetas[1]
    if(dim(R_X)[2]==2) rbetas<-c(rbetas,0)
    }else{
    rbetas <- rep(0,ncol(R_X))
    }
    if(home.field){
    G <- thetas[2]
    }else{
     G <- thetas[1]
    }
    G<-G*Diagonal(length(teams))
    #update.eta returns new var.eta with eta and likelihood as attr()
    new.eta <- update.eta(eta = eta.hat, R_X = R_X, R_Y = R_Y, R_Z = R_Z,rbetas = rbetas, G = G, n_eta = n_eta, cons.logLik = cons.logLik)
    eta <- attr(new.eta, "eta")
    eta.hat<-eta
    var.eta <- var.eta.hat <- new.eta
    eta.hat <- as.vector(eta)
    temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
   # temp_mat_R <- attr(new.eta, "h.inv") + tcrossprod(eta.hat,
    #        eta.hat)
    rm(new.eta)

   if(home.field) score.r <- Sc.rbetas.f.first.order(eta = eta, rbetas = rbetas, R_X = R_X, R_Y = R_Y, R_Z = R_Z)[1]          
    
     gam_t_sc <- list()
        index1 <- 0
        score.G <- Matrix(0, 0, 0)
         gam_t_sc <- matrix(0, 1,1)
         index2 <- c(1)
         for (k in 1:nrow(G)) {
                gam_t_sc <- gam_t_sc + temp_mat[(index2):(index2 + 
                  0), (index2):(index2 + 0)]
                index2 <- index2 + 1
            }
            gam_t <- G[1, 1]
            sv_gam_t <- chol2inv(chol(gam_t))
        der <- -0.5 * (nrow(G) * sv_gam_t - sv_gam_t %*% 
                gam_t_sc %*% sv_gam_t)
            if (is.numeric(drop(sv_gam_t))) {
                score.eta.t <- der
            }
            else {
                score.eta.t <- 2 * der - diag(diag(der))
            }
            
           # for (k in 1:nteams) {
           #     score.G <- bdiag(score.G, score.eta.t)
           # }
        
      score.G<-ltriangle(score.eta.t)   
                                             

     if(home.field){      
    -c(score.r,  score.G)
    }else{
     -c( score.G)
    }
}

#-------------------------------
#Data format

#Z_mat$home<-as.character(Z_mat$home)
#Z_mat$away<-as.character(Z_mat$away)
teams <- sort(unique(c(Z_mat$home,Z_mat$away)))
#R_mat contains data for missing data mechanism
R_mat <- Z_mat
#number of measured scores
Nr <- length(R_mat$home_win)

R_RE_mat <- Matrix(0,Nr,length(teams))
colnames(R_RE_mat)<-teams
for(i in 1:length(teams)){
R_RE_mat[R_mat$home==teams[i],i]<-rep(1,length(R_RE_mat[R_mat$home==teams[i],i]))
R_RE_mat[R_mat$away==teams[i],i]<-rep(-1,length(R_RE_mat[R_mat$away==teams[i],i]))
}
#populate the X matrix for missing data model
R_X_mat <- sparse.model.matrix(r_fixed_effects, R_mat, drop.unused.levels = TRUE)
#drop 0 columns from X_mat
R_X_mat <- R_X_mat[, !(colSums(abs(R_X_mat)) == 0), drop = FALSE]
    if (rankMatrix(R_X_mat,method = 'qrLINPACK')[1] != dim(R_X_mat)[2]) {
        stop("WARNING: Fixed-effects design matrix not full-rank")
    }
n_eta <- length(teams)
n_rbeta <- dim(R_X_mat)[2]

#The results are loaded into the Z matricies
#Z[[]] and R_Z[[]] below -----------------------
FE.count<-0
R_X <- Matrix(R_X_mat)

R_Y <- as.vector(R_mat$home_win)
R_Z <- Matrix(R_RE_mat)
t_R_Z <- t(R_Z)
cross_R_Z <- crossprod(R_Z)

if(any(R_Y==2)){
tie.indx<-which(R_Y==2)
X.tie<-R_X[tie.indx,,drop=FALSE]
Z.tie<-R_Z[tie.indx,,drop=FALSE]
R_Y<-R_Y[-tie.indx]
R_X<-R_X[-tie.indx,,drop=FALSE]
R_Z<-R_Z[-tie.indx,,drop=FALSE]
for(i in 1:length(tie.indx)){
R_Y<-c(R_Y,1,0)
R_X<-rbind(R_X,X.tie[i,,drop=FALSE],X.tie[i,,drop=FALSE])
R_Z<-rbind(R_Z,Z.tie[i,,drop=FALSE],Z.tie[i,,drop=FALSE])
}
t_R_Z <- t(R_Z)
cross_R_Z <- crossprod(R_Z)
Nr <- length(R_Y)
}

#initialize parameters
eta.hat <- trc.y1 <- numeric(n_eta)
var.eta.hat <- Matrix(0, n_eta, n_eta)
G <- Diagonal(n_eta)   
cons.logLik <- 0.5 * n_eta * log(2 * pi)
#Partition desigin matrices by year and
#calculate initial parameter values
# from data

rbetas <- rep(0,ncol(R_X_mat))
#these next few lines are used to populate
#comp.list, a list of components needed for
#the E-step update (not all n_eta^2 are needed)
dummy_mat <- Diagonal(n_eta)
comp.list <- which(abs(dummy_mat) > 1e-10, arr.ind = TRUE)
comp.list <- unique(comp.list)
comp.list[comp.list[,1]<=comp.list[,2],]
iter <- control$iter.EM
r.mat <- Matrix(0, iter, n_rbeta)
time.mat <- Matrix(0, iter, 1)
G.mat <- Matrix(0, iter, 1)
lgLik <- numeric(iter)
L1.conv <- FALSE
L2.conv <- FALSE
L1.conv.it <- 0
#save(list = ls(all = TRUE), file = "test31518.RData",envir=environment())   
#Begin EM algorithm
for (it in 1:iter) {
    ptm <- proc.time()
    #rm(var.eta.hat)

    new.eta <- update.eta(eta = eta.hat, R_X = R_X, R_Y = R_Y, R_Z = R_Z,rbetas = rbetas, G = G, n_eta = n_eta, cons.logLik = cons.logLik)

    #save parameter values in matrix
    r.mat[it, ] <- c(rbetas)
    lgLik[it] <- attr(new.eta, "likelihood")
    trc.y1 <- numeric(n_eta)
    trc.y2 <- Matrix(0, n_eta, n_eta)
    G.mat[it, ] <- reduce.G(G)
    eta <- as.vector(attr(new.eta, "eta"))
    var.eta <- new.eta
    rm(new.eta)
    thets1 <- c(G.mat[it - 1, ])
    thets2 <- c(G.mat[it, ])
    # print results if verbose
    if ((control$verbose) & (it > 1)) {
        cat("\n\niter:", it, "\n")
        cat("log-likelihood:", lgLik[it], "\n")
        cat("max % change in parm.", 100 * max(round(abs((thets2 - thets1)/thets1 * as.numeric(abs(thets1) > .00001)), 5)), "%\n")
        cat("b.mean:", round(rbetas, 4), "\n")
        cat("G:", reduce.G(G),"\n")
    }
    #check for convergence of first order Laplace approximation
    if ((it > 5) & (L1.conv == FALSE)) {
        check.parm1 <- max(abs(thets2 - thets1) * as.numeric(abs(thets1) <= .00001)) < control$tol1
        check.parm2 <- max(abs((thets2 - thets1)/thets1 * as.numeric(abs(thets1) > .00001))) < control$tol1
        if (check.parm1 & check.parm2) {
            L1.conv <- TRUE
            L1.conv.it <- it
            eblup.L1 <- round(cbind(eta, sqrt(diag(var.eta))), 4)
            if (control$verbose) {
                cat("\n First order Laplace has converged \n")
                cat("log-likelihood:", lgLik[it], "\n")
                cat("max % change in parm.", 100 * round(max(abs(thets2 - thets1)/abs(thets1) * as.numeric(abs(thets1) > .00001)), 5), "%\n")
                        cat("b.mean:", round(rbetas, 4), "\n")
        cat("G:", reduce.G(G),"\n")
            }
             
#first.order.eblup<-cbind(1:120,eblup[order(eblup[,1],decreasing=TRUE),])
            if(first.order) break
        }
    }
    #E-step
    #Fully exponential corrections, calculated after 1st order
    #Laplace approximation converges
          if (L1.conv == TRUE) {
              Svar <- as.matrix(R_Z %*% var.eta)
             temp.trc.D <- as.vector(rder4((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta)))
          temp.trc.C <- (-1)^(1 - R_Y) * rder3((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta))
          trc.y1 <- rep(0,n_eta)
          check.parmFE1 <- max(abs(thets2 - thets1) * as.numeric(abs(thets1) <= 0.0001)) < control$tolFE
          check.parmFE2 <- max(abs((thets2 - thets1)/thets1 * as.numeric(abs(thets1) > 0.0001))) < control$tolFE
          #calculate corrections to var.eta.hat after only calculating corrections to eta.hat for a while

              #compute the trace corrections for var.eta

              #The function trace.calc2 is used to calculate the trace corrections
              Trace.calc2 <- function(comp.k,comp.l,dsig.dc.indx5, indx.5) {
                  if (comp.k !=comp.l) {
                    if(comp.k!=indx.5) comp2<-comp.k
                    if(comp.l!=indx.5) comp2<-comp.l
                    q<-var.eta%*%(dsig.dc.indx5%*%var.eta[,comp2])
                    d2.sig<- t_R_Z%*%(R_Z*(as.vector(R_Z%*%q)*temp.trc.C+ Svar[,indx.5]*Svar[,comp2]*temp.trc.D))
                    dsig.dc.comp2<- t_R_Z%*%(R_Z*Svar[,comp2]*temp.trc.C)
                    res<-sum(diag(var.eta%*%d2.sig)) + sum((t(dsig.dc.indx5)%*%var.eta)*(var.eta%*%dsig.dc.comp2))
                  }else{

                    q<-var.eta%*%(dsig.dc.indx5%*%var.eta[,indx.5])
                    d2.sig<- t_R_Z%*%(R_Z*(as.vector(R_Z%*%q)*temp.trc.C+ Svar[,indx.5]*Svar[,indx.5]*temp.trc.D))
                    res<-sum(diag(var.eta%*%d2.sig)) + sum((t(dsig.dc.indx5)%*%var.eta)*(var.eta%*%dsig.dc.indx5))

                  }
                   res

              }
              temp.comp <- comp.list
              length.comp <- dim(comp.list)[1]
              var.corrections <- matrix(0, 0, 3)
              for (indx.5 in n_eta:1) {
                  ptm.cor <- proc.time()[3]
                  guide1 <- unique(which(temp.comp == indx.5, arr.ind = TRUE)[, 1])
                  guide2 <- comp.list[guide1, , drop = FALSE]
                  if (length(guide1) > 0) {
                  if(indx.5==n_eta)cat("Calculating FE corrections for random effects vector...\n")
                    dsig.dc.indx5<- t_R_Z%*%(R_Z*Svar[,indx.5]*temp.trc.C)
                    trc.y1[indx.5]<- sum(diag(var.eta%*%dsig.dc.indx5))
              if (check.parmFE1 & check.parmFE2  | FE.count > 0) {
              if(indx.5==n_eta)cat("Calculating FE corrections for random effects covariance matrix...\n")
              FE.count <- FE.count + 1
                    res.temp <- cbind(guide2, mapply(Trace.calc2, comp.k = guide2[, 1], comp.l = guide2[, 2], MoreArgs = list(dsig.dc.indx5, indx.5 = indx.5)))
                    temp.comp <- temp.comp[-(guide1), , drop = FALSE]
                    var.corrections <- rbind(var.corrections, res.temp)
                    }
   
                  }
                #  cat(1-dim(temp.comp)[1]/length.comp,'%,',indx.5,',',proc.time()[3]-ptm.cor,'\n')
              }
              rm(indx.5)
                        if (check.parmFE1 & check.parmFE2  | FE.count > 0) {
       
              trc.y2 <- sparseMatrix(i = var.corrections[, 1], j = var.corrections[, 2], x = var.corrections[, 3], dims = c(n_eta, n_eta))
              trc.y2 <- trc.y2 + t(trc.y2) - diag(diag(trc.y2))
              }
              rm(Svar,dsig.dc.indx5)
          
      }
    
    
    eta.hat <- as.vector(eta + 0.5 * trc.y1)
    flush.console()
    ptm.rbet <- proc.time()[3]

if(first.order==FALSE&L1.conv==TRUE){
rbetasn <- update.rbetas(eta = eta, var.eta = var.eta, rbetas = rbetas, R_X=R_X, R_Y=R_Y, R_Z=R_Z, n_eta=n_eta)
}else{
rbetasn <- update.rbetas.first.order(eta = eta, rbetas = rbetas, R_X=R_X, R_Y=R_Y, R_Z=R_Z, n_eta=n_eta)
}
    rm(eta)
    var.eta <- var.eta + 0.5 * trc.y2
    rm(trc.y2)
    var.eta.hat <- var.eta
    eblup <- as.matrix(cbind(eta.hat, sqrt(diag(var.eta.hat))))
    colnames(eblup) <- c("eblup", "std. error")
    rownames(eblup) <- teams
    rm(var.eta)
              if(first.order==TRUE&L1.conv==TRUE){
        
                if (control$verbose) {
                cat("\nFirst-order Model has converged!\n")
                cat("\n\niter:", it, "\n")
                cat("log-likelihood:", lgLik[it], "\n")
                cat("max % change in parm.", 100 * round(max(abs((thets2 - thets1)/thets1) * as.numeric(abs(thets1) > 0.01)), 5), "%\n")
                        cat("b.mean:", round(rbetas, 4), "\n")
        cat("G:", reduce.G(G),"\n")
            }
            break
                
      }
    
    
    # check convergence of Fully exponential model
    if ((L1.conv == TRUE) & (it > L1.conv.it)) {
        check.parm1 <- max((thets2 - thets1) * as.numeric(abs(thets1) <= .00001)) < control$tol2
        check.parm2 <- max(abs((thets2 - thets1)/thets1) * as.numeric(abs(thets1) > .00001)) <control$tol2
        if (check.parm1 & check.parm2) {
            L2.conv <- TRUE
            if (control$verbose) {
                cat("\nFully Exponential Model has converged!\n")
                cat("\n\niter:", it, "\n")
                cat("log-likelihood:", lgLik[it], "\n")
                cat("max % change in parm.", 100 * round(max(abs((thets2 - thets1)/thets1) * as.numeric(abs(thets1) > .00001)), 5), "%\n")
                        cat("b.mean:", round(rbetas, 4), "\n")
        cat("G:", reduce.G(G),"\n")
            }
            break
        }
    }
    # M-step
    #The following steps update the G
    # matrix
    
    temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
    Gn<-Diagonal(n_eta)*mean(diag(temp_mat))

    rbetas <- rbetasn
    if(home_field){
     rbetas[colnames(R_X)=="neutral.site"]<-0
     }else{
      rbetas<-numeric(length(rbetas))
     }
    G <- Gn
    rm(Gn)
    it.time <- (proc.time() - ptm)[3]
    time.mat[it, ] <- c(it.time)
    cat("Iteration", it, "took", it.time, "\n")
    eblup <- cbind(eta.hat, sqrt(diag(var.eta.hat)))
    colnames(eblup) <- c("eblup", "std. error")
    rownames(eblup) <- teams
}  #end EM


if(home.field){
thetas <- c( rbetas[1], G[1,1])
names(thetas)<-c("Binary mean","G[1,1]")
}else{
thetas <- c(G[1,1])
names(thetas)<-c("G[1,1]")
}

Hessian<-NULL
if(control$Hessian){
cat("\nCalculating Hessian with a central difference approximation...\n")
flush.console()
Hessian <- symmpart(jacobian(Score, thetas, method="simple"))
rownames(Hessian)<-colnames(Hessian)<-names(thetas)
#std_errors <- c(sqrt(diag(solve(Hessian))))
if(inherits(try(chol(Hessian),silent=TRUE), "try-error")) cat("\nWarning: Hessian not positive-definite\n")
}

model.matrices=list(X=R_X,Y=R_Y,Z=R_Z,G=G,eta=eta.hat,beta=rbetas)

            G.res<-as.matrix(G[1,1])
colnames(G.res)<-c("Rating")
G.res.cor<-cov2cor(G.res)
 res<-list(n.ratings.mov=NULL,n.ratings.offense=NULL,n.ratings.defense=NULL,p.ratings.offense=NULL,p.ratings.defense=NULL,b.ratings=eblup[,1],n.mean=NULL,p.mean=NULL,b.mean=rbetas[1],G=G.res,G.cor=G.res.cor,R=NULL, R.cor=NULL,home.field=home.field,Hessian=Hessian,parameters=thetas,model.matrices=model.matrices)
}
