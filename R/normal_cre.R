normal_cre <-
function(Z_mat=Z_mat, first.order=first.order,home.field=home.field, control=control){

if(home.field&!control$OT.flag){
y_fixed_effects <- formula(~Location+0)
}else if(home.field&control$OT.flag){
y_fixed_effects <- formula(~Location+as.factor(OT)+0)   
}else if(!home.field&control$OT.flag){
y_fixed_effects <- formula(~as.factor(OT)+0) 
}else{
y_fixed_effects <- formula(~1)
}


#cutoff<-100
home_field<-home.field


#football_EM_twovar_cont<-function(football_data,home_field=FALSE, fbs_list=fbs_list_o,control=control,cutoff=100){
#Z_mat<-football_data

X<-NULL
      H.eta <- function(sigmas, cross_Z_j, Sig.mat, G.inv, nyear, 
          n_eta,sigmas2, cross_R_Z_j, Sig.mat2,R_R.inv) {
          h.eta <- G.inv
          

        h.r<-crossprod(R_Z, R_R.inv) %*% R_Z
        
          (h.eta + h.r)
      }
      ltriangle <- function(x) {
          if (!is.null(dim(x)[2])) {
              resA <- as.vector(x[lower.tri(x, diag = TRUE)])
              resA
          }
          else {
              nx <- length(x)
              d <- 0.5 * (-1 + sqrt(1 + 8 * nx))
              resB <- .symDiagonal(d)
              resB[lower.tri(resB, diag = TRUE)] <- x
              if (nx > 1) {
                  resB <- resB + t(resB) - diag(diag(resB))
              }
              as(resB, "sparseMatrix")
          }
      }

 reduce.G<-function(G) ltriangle(as.matrix(G[1:2,1:2]))
#-------------------
   update.eta <- function(X, Y, Z, cross_Z_j, Sig.mat, 
          ybetas, sigmas, G, nyear, n_eta, cons.logLik,R_X, R_Y, R_Z, 
              cross_R_Z_j, Sig.mat2, ybetas2, 
              sigmas2,R_R.inv) {
          G.chol <- chol(G)
          G.inv <- chol2inv(G.chol)
          H <- H.eta(sigmas = sigmas, cross_Z_j = cross_Z_j, Sig.mat = Sig.mat, 
              G.inv = G.inv, nyear = nyear, n_eta = n_eta, sigmas2=sigmas2,cross_R_Z_j=cross_R_Z_j,Sig.mat2=Sig.mat2,R_R.inv)
          chol.H <- chol(H)
          var.eta <- as.matrix(solve(H))
          rm(H)
          eta<-var.eta%*%t(R_Z)%*%R_R.inv%*%(R_Y-R_X%*%ybetas2)
          log.p.eta <- -(length(eta)/2) * log(2 * pi) - sum(log(diag(G.chol))) - 
              0.5 * crossprod(eta, G.inv) %*% eta
          #log.p.y <- sum(dnorm(Y, as.vector(X %*% ybetas + Z %*% 
           #   eta), as.vector(Sig.mat %*% sigmas), log = TRUE))
          log.p.r <- -(Nr/2) * log(2 * pi) + sum(log(diag(chol(R_R.inv)))) - 
            0.5 * crossprod(R_Y - R_X %*% ybetas2 - R_Z %*% eta, R_R.inv) %*% 
                (R_Y - R_X %*% ybetas2 - R_Z %*% eta)
          res <- var.eta
          attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta + 
              log.p.r - 0.5 * (2 * sum(log(diag(chol.H)))))
          attr(res, "eta") <- eta
          res
      }
      Score <- function(thetas, eta, ybetas, X, Y, Z, cross_Z_j, 
          X_j, cross_X_j, Y_j, Z_j, Sig.mat, year.count, n_ybeta, 
          nyear, n_eta, nstudent, nteacher, Kg, cons.logLik, con) {
          sigmas_sq <- pmax(thetas[(1):(nyear)], 1e-08)
          G <- thetas[seq(nyear + 1, length(thetas))]
          G <- reduce.G(G = G, nstudent = nstudent, nyear = nyear, 
              nteacher = nteacher, Kg = Kg)
          new.eta <- update.eta(X = X, Y = Y, Z = Z, 
              cross_Z_j = cross_Z_j, Sig.mat = Sig.mat, ybetas = ybetas, 
              sigmas = sqrt(sigmas_sq), G = G, nyear = nyear, n_eta = n_eta, 
              cons.logLik = cons.logLik)
          eta <- attr(new.eta, "eta")
          eta.hat <- eta
          var.eta.hat <- new.eta
          rm(new.eta)
          score.sigmas_sq <- as.vector(-1/(2 * sigmas_sq) * year.count + 
              1/(2 * sigmas_sq^2) * ((t(Sig.mat) %*% (as.vector(Y - 
                  X %*% ybetas) * as.vector(Y - X %*% ybetas - 
                  2 * Z %*% eta.hat))) + sapply(cross_Z_j, function(x) sum(diag(x %*% 
                  var.eta.hat))) + sapply(cross_Z_j, function(x) as.numeric(crossprod(eta.hat, 
                  x) %*% eta.hat))))
          temp_mat <- G
          gam_stu <- temp_mat[1, 1]
          sv_gam_stu <- 1/gam_stu
          gam_t <- list()
          sv_gam_t <- list()
          index1 <- nstudent
          for (j in 1:nyear) {
              gam_t[[j]] <- matrix(0, Kg[j], Kg[j])
              temp_mat_j <- temp_mat[(index1 + 1):(index1 + nteacher[j] * 
                  Kg[j]), (index1 + 1):(index1 + nteacher[j] * 
                  Kg[j])]
              gam_t[[j]] <- temp_mat_j[1:Kg[j], 1:Kg[j]]
              sv_gam_t[[j]] <- chol2inv(chol(gam_t[[j]]))
              index1 <- index1 + nteacher[j] * Kg[j]
          }
          rm(j)
          score_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
          gam_stu_sc <- sum(diag(score_mat)[1:nstudent])
          score.eta.stu <- drop(-0.5 * (nstudent * sv_gam_stu - 
              sv_gam_stu %*% gam_stu_sc %*% sv_gam_stu))
          score.G <- score.eta.stu * diag(nstudent)
          gam_t_sc <- list()
          index1 <- nstudent
          for (j in 1:nyear) {
              gam_t_sc[[j]] <- matrix(0, Kg[j], Kg[j])
              score_mat_j <- score_mat[(index1 + 1):(index1 + nteacher[j] * 
                  Kg[j]), (index1 + 1):(index1 + nteacher[j] * 
                  Kg[j])]
              index2 <- c(1)
              for (k in 1:nteacher[j]) {
                  gam_t_sc[[j]] <- gam_t_sc[[j]] + score_mat_j[(index2):(index2 + 
                    Kg[j] - 1), (index2):(index2 + Kg[j] - 1)]
                  index2 <- index2 + Kg[j]
              }
              index1 <- index1 + nteacher[j] * Kg[j]
              der <- -0.5 * (nteacher[j] * sv_gam_t[[j]] - sv_gam_t[[j]] %*% 
                  gam_t_sc[[j]] %*% sv_gam_t[[j]])
              if (is.numeric(drop(sv_gam_t[[j]]))) {
                  score.eta.t <- der
              }
              else {
                  score.eta.t <- 2 * der - diag(diag(der))
              }
              for (k in 1:nteacher[j]) {
                  score.G <- bdiag(score.G, score.eta.t)
              }
          }
          rm(j, k)
          -c(score.sigmas_sq, reduce.G(G = score.G, nstudent = nstudent, 
              nyear = nyear, nteacher = nteacher, Kg = Kg))
      }
    update.ybeta <- function(X, Y, Z, R_inv, eta.hat) {
        A.ybeta <- crossprod(X, R_inv) %*% X
        B.ybeta <- crossprod(X, R_inv) %*% (Y - Z %*% eta.hat)
        as.vector(solve(A.ybeta, B.ybeta))
    }
#Data format

Z_mat$home<-as.character(Z_mat$home)
Z_mat$away<-as.character(Z_mat$away)
Z_mat$year<-rep(1,dim(Z_mat)[1])
teams <- sort(unique(c(Z_mat$home,Z_mat$away)))
nteams<-length(teams)
teamsfbs<-teams
nfbs<-length(teamsfbs)


J_Y<-c(t(cbind(Z_mat$Score.For,Z_mat$Score.Against)))
#number of measured scores
Nr <- length(Z_mat$home_win)

R_RE_mat <- Matrix(0,Nr,length(teams))
J_RE_mat <- Matrix(0,2*Nr,2*length(teams))
colnames(R_RE_mat)<-teams
colnames(J_RE_mat)<-rep(teams,each=2)
for(i in 1:length(teams)){
R_RE_mat[Z_mat$home==teams[i],i]<-rep(1,length(R_RE_mat[Z_mat$home==teams[i],i]))
R_RE_mat[Z_mat$away==teams[i],i]<-rep(-1,length(R_RE_mat[Z_mat$away==teams[i],i]))
}

#offense then defense
  joffense<-c(t(cbind(Z_mat$home,Z_mat$away)))
  jdefense<-c(t(cbind(Z_mat$away,Z_mat$home)))
  J_mat<-cbind(as.numeric(J_Y),joffense,jdefense)
  J_mat<-as.data.frame(J_mat)
  templ<-rep(Z_mat$neutral.site,each=2)
  templ2<-rep("Neutral Site",length(templ))
  for(i in 1:(2*Nr)){
  if(templ[i]==0&i%%2==1)
  templ2[i]<-"Home"
  if(templ[i]==0&i%%2==0)
  templ2[i]<-"Away"
  }
  J_mat<-cbind(J_mat,templ2)
  colnames(J_mat)<-c("J_Y","offense","defense","Location")
  if(control$OT.flag){
   J_mat<-cbind(J_mat,rep(Z_mat$OT,each=2))
   colnames(J_mat)<-c("J_Y","offense","defense","Location","OT")
   }

J_mat<-as.data.frame(J_mat)
#J_mat$fcs<-as.factor(J_mat$fcs)
J_mat$J_Y<-as.numeric(J_Y)

jreo<-sparse.model.matrix(as.formula(~offense+0),data=J_mat)
jred<--1*sparse.model.matrix(as.formula(~defense+0),data=J_mat)

J_RE_mat[,seq(1,2*length(teams),by=2)]<-jreo
J_RE_mat[,seq(2,2*length(teams),by=2)]<-jred




J_X_mat <- sparse.model.matrix(y_fixed_effects, J_mat, drop.unused.levels = TRUE)

#populate the X matrix for missing data model
#R_X_mat <- sparse.model.matrix(r_fixed_effects, Z_mat, drop.unused.levels = TRUE)
#drop 0 columns from X_mat
#R_X_mat <- R_X_mat[, !(colSums(abs(R_X_mat)) == 0), drop = FALSE]
#if (rankMatrix(R_X_mat)[1] != dim(R_X_mat)[2]) {
#    cat("WARNING: Fixed-effects design matrix for missing-data model not full-rank", "\n")
#    break
#}
n_eta <- 2*length(teams)
#n_ybeta<-n_rbeta <- dim(R_X_mat)[2]
n_jbeta <- dim(J_X_mat)[2]
#X_mat<-R_X_mat

Sig.mat <- as.matrix(rep(1,Nr))
Sig.mat2 <- as.matrix(rep(1,2*Nr))
nyear<-1
#RE_mat<-R_RE_mat
#The results are loaded into the Z matricies
#Z[[]] and R_Z[[]] below -----------------------
FE.count<-0

R_X <- Matrix(J_X_mat)
R_Y <- as.numeric(as.vector(J_mat$J_Y))
#R_Y<-pmin(abs(R_Y),cutoff)*sign(R_Y)
#Z <- Matrix(new_yz)
Z<-c(NULL)
R_Z <- Matrix(J_RE_mat)
t_R_Z <- t(R_Z)


              #from here on, R means J
#X <- Matrix(R_X_mat)
Y <- as.vector(Z_mat$Score.For-Z_mat$Score.Against)
#Y<-pmin(abs(Y),cutoff)*sign(Y)

#t_Z <- t(Z)
#cross_Z <- crossprod(Z)
cross_R_Z <- crossprod(R_Z)
cross_Z_j <- list()
cross_R_Z_j <- list()
      X_j <- list(NULL)
      R_X_j <- list(NULL)
      cross_X_j <- list(NULL)
      cross_R_X_j <- list(NULL)
      Y_j <- list(NULL)
      R_Y_j <- list(NULL)
      Z_j <- list(NULL)
      R_Z_j <- list(NULL)
      for (j in 1:nyear) {
#          cross_Z_j[[j]] <- crossprod(Matrix(new_yz[Z_mat$year == 
 #             j, ]))
              cross_R_Z_j[[j]] <- crossprod(Matrix(J_RE_mat[Z_mat$year == 
              j, ]))
 #         X_j[[j]] <- X_mat[Z_mat$year == j, ]
          R_X_j[[j]] <- J_X_mat
 #         Y_j[[j]] <- as.vector(Y[Z_mat$year == j ])
          R_Y_j[[j]] <- as.vector(J_Y[Z_mat$year == j ])
 #         Z_j[[j]] <- new_yz[Z_mat$year == j, ]
          R_Z_j[[j]] <- J_RE_mat[Z_mat$year == j, ]
  #        cross_X_j[[j]] <- crossprod(X_j[[j]])
           cross_R_X_j[[j]] <- crossprod(R_X_j[[j]])
      }
#initialize parameters
eta.hat <- numeric(n_eta)
var.eta.hat <- Matrix(0, n_eta, n_eta)
G <- 100*Diagonal(n_eta)
R_R<-R_R.inv<-Diagonal(nrow(J_mat))
cons.logLik <- 0.5 * n_eta * log(2 * pi)
#Partition desigin matrices by year and
#calculate initial parameter values
# from data
         sigmas <- c(rep(0, nyear))
         sigmas2 <- c(rep(0, nyear))
ybetas<-0
ybetas2 <- update.ybeta(X=R_X, Y=R_Y, Z=R_Z, R_inv=R_R.inv, eta.hat=eta.hat)
names(ybetas2)<-colnames(J_X_mat)
#these next few lines are used to populate
#comp.list, a list of components needed for
#the E-step update (not all n_eta^2 are needed)

year.count<-Nr
iter <- control$iter.EM
r.mat <- Matrix(0, iter, length(ybetas2))
time.mat <- Matrix(0, iter, 1)
G.mat <- Matrix(0, iter, 3)
lgLik <- numeric(iter)
L1.conv <- FALSE
L2.conv <- FALSE
L1.conv.it <- 0
#Begin EM algorithm
for (it in 1:iter) {
    ptm <- proc.time()
    rm(var.eta.hat)

    new.eta <- update.eta(X = X, Y = Y, Z = Z, 
              cross_Z_j = cross_Z_j, Sig.mat = Sig.mat, ybetas = ybetas, 
              sigmas = sigmas, G = G, nyear = nyear, n_eta = n_eta, 
              cons.logLik = cons.logLik,R_X = R_X, R_Y = R_Y, R_Z = R_Z, 
              cross_R_Z_j = cross_R_Z_j, Sig.mat2 = Sig.mat2, ybetas2 = ybetas2, 
              sigmas2 = sigmas2,R_R.inv=R_R.inv)

    #save parameter values in matrix
    r.mat[it, ] <- c(ybetas2)
    lgLik[it] <- attr(new.eta, "likelihood")
    trc.y1 <- numeric(n_eta)
    trc.y2 <- Matrix(0, n_eta, n_eta)
    G.mat[it, ] <- reduce.G(G)
    eta <- as.vector(attr(new.eta, "eta"))
    var.eta <- new.eta
        eta.hat <- as.vector(eta)
        var.eta.hat <- var.eta
    rm(new.eta)
    thets1 <- c(r.mat[it - 1, ], G.mat[it - 1, ])
    thets2 <- c(r.mat[it, ], G.mat[it, ])
    # print results if verbose
    if ((control$verbose) & (it > 1)) {
        cat("\n\niter:", it, "\n")
        cat("log-likelihood:", lgLik[it], "\n")
                    cat("change in loglik:", sprintf("%.7f", lgLik[it] - 
                      lgLik[it - 1]), "\n")
        cat("n.mean:", round(ybetas2, 4), "\n")
        cat("G:", reduce.G(G),"\n")
  
    }
    if (it > 5) {
              check.lik <- abs(lgLik[it] - lgLik[it - 1])/abs(lgLik[it] + 
                  control$tol1) < control$tol1
              if (check.lik) {
                  conv <- TRUE
                  if (control$verbose) {
                    cat("\n\n Algorithm converged.\n")
                    cat("\n\niter:", it, "\n")
                    cat("log-likelihood:", sprintf("%.7f", lgLik[it]), 
                      "\n")
                    cat("change in loglik:", sprintf("%.7f", lgLik[it] - 
                      lgLik[it - 1]), "\n")
                    cat("n.mean:", round(ybetas2, 4), "\n")
                  
                             cat("G:", reduce.G(G),"\n")
                      flush.console()
                    }
                    rm(j)
                     break
                  }
                 
              }

          


          

    #Fully exponential corrections, calculated after 1st order
    #Laplace approximation converges
  

    flush.console()
    ptm.rbet <- proc.time()[3]
   # cat("starting update.ybetas", "\n")

    #cat("finished update.ybetas", proc.time()[3] - ptm.rbet, "\n")
    rm(eta)
    
  

    eblup <- as.matrix(cBind(eta.hat, sqrt(diag(var.eta.hat))))
    colnames(eblup) <- c("eblup", "std. error")
    rownames(eblup) <- rep(teams,each=2)
    rm(var.eta)

    # M-step
    #The following steps update the G
    # matrix
    
    temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
    gt1<-gt2<-matrix(0,2,2)
    for(i in 1:length(teams)){
     gt1<-gt1+temp_mat[(2*(i-1)+1):(2*i),(2*(i-1)+1):(2*i)]
    }
    gt1<-gt1/nfbs

    Gn<-kronecker(Diagonal(nfbs),symmpart(gt1))
                        sigup<-matrix(0,2,2)
for(i in 1:(nrow(J_mat)/2)){
yb<-R_Y[(2*i-1):(2*i)]
xb<-R_X[(2*i-1):(2*i),,drop=FALSE]
zb<-R_Z[(2*i-1):(2*i),]
if(home.field){
yxb<-yb-xb%*%ybetas2
}else{
yxb<-as.matrix(yb-rep(ybetas2,2))
}
sigup<- suppressWarnings(sigup+suppressMessages(tcrossprod(yxb))-yxb%*%t(zb%*%eta.hat)- (zb%*%eta.hat)%*%t(yxb)+zb%*%temp_mat%*%t(zb))
}
sigup<-symmpart(sigup/(nrow(J_mat)/2))
ybetas2 <- update.ybeta(X=R_X, Y=R_Y, Z=R_Z, R_inv=R_R.inv, eta.hat=eta.hat)

    #if(home_field) ybetas2[colnames(R_X)=="LocationNeutral Site"]<-0

    R_R<-suppressMessages(kronecker(Diagonal(nrow(J_mat)/2),sigup))
    R_R.inv<-suppressMessages(kronecker(Diagonal(nrow(J_mat)/2),solve(sigup)))
    G <- Gn
    rm(Gn)
    it.time <- (proc.time() - ptm)[3]
    time.mat[it, ] <- c(it.time)
    cat("Iteration", it, "took", it.time, "\n")
    eblup <- cBind(eta.hat, sqrt(diag(var.eta.hat)))
    colnames(eblup) <- c("eblup", "std. error")
    rownames(eblup) <- rep(teams,each=2)
}  #end EM

G.res<-as.matrix(G[1:2,1:2])
colnames(G.res)<-c("Offense","Defense")
G.res.cor<-cov2cor(G.res)
R.res<-as.matrix(R_R[1:2,1:2])
colnames(R.res)<-c("Home","Away")
R.res.cor<-cov2cor(R.res)
if(!home.field) ybetas2<-ybetas2[1]
names(ybetas2)<-colnames(J_X_mat)
  
   res<-list(n.ratings.offense=eblup[seq(1,2*nteams,by=2),1],n.ratings.defense=eblup[seq(2,2*nteams,by=2),1],p.ratings.offense=NULL,p.ratings.defense=NULL,b.ratings=NULL,n.mean=ybetas2,p.mean=NULL,b.mean=NULL,G=G.res,G.cor=G.res.cor,R=R.res,R.cor=R.res.cor,home.field=home.field)
}
