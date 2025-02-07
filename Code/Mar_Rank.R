library(mvtnorm)
library(pracma)
library(tmvtnorm)
#Estimating rank of MAR Data matrix

#Simulates FA/PCA matrix
gen.mat <- function(N, p, r, load.mean =rep(6,r), opt.errdist= "norm", opt.mod = "FA", diag.variability = FALSE, opt.missing = FALSE, missing.prob=0 ){

	lam.mat   <- matrix(rnorm(N*r,0,1), N, r)
	if( diag.variability == FALSE){

		id.mat    <- diag(runif(p, 1, 1),p, p)

	}else{

		id.mat    <- diag(runif(p, 1, 5),p, p)
	}

	if(opt.errdist == "norm"){
		noise.mat <- rmvnorm(N, runif(p,0,1) , 5*id.mat)
	}else if( opt.errdist =="t"){
		noise.mat <- rmvt(N, df=10, sigma= diag(p)*10, delta=rep(0,p) , type="shifted")
	}else if( opt.errdist =="truncnorm"){
	  noise.mat <- rtmvnorm(N, mean=rep(0, p), sigma=id.mat, lower=rep(-1, p), upper=rep(1, p), algorithm = "gibbs")
	}else if( opt.errdist =="uniform"){
	  noise.mat <- matrix(runif(N*p, min=-1, max=1), N, p)
	}else if(opt.errdist == "mixnorm"){
	  # Parameters for the two normals
	  mu1 <- rep(0, p)  # mean vector for the first normal
	  sd1 <- id.mat     # covariance matrix for the first normal
	  mu2 <- rep(5, p)  # mean vector for the second normal
	  sd2 <- diag(runif(p, 1, 5),p, p)     # covariance matrix for the second normal
	  lambda <- 0.5     # mixing proportion

	  # For each observation, decide which normal to sample from
	  which_norm <- rbinom(N, 1, lambda)

	  noise.mat <- matrix(0, N, p)
	  for(i in 1:N){
	    if(which_norm[i] == 1) {
	      noise.mat[i,] <- rmvnorm(1, mu1, sd1)
	    } else {
	      noise.mat[i,] <- rmvnorm(1, mu2, sd2)
	    }
	  }
	}

	if( opt.mod == "FA"){

		if(r==1 ){
			load.mat <- matrix(rnorm(p*r,0, 1),r, p)
			load.mat <- load.mat/sqrt( sum(load.mat^2))

			load.mat <- load.mean*load.mat
			x        <- lam.mat %*%load.mat + noise.mat%*%id.mat

		}else{

			load.mat       <- matrix(0,r, p)
			diag(load.mat) <- rep(1, r, 1)
			r.dum          <- r + 1
			load.mat[,r.dum : p] <- 1

			orth.mat <- orth(load.mat)
			col.ind  <- rep(NA, p)
			col.ind[1:r] <- seq(1,r,1)
			p.rem        <- p -r
			r.dum        <- r+ 1
			col.ind[r.dum:p] <- sample( seq(1,r,1),p.rem, replace=TRUE)
		    load.mat         <- orth.mat[,col.ind]
			load.mat <- diag(load.mean)%*%load.mat
			x        <- lam.mat %*%load.mat + noise.mat%*%id.mat
		}

	}else{
		if(r==1 ){
			load.mat       <- matrix(0,r, p)
			diag(load.mat) <- rep(1,r)
			col.ind  <- rep(NA, p)
			col.ind[1:r] <- seq(1,r,1)
			p.rem        <- p -r
			r.dum        <- r+ 1
			col.ind[r.dum:p] <- sample( seq(1,r,1),p.rem, replace=TRUE)
		    load.mat         <- load.mat[,col.ind]

			load.mat <- load.mean*load.mat
			x        <- lam.mat %*%load.mat + noise.mat%*%id.mat

		}else{

			load.mat       <- matrix(0,r, p)
			diag(load.mat) <- rep(1,r)

			col.ind  <- rep(NA, p)
			col.ind[1:r] <- seq(1,r,1)
			p.rem        <- p -r
			r.dum        <- r+ 1
			col.ind[r.dum:p] <- sample( seq(1,r,1),p.rem, replace=TRUE)
		    load.mat         <- load.mat[,col.ind]

			load.mat <- diag(load.mean)%*%load.mat
			x        <- lam.mat %*%load.mat + noise.mat%*%id.mat
		}

	}

	if(opt.missing){
	 censor.mat <- matrix(1, N, p)
	 censor.mat <- matrix(rbinom(N*p, 1, 1-missing.prob), N, p)
	 x          <- x*censor.mat
	}

	list(x=x, lam.mat = lam.mat, mu = lam.mat %*%load.mat, Sigma = id.mat, signal = load.mean)
}


sol <- function(v,t,w, gam){
	z <- -1/v + gam*sum(t*w/(1+t*v))
}

upper_edge <- function(t,w, gam){
	max.phi <- max(t)
	low     <- -1/max.phi
	z <- optimize(sol, c(low,0),maximum=FALSE,t=t,w=w, gam=gam, tol= 1e-16)$objective
z
}

dpa <- function(mat, opt.scale= FALSE, scale.fac=1.05, K_m = 10, missing.prob= 0){

  N    <- nrow(mat)
  p    <- ncol(mat)
  gam  <- p/N

  out <- 0

  #############################
  # Sequential Test Procedure
  #############################

  if(opt.scale){
    mat <- scale(mat, center=TRUE,scale=TRUE)
  }

  obj <- svd(mat)

  for(k in 1:K_m){
    if (k == 0){
      mat <- mat
    }else{
      mat <- mat - obj$d[k]*(obj$u[,k])%*%(t(obj$v[,k]))
    }
    X   <- mat
    S   <- t(X)%*%X/N
    D   <- diag(S)

    lambda <- svd(S)$d[1]

    num_clus <- k
    t        <- D[1:num_clus]
    w       <- rep(1,num_clus)/num_clus
    z.alpha	<- upper_edge(t,w, gam)

    if( lambda <= scale.fac*z.alpha){
      out = k
      break
    }

  }

  out
}

dpa_mar <- function(mat, opt.scale= FALSE, scale.fac=1.05, K_m = p, missing.prob= 0){
  N    <- nrow(mat)
  p    <- ncol(mat)
  gam  <- p/N
  out <- 0
  
  # Replace NAs with 0 for missing values
  mat[is.na(mat)] <- 0
  
  if(opt.scale){
    mat <- scale(mat, center=TRUE,scale=TRUE)
  }
  
  obj <- svd(mat)
  
  for(k in 0:K_m){
    if (k == 0){
      curr_mat <- mat
    }else{
      curr_mat <- mat - obj$d[k]*(obj$u[,k])%*%(t(obj$v[,k]))
    }
    X <- curr_mat
    S <- t(X)%*%X/N
    D <- diag(S)
    
    lambda <- svd(S)$d[1]
    
    num_clus <- k
    t <- D[1:num_clus]
    w <- rep(1,num_clus)/num_clus
    
    # Handle the case when missing.prob is 0
    if(missing.prob > 0) {
      z.alpha <- upper_edge(t,w, gam)*(1/(missing.prob))
      z.alpha <- z.alpha + gam/(missing.prob*(1-missing.prob))
    } else {
      z.alpha <- upper_edge(t,w, gam)
    }
    
    if(lambda <= scale.fac*z.alpha){
      out = k
      break
    }
  }
  out
}

dpa_mar_3 <- function(mat, opt.scale= FALSE, scale.fac=1.05, K_m = p, missing.prob= 0){

  N    <- nrow(mat)
  p    <- ncol(mat)
  gam  <- p/N

  out <- 0

  #############################
  # Sequential Test Procedure
  #############################
  update.dpa = FALSE
  res.dpa<- 0
  var <- rep(0, N)
  gam <- p/N

  if(opt.scale){
    mat <- scale(mat, center=TRUE, scale=TRUE)
  }

  obj <- svd(mat)

  for(k in 0:K_m){
    if (k == 0){
      mat <- mat
    }else{
      mat <- mat - obj$d[k]*(obj$u[,k])%*%(t(obj$v[,k]))
    }
    X   <- mat
    S   <- t(X)%*%X/N
    D   <- diag(S)

    lambda <- svd(S)$d[1]

    num_clus <- k
    t        <- D[1:num_clus]
    w       <- rep(1,num_clus)/num_clus
    z.alpha	<- upper_edge(t,w, gam)*(1/(missing.prob))
    z.alpha <- z.alpha + gam * ((missing.prob)/(1- missing.prob))

    if( lambda <= scale.fac*z.alpha){
      out = k
      break
    }

  }

  out
}

dpa_mar_2 <- function(mat, opt.scale= FALSE, scale.fac=1.05, K_m = p, missing.prob= 0){

mat <- cov(mat)
N    <- nrow(mat)
p    <- ncol(mat)
gam  <- p/N

out <- 0

#############################
# Sequential Test Procedure
#############################
update.dpa = FALSE
res.dpa<- 0
var <- rep(0, N)
gam <- p/N

if(opt.scale){
	mat <- scale(mat, center=TRUE, scale=TRUE)
}


obj <- svd(mat)

for(k in 0:K_m){
	if (k == 0){
		mat <- mat
	}else{
		mat <- mat - obj$d[k]*(obj$u[,k])%*%(t(obj$v[,k]))
	}
	X   <- mat
	S   <- t(X)%*%X/N
	D   <- diag(S)

	lambda <- svd(S)$d[1]

	num_clus <- k
	t        <- D[1:num_clus]
	w       <- rep(1,num_clus)/num_clus
	z.alpha	<- upper_edge(t,w, gam)*(1/(missing.prob))

	#print( paste(k,z.alpha, lambda, upper_edge(t,w, gam), upper_edge(t,w, gam)*1/missing.prob, sep=","))
	if( lambda <= scale.fac*z.alpha){
		out = k
		break
	}

}
out

}

dpa_mar_1 <- function(mat, opt.scale= FALSE, scale.fac=1.05, K_m = p, missing.prob= 0){

mat <- cor(mat)
N    <- nrow(mat)
p    <- ncol(mat)
gam  <- p/N

out <- 0

#############################
# Sequential Test Procedure
#############################
update.dpa = FALSE
res.dpa<- 0
var <- rep(0, N)
gam <- p/N

if(opt.scale){
	mat <- scale(mat, center=TRUE, scale=TRUE)
}


obj <- svd(mat)

for(k in 0:K_m){
	if (k == 0){
		mat <- mat
	}else{
		mat <- mat - obj$d[k]*(obj$u[,k])%*%(t(obj$v[,k]))
	}
	X   <- mat
	S   <- t(X)%*%X/N
	D   <- diag(S)

	lambda <- svd(S)$d[1]

	num_clus <- k
	t        <- D[1:num_clus]
	w       <- rep(1,num_clus)/num_clus
	z.alpha	<- upper_edge(t,w, gam)*(1/(missing.prob))

	#print( paste(k,z.alpha, lambda, upper_edge(t,w, gam), upper_edge(t,w, gam)*1/missing.prob, sep=","))
	if( lambda <= scale.fac*z.alpha){
		out = k
		break
	}

}
out
}


result <- matrix(0,10,4)
for (i in 1:10){
#Simulates FA/PCA matrix
N <- 500
p <- 100
mat.s <-  gen.mat(N, p, r=6, load.mean =rep(12,6), opt.errdist= "norm", opt.mod = "FA", diag.variability = FALSE, opt.missing = TRUE, missing.prob=0.4)
mat <- mat.s$x
phat <- length(which(mat==0))/(N*p)

#out_dpa <- dpa(mat, opt.scale = TRUE)
result[i,1] <- dpa_mar(mat, opt.scale = TRUE, missing.prob = phat)
result[i,2] <- dpa_mar_1(mat, opt.scale = TRUE, missing.prob = phat)
result[i,3] <- dpa(mat, opt.scale = TRUE, missing.prob = phat)
result[i,4] <- dpa_mar_3(mat, opt.scale = TRUE, missing.prob = phat)

}

result	

N <- 500
p <- 100
mat.s <-  gen.mat(N, p, r=3, load.mean =rep(6,3), opt.errdist= "norm", opt.mod = "FA", diag.variability = FALSE, opt.missing = TRUE, missing.prob=0.5)
mat <- mat.s$x
phat <- length(which(mat==0))/(N*p)

dpa_mar(mat, opt.scale = TRUE, missing.prob = phat)
dpa_mar_3(mat, opt.scale = TRUE, missing.prob = phat)

out_2 <- dpa(mat, opt.scale = TRUE, missing.prob = 1 - phat)
out_2
out_0 <- dpa(mat, opt.scale = TRUE, missing.prob = 1)
out_0