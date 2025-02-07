# Required libraries
library(mvtnorm)
library(pracma)
library(tmvtnorm)

#Simulates FA/PCA matrix
gen.mat <- function(N, p, r, load.mean =rep(6,r), opt.errdist= "norm", opt.mod = "FA", diag.variability = FALSE, opt.missing = FALSE, missing.prob=0 ){
  # ... [previous gen.mat implementation] ...
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
    mu1 <- rep(0, p)
    sd1 <- id.mat
    mu2 <- rep(5, p)
    sd2 <- diag(runif(p, 1, 5),p, p)
    lambda <- 0.5
    
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

# Original DPA-MAR-3 function
dpa_mar_3 <- function(mat, opt.scale=FALSE, scale.fac=1.05, K_m=ncol(mat), missing.prob=0) {
  N <- nrow(mat)
  p <- ncol(mat)
  gam <- p/N
  out <- 0
  
  if(opt.scale) {
    mat <- scale(mat, center=TRUE, scale=TRUE)
  }
  
  obj <- svd(mat)
  curr_mat <- mat
  
  for(k in 0:K_m) {
    if(k > 0) {
      curr_mat <- curr_mat - obj$d[k]*(obj$u[,k])%*%(t(obj$v[,k]))
    }
    
    S <- t(curr_mat)%*%curr_mat/N
    D <- diag(S)
    lambda <- svd(S)$d[1]
    
    num_clus <- k
    t <- D[1:max(1,num_clus)]
    w <- rep(1,max(1,num_clus))/max(1,num_clus)
    z.alpha <- upper_edge(t,w,gam)*(1/missing.prob)
    z.alpha <- z.alpha + gam*((missing.prob)/(1-missing.prob))
    
    if(lambda <= scale.fac*z.alpha) {
      out <- k
      break
    }
  }
  return(out)
}

# New version using limiting distribution theory
# Key differences from dpa_mar_3:
# 1. Threshold calculation: Uses (1/p0)*(1 + sqrt(gamma))^2 - ((1-p0)/p0) instead of upper_edge function
# 2. Missing value handling: Explicitly calculates observed probability p0 from data
# 3. Scaling: Implements custom scaling for missing values instead of using scale()
dpa_mar_limit <- function(mat, opt.scale=FALSE, scale.fac=1.05, K_m=ncol(mat), missing.prob=0) {
  N <- nrow(mat)
  p <- ncol(mat)
  gamma <- p/N
  out <- 0
  
  # Different from dpa_mar_3: Calculate observed probability from data
  mat[mat == 0] <- NA
  p0 <- 1 - mean(is.na(mat))
  
  # Different from dpa_mar_3: Custom scaling for missing values
  if(opt.scale) {
    mat_scaled <- mat
    for(j in 1:ncol(mat)) {
      col_mean <- mean(mat[,j], na.rm=TRUE)
      col_sd <- sd(mat[,j], na.rm=TRUE)
      mat_scaled[,j] <- (mat[,j] - col_mean) / col_sd
    }
    mat <- mat_scaled
  }
  
  mat[is.na(mat)] <- 0
  curr_mat <- mat
  obj <- svd(curr_mat)
  
  for(k in 0:K_m) {
    if(k > 0) {
      curr_mat <- curr_mat - obj$d[k]*(obj$u[,k])%*%(t(obj$v[,k]))
    }
    
    S <- t(curr_mat)%*%curr_mat/N
    lambda <- eigen(S)$values[1]
    
    # Different from dpa_mar_3: Simpler threshold based on limiting distribution
    z.alpha <- (1/p0)*(1 + sqrt(gamma))^2 - ((1-p0)/p0)
    
    if(lambda <= scale.fac*z.alpha) {
      out <- k
      break
    }
  }
  return(out)
}

# Helper functions for dpa_mar_3
# sol: Calculates the objective function for finding the upper edge
# Input:
#   v: variable to optimize
#   t: eigenvalues
#   w: weights
#   gam: gamma = p/N ratio
sol <- function(v,t,w,gam) {
  -1/v + gam*sum(t*w/(1+t*v))
}

# upper_edge: Finds the upper edge of the limiting spectral distribution
# Used by dpa_mar_3 for threshold calculation
# Input:
#   t: eigenvalues
#   w: weights
#   gam: gamma = p/N ratio
upper_edge <- function(t,w,gam) {
  max.phi <- max(t)
  low <- -1/max.phi
  optimize(sol, c(low,0), maximum=FALSE, t=t, w=w, gam=gam, tol=1e-16)$objective
}

# Run comparison
N <- 500
p <- 100
missing_probs <- c(0.2, 0.4, 0.6)
ranks <- c(3, 5, 8)
signals <- c(3, 6, 12, 24)  
num_trials <- 5

results <- data.frame()
detailed_results <- list()

for(mp in missing_probs) {
  for(r in ranks) {
    for(s in signals) {
      cat("\nTesting: missing_prob =", mp, "rank =", r, "signal =", s, "\n")
      
      trial_results <- list()
      
      for(i in 1:num_trials) {
        # Generate data
        mat.s <- gen.mat(N, p, r=r, 
                        load.mean=rep(s,r), 
                        opt.errdist="norm", 
                        opt.mod="FA", 
                        diag.variability=FALSE, 
                        opt.missing=TRUE, 
                        missing.prob=mp)
        mat <- mat.s$x
        
        # Store singular values
        singular_values <- svd(mat)$d
        singular_values <- singular_values[1:min(20, length(singular_values))]
        
        # Run both methods
        rank_mar3 <- dpa_mar_3(mat, opt.scale=TRUE, missing.prob=mp)
        rank_limit <- dpa_mar_limit(mat, opt.scale=TRUE, missing.prob=mp)
        
        # Store detailed results for this trial
        trial_results[[i]] <- list(
          singular_values = singular_values,
          matrix = mat,
          rank_mar3 = rank_mar3,
          rank_limit = rank_limit
        )
        
        # Store summary results
        results <- rbind(results, data.frame(
          missing_prob = mp,
          true_rank = r,
          signal = s,
          trial = i,
          rank_mar3 = rank_mar3,
          rank_limit = rank_limit
        ))
      }
      
      # Store detailed results for this configuration
      detailed_results[[paste("mp", mp, "r", r, "s", s, sep="_")]] <- trial_results
    }
  }
}

# Create timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")

# Print summary
cat("\nResults Summary:\n")
summary_df <- data.frame()

for(mp in missing_probs) {
  for(r in ranks) {
    for(s in signals) {
      subset <- results[results$missing_prob == mp & 
                       results$true_rank == r & 
                       results$signal == s,]
      
      cat(sprintf("\nmissing_prob = %.1f, true_rank = %d, signal = %.1f\n", 
                  mp, r, s))
      cat("MAR3 average rank:", mean(subset$rank_mar3), "\n")
      cat("Limit average rank:", mean(subset$rank_limit), "\n")
      cat("MAR3 accuracy:", mean(subset$rank_mar3 == r), "\n")
      cat("Limit accuracy:", mean(subset$rank_limit == r), "\n")
      
      # Add to summary dataframe
      summary_df <- rbind(summary_df, data.frame(
        missing_prob = mp,
        true_rank = r,
        signal = s,
        mar3_avg_rank = mean(subset$rank_mar3),
        limit_avg_rank = mean(subset$rank_limit),
        mar3_accuracy = mean(subset$rank_mar3 == r),
        limit_accuracy = mean(subset$rank_limit == r)
      ))
    }
  }
}

# Create experiment_data directory if it doesn't exist
dir.create("../experiment_data", showWarnings = FALSE)

# Save all results
save(results, detailed_results, summary_df, 
     file=sprintf("../experiment_data/dpa_mar_comparison_%s.RData", timestamp))

# Save summary as CSV
write.csv(summary_df, 
          file=sprintf("../experiment_data/dpa_mar_comparison_summary_%s.csv", timestamp),
          row.names=FALSE)

# Print overall performance by signal strength
cat("\n\nOverall Performance by Signal Strength:\n")
for(s in signals) {
  subset <- summary_df[summary_df$signal == s,]
  cat(sprintf("\nSignal Strength = %.1f:\n", s))
  cat("MAR3 average accuracy:", mean(subset$mar3_accuracy), "\n")
  cat("Limit average accuracy:", mean(subset$limit_accuracy), "\n")
}
